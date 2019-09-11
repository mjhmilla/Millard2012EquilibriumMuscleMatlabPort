% -------------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  %
% See http:%opensim.stanford.edu and the NOTICE file for more information.  %
% OpenSim is developed at Stanford University and supported by the US        %
% National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    %
% through the Warrior Web program.                                           %
%                                                                            %
% Copyright (c) 2005-2012 Stanford University and the Authors                %
% Author(s): Matthew Millard                                                 %
%                                                                            %
% Licensed under the Apache License, Version 2.0 (the 'License'); you may    %
% not use this file except in compliance with the License. You may obtain a  %
% copy of the License at http:%www.apache.org/licenses/LICENSE-2.0.         %
%                                                                            %
% Unless required by applicable law or agreed to in writing, softNware        %
% distributed under the License is distributed on an 'AS IS' BASIS,          %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   %
% See the License for the specific language governing permissions and        %
% limitations under the License.                                             %
% -------------------------------------------------------------------------- %
%
% Derivative work
% Authors(s): Millard
% Updates   : ported to code to Matlab
%
% If you use this code in your work please cite this paper
%
%  Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). 
%    Flexing computational muscle: modeling and simulation of 
%    musculotendon dynamics. Journal of biomechanical engineering, 
%    135(2), 021005.
%%
function normMuscleCurves = createDefaultNormalizedMuscleCurves(...
                muscleName,tendonStrainAtOneNormForceInput,flag_updateCurves, flag_plotCurves)
%%
% This function constructs the default normalized muscle fiber and tendon
% characteristic curves that are used in the Millard2012EquilibriumMuscle
% model in OpenSim and puts them into a struct.
%
% @param muscleName: the name of the muscle. This string is used to give
%                    each of the curves a specific name
%
% @param flag_plotCurves: 0: nothing will happen
%                         1: each curves value, 1st, and 2nd derivative
%                         will be plotted over the full Beziercurve domain 
%                         and a little beyond so that you can see
%                         extrapolation. If the curve has an integral then
%                         its value too will be plotted
%                         
% @return normMuscleCurves: a struct containing the following fields
%
%          .activeForceLengthCurve
%          .fiberForceLengthCurve
%          .fiberForceVelocityCurve
%          .tendonForceLengthCurve
%
%          .activeForceLengthCurveHACK
%          .fiberForceVelocityInverseCurveHACK
%                       
%
%%
normMuscleCurves = [];

if(exist('normMuscleCurves.mat','file') == 2 && flag_updateCurves == 0)
    tmp = load('normMuscleCurves.mat');
    normMuscleCurves = tmp.normMuscleCurves;
else    
    %%
    %Active force length curve
    %%
    lce0 = 0.47-0.0259;
    lce1 = 0.73;
    lce2 = 1.0;
    lce3 = 1.8123;
    minActiveForceLengthValue = 0;
    curviness = 1.0;
    computeIntegral = 0;
    plateauSlope = 0.8616;

    activeForceLengthCurve = createFiberActiveForceLengthCurve( lce0,...
                                                                lce1, ...
                                                                lce2, ...
                                                                lce3, ...
                                                                minActiveForceLengthValue,...
                                                                plateauSlope, ...
                                                                curviness, ...
                                                                computeIntegral, ...
                                                                muscleName);  

    normMuscleCurves.activeForceLengthCurve = activeForceLengthCurve;

    %%
    %Fiber Passive Force
    %%
    eZero = 0; 
    eIso  = 0.7;
    kLow  = 0.2;
    kIso  = 2/(eIso-eZero);
    curviness = 0.75;
    computeIntegral = 1;

    fiberForceLengthCurve = createFiberForceLengthCurve(eZero,...
                                                        eIso,...
                                                        kLow,...
                                                        kIso,...
                                                        curviness,...
                                                        computeIntegral,...
                                                        muscleName);        
    normMuscleCurves.fiberForceLengthCurve = fiberForceLengthCurve;

    %%
    %Fiber Force Velocity Curve
    %%                                                
    fmaxE         = 1.4;
    dydxC         = 0;
    dydxNearC     = 0.15;
    dydxIso       = 5.0;
    dydxE         = 0.1;
    dydxNearE     = 0.15;
    flag_smoothenNonZeroDyDxC = 0;
    flag_usingOctave = 0;
    
    fvAtHalfVMax = 0.15; 
%   A fast/slow twitch force-velocity curve can be created by adjusting
%   fvAtHalfVMax:
%
%     Slow twitch: fvAtHalfVmax = 0.15
%     Fast twitch: fvAtHalfVmax = 0.22
%
%   These pameters come from Fig. 3a (fast) and Fig. 3b of Rantunga
%
%   Ranatunga, K. W. (1984). The force-velocity relation of rat
%   fast-and slow-twitch muscles examined at different temperatures.
%   The Journal of Physiology, 351, 517.

    concCurviness = 0.7;
    eccCurviness  = 0.9;
    computeIntegral = 0;

%     fiberForceVelocityCurve =  createFiberForceVelocityCurve2012(...
%                                             fmaxE, ...
%                                             dydxC,dydxNearC, ...
%                                             dydxIso, ...
%                                             dydxE, dydxNearE,...
%                                             concCurviness, eccCurviness,...
%                                             computeIntegral, muscleName);

     fiberForceVelocityCurve =  createFiberForceVelocityCurve2018(...
                                             fmaxE, ...
                                             dydxE, ...
                                             dydxC, ...
                                             flag_smoothenNonZeroDyDxC,...
                                             dydxNearE,...
                                             fvAtHalfVMax,...
                                             eccCurviness,...
                                             muscleName);


    normMuscleCurves.fiberForceVelocityCurve = fiberForceVelocityCurve;
    
    %Now make an invertable version of this curve - the slopes at the
    %end must be finite. 
    fiberForceVelocityCurveHACK =  createFiberForceVelocityCurve2018(...
                                             fmaxE, ...
                                             dydxNearE, ...
                                             dydxNearC, ...
                                             flag_smoothenNonZeroDyDxC,...
                                             dydxNearE,...
                                             fvAtHalfVMax,...
                                             eccCurviness,...
                                             muscleName);


    normMuscleCurves.fiberForceVelocityCurveHACK = fiberForceVelocityCurveHACK;    
    
    %Finally invert the curve. This is used to create a very accurate
    %initial guess for the Newton routine that is used to solve for the 
    %fiber velocity. This is quite useful during very slow eccentric
    %contractions where the force-velocity curve has a very sharp corner.
    fiberForceVelocityInverseCurveHACK = createInverseBezierCurve(fiberForceVelocityCurveHACK);
    
    normMuscleCurves.fiberForceVelocityInverseCurveHACK = fiberForceVelocityInverseCurveHACK;
    
    %%
    %Tendon Force Length Curve
    %%
    eIso            = 0.049;
    if(tendonStrainAtOneNormForceInput >= 0)
      eIso = tendonStrainAtOneNormForceInput;
    end
    
    kIso            = 1.375/eIso;
    fToe            = 2.0/3.0;
    curviness       = 0.5;
    computeIntegral = 1;

    tendonForceLengthCurve = ...
              createTendonForceLengthCurve( eIso, kIso, ...
                                            fToe, curviness, ...
                                            computeIntegral, ...
                                            muscleName);

    normMuscleCurves.tendonForceLengthCurve = tendonForceLengthCurve;

    %%
    % HACKED 'Classic' muscle modeling curves
    %
    % Almost all classic elastic tendon Hill-type muscle models take the
    % tendon-fiber equilibrium force equation:  
    %
    %    a(fl(lce)*fv(dlce/dt) + fpe(lce))*cos(alpha) - ft(lt) = 0;
    %
    % and massage it into an ode:
    %
    %    dlce/dt = fv^-1 ( [ (ft / a*cos(alpha)) - fpe] / fl)
    %
    % Which goes singular when ever:
    %
    %  a -> 0
    %  cos(alpha) -> 0
    %  fl -> 0
    %  fv^-1 -> inf as dlce/dt -> vmax
    %
    % To use this model without causing singularities means
    %
    %  a > 0      : the muscle cannot turn off          -> not physically true
    %  alpha < 90 : the pennation angle cannot go to 90 -> probably correct
    %  fl > 0     : the fiber can always generate force -> not physically true
    %  fv^-1      : the fiber can generate force at all
    %               velocities, and can even generate
    %               compressive forces at high shortening
    %               velocities.                          -> not true.
    %               As long as the fiber is connected to
    %               an elastic tendon with a force-length
    %               curve that does not go negative, then
    %               the fiber cannot apply compressive
    %               forces to the model. In the case of a
    %               rigid tendon however, there is a risk
    %               that compressive force are applied at
    %               very high shortening velocities (beyond
    %               vmax in shortening).
    %%

    %%
    %*HACKED Fiber Active Force Length Curve
    %%

    lce0 = 0.47-0.0259;
    lce1 = 0.73;
    lce2 = 1.0;
    lce3 = 1.8123;
    minActiveForceLengthValue = 0.1; %*Here's the hack: minimum value > 0
    curviness = 1.0;                 % in many papers this value really 
    computeIntegral = 0;             % is 0.1! This is huge!
    plateauSlope = 0.8616;

    activeForceLengthCurveHack = createFiberActiveForceLengthCurve( lce0,...
                                                                lce1, ...
                                                                lce2, ...
                                                                lce3, ...
                                                                minActiveForceLengthValue,...
                                                                plateauSlope, ...
                                                                curviness, ...
                                                                computeIntegral, ...
                                                                muscleName);

    normMuscleCurves.activeForceLengthCurveHACK = activeForceLengthCurveHack;
 
    %%
    %Generate plots of the curves using the function 
    % calcNormalizedMuscleCurveDerivative
    %%

    save('normMuscleCurves.mat','normMuscleCurves');
end

if(flag_plotCurves ==1)
    
    curveParamVector(1).struct = normMuscleCurves.activeForceLengthCurve;
    curvePlotInfo(1).xlabel = 'Norm. Fiber Length';
    curvePlotInfo(1).ylabel = 'Norm. Fiber Force';
    curvePlotInfo(1).xSymbol = 'lceN';
    curvePlotInfo(1).ySymbol = 'flN';
    curvePlotInfo(1).title = 'Active Force Length Curve';
            
    curveParamVector(2).struct = normMuscleCurves.fiberForceLengthCurve;
    curvePlotInfo(2).xlabel = 'Norm. Fiber Length';
    curvePlotInfo(2).ylabel = 'Norm. Fiber Force';
    curvePlotInfo(2).xSymbol = 'lceN';
    curvePlotInfo(2).ySymbol = 'fpeN';    
    curvePlotInfo(2).title = 'Force Length Curve';
        
    curveParamVector(3).struct = normMuscleCurves.fiberForceVelocityCurve;
    curvePlotInfo(3).xlabel = 'Norm. Fiber Velocity';
    curvePlotInfo(3).ylabel = 'Norm. Fiber Force';
    curvePlotInfo(3).xSymbol = 'dlceN/dt';
    curvePlotInfo(3).ySymbol = 'fvN';        
    curvePlotInfo(3).title = 'Force Velocity Curve';
    
    curveParamVector(4).struct = normMuscleCurves.tendonForceLengthCurve;
    curvePlotInfo(4).xlabel = 'Norm. Tendon Length';
    curvePlotInfo(4).ylabel = 'Norm. Tendon Force';
    curvePlotInfo(4).xSymbol = 'ltN';
    curvePlotInfo(4).ySymbol = 'ftN';            
    curvePlotInfo(4).title = 'Tendon Force Length Curve';

    curveParamVector(5).struct = normMuscleCurves.activeForceLengthCurveHACK;
    curvePlotInfo(5).xlabel = 'Norm. Fiber Length';
    curvePlotInfo(5).ylabel = 'Norm. Fiber Force';
    curvePlotInfo(5).xSymbol = 'lceN';
    curvePlotInfo(5).ySymbol = 'flN';                
    curvePlotInfo(5).title = 'Approx. Active Force Length Curve';
    
    curveParamVector(6).struct = normMuscleCurves.fiberForceVelocityCurveHACK;
    curvePlotInfo(6).xlabel = 'Norm. Fiber Velocity';
    curvePlotInfo(6).ylabel = 'Norm. Fiber Force';
    curvePlotInfo(6).xSymbol = 'dlceN/dt';
    curvePlotInfo(6).ySymbol = 'fvN';                    
    curvePlotInfo(6).title = 'Approx. Force Velocity Curve';
    
    curveParamVector(7).struct = normMuscleCurves.fiberForceVelocityInverseCurveHACK;
    curvePlotInfo(7).xlabel = 'Norm. Fiber Force';
    curvePlotInfo(7).ylabel = 'Norm. Fiber Velocity';
    curvePlotInfo(7).xSymbol = 'fvN';
    curvePlotInfo(7).ySymbol = 'dlceN/dt';                        
    curvePlotInfo(7).title = 'Approx. Force Velocity Inverse Curve';
    
    
    for i=1:1:length(curveParamVector)
        fig(i) = figure;

        curveSample = calcBezierYFcnXCurveSampleVector(...
                                        curveParamVector(i).struct, 100);


        xmin = min(curveSample.x);
        xmax = max(curveSample.x);

        ymin = min(curveSample.y);
        ymax = max(curveSample.y);
        yDelta = ymax-ymin;


        xV   = curveSample.x;
        yV   = curveSample.y;
        y1V  = curveSample.dydx;
        y2V  = curveSample.d2ydx2;

        %xStr = 'Norm. Fiber Length (lceN)';
        %yStr = 'Norm. Active Force (flN)';
        titleStr = curvePlotInfo(i).title;%curveParamVector(i).struct.name;


        subplot(2,2,1);
            plot(xV,yV,'b');
            xlabel(curvePlotInfo(i).xSymbol);
            ylabel(curvePlotInfo(i).ySymbol);
            title(titleStr);
            grid on;
            hold on;
                xlim([xmin,xmax]);
                ylim([ymin-0.1*yDelta,ymax+0.1*yDelta]);

        subplot(2,2,2);        
            plot(xV,y1V,'r');
            xlabel(curvePlotInfo(i).xSymbol);
            ylabel(['d/',curvePlotInfo(i).xSymbol,' ',curvePlotInfo(i).ySymbol]);
            title(['d/',curvePlotInfo(i).xSymbol,' ',titleStr]);

            grid on;
            hold on;

                xlim([xmin,xmax]);                                

        subplot(2,2,3);        
            plot(xV,y2V,'m');
            xlabel(curvePlotInfo(i).xSymbol);
            ylabel(['d^2/',curvePlotInfo(i).xSymbol,'^2 ',curvePlotInfo(i).ySymbol]);
            title(['d^2/',curvePlotInfo(i).xSymbol,'^2 ',titleStr]);
            grid on;
            hold on;
            xlim([xmin,xmax]);       

        if(isempty(curveParamVector(i).struct.integral)==0)  
            intYdx = curveSample.intYdx;

            subplot(2,2,4)
                plot(xV, intYdx,'g');
                xlabel(curvePlotInfo(i).xSymbol);
                ylabel(['intYdx ',curvePlotInfo(i).ySymbol]);
                title(['intYdx ',titleStr]);
                grid on;
                hold on;
                xlim([xmin,xmax]);  

        end

    end
end
    





            