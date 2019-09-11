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
function fiberForceVelocityInverseCurve = ...    
    createFiberForceVelocityInverseCurve2012(fmaxE, ...
                                         dydxC, dydxNearC, ...
                                         dydxIso,...
                                         dydxE, dydxNearE,...
                                         concCurviness, eccCurviness,...
                                         computeIntegral, muscleName)
%%

%This function will generate a C2 continuous (continuous to the 2nd
%derivative) inverse curve that the function 
%createFiberForceVelocityCurve generates. The inverse force velocity 
%curve is required by every equilibrium muscle model in order to compute
%the derivative of fiber velocity. To generate the inverse force velocity
%curve simply call this function with EXACTLY the same parameter values
%that you used to generate the force velocity curve. See the parameter
%descriptions for createFiberForceVelocityCurve, as the parameters for
%the inverse function are identical. The curve name should be different,
%however, because this is an inverse curve 
%(e.g. 'bicep_fiberForceVelocityInverseCurve')
%
%
%\image html fig_SmoothSegmentedFunctionFactory_fvInvCurve.png
%%
fiberForceVelocityInverseCurve = [];
fiberForceVelocityInverseCurve.name = ...
    sprintf('%s.%s',muscleName,'fiberForceVelocityInverse');


%Ensure that the inputs are within a valid range
assert( fmaxE > 1.0, ...
    sprintf('%s: fmaxE must be greater than 1',...
    fiberForceVelocityInverseCurve.name)); 

assert( (dydxC >= eps^0.5 && dydxC < 1), ...
    sprintf('%s: dydxC must be greater than eps^0.5 and less than 1',...
    fiberForceVelocityInverseCurve.name));

assert( (dydxNearC >= dydxC && dydxNearC < 1), ...
    sprintf('%s: dydxNearC must be greater than 0 and less than 1',...
        fiberForceVelocityInverseCurve.name));
    
assert( dydxIso > 1, ...
    sprintf('%s: dydxIso must be greater than or equal to 1',...
    fiberForceVelocityInverseCurve.name));

assert( (dydxE >= eps^0.5 && dydxE < (fmaxE-1)),  ...   
    sprintf('%s: dydxE must be greater than or equal to eps^0.5 and less than fmaxE-1 (%f)',...
    fiberForceVelocityInverseCurve.name,(fmaxE-1)));

assert( (dydxNearE >= dydxE && dydxNearE < (fmaxE-1)), ...   
    sprintf('%s: dydxNearE must be greater than or equal to dydxE and less than fmaxE-1 (%f)',...
    fiberForceVelocityInverseCurve.name,(fmaxE-1)));

assert( (concCurviness <= 1.0 && concCurviness >= 0),  ...  
    sprintf('%s: concCurviness must be between 0 and 1',...
    fiberForceVelocityInverseCurve.name));

assert( (eccCurviness <= 1.0 && eccCurviness >= 0),    ...
    sprintf('%s: eccCurviness must be between 0 and 1',...
    fiberForceVelocityInverseCurve.name));


%Translate the users parameters into Bezier point locations
cC = scaleCurviness(concCurviness);
cE = scaleCurviness(eccCurviness);

%Compute the concentric control point locations
xC   = -1;
yC   = 0;

xNearC = -0.9;
yNearC = yC + 0.5*dydxNearC*(xNearC-xC) + 0.5*dydxC*(xNearC-xC);

xIso = 0;
yIso = 1;

xE   = 1;
yE   = fmaxE;

xNearE = 0.9;
yNearE = yE + 0.5*dydxNearE*(xNearE-xE) + 0.5*dydxE*(xNearE-xE);


concPts1 = ...
    calcQuinticBezierCornerControlPoints(xC,     yC,    dydxC, 0,...
                                     xNearC, yNearC,dydxNearC, 0, cC);
                                        
concPts2 = ...
    calcQuinticBezierCornerControlPoints(xNearC,yNearC,dydxNearC, 0,...
                                           xIso,  yIso,  dydxIso, 0, cC);
                                       
eccPts1 = ...
    calcQuinticBezierCornerControlPoints(xIso,   yIso,   dydxIso, 0,... 
                                       xNearE, yNearE, dydxNearE, 0, cE);

eccPts2 = ...
    calcQuinticBezierCornerControlPoints(xNearE, yNearE, dydxNearE, 0,... 
                                             xE,     yE,     dydxE, 0, cE);

                                        
xpts = [concPts1(:,1) concPts2(:,1) eccPts1(:,1) eccPts2(:,1)];
ypts = [concPts1(:,2) concPts2(:,2) eccPts1(:,2) eccPts2(:,2)];

%Note the xy inversion! 
fiberForceVelocityInverseCurve.xpts = ypts;
fiberForceVelocityInverseCurve.ypts = xpts;

fiberForceVelocityInverseCurve.xEnd = [yC yE];
fiberForceVelocityInverseCurve.yEnd = [xC xE];
fiberForceVelocityInverseCurve.dydxEnd= [1/dydxC 1/dydxE];
fiberForceVelocityInverseCurve.d2ydx2End=[0,0];

fiberForceVelocityInverseCurve.integral = [];

if(computeIntegral == 1)
    fiberForceVelocityInverseCurve.integral = ...
        createCurveIntegralStructure(fiberForceVelocityInverseCurve, 1000, 1e-12,NaN);    
end


        
