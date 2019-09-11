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
function mtInfo = calcMillard2012DampedEquilibriumMuscleInfo(...
                                            activationState,...                                            
                                            pathState, ... 
                                            muscleState,...                                                                                       
                                            muscleArchitecture,...
                                            normMuscleCurves,...
                                            modelConfig)
%%
% This function calculates the kinematics, forces, powers, stored elastic
% energy of the elements in a lumped parameter muscle model. This muscle
% model can be configured such that the tendon is rigid or elastic, and
% also so that the classic state equations are used or the updated damped
% equilibrium equations are used to model musculotendon dynamics. This is a
% Matlab port of the model described in 
%
%  Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). 
%    Flexing computational muscle: modeling and simulation of 
%    musculotendon dynamics. Journal of biomechanical engineering, 
%    135(2), 021005.
% 
% @param activationState: scalar of the muscle activation [0,1]
%
% @param pathState: 2x1 vector 
%                   pathState(1) : lengthening velocity of the path in m/s
%                   pathState(2) : length of the path in m
%
% @param muscleState: this is dependent on the settings in modelConfig:
%        If modelConfig.useElasticTendon is 0: this field is ignored
%
%        If modelConfig.useElasticTendon is 0: then
%           muscleState is a scalar of fiberLengthAlongTendon
%
%        If muscleState is a 2x1 vector of 
%          [fiberVelocityAlongTendon; fiberLengthAlongTendon]
%           then the model is put into initialization mode and the value of
%           the equilibrium equation is evaluated, rather than using it to
%           solve for the fiber velocity. This is ugly, I know.
%
% @param muscleArchtecture: a structure containing the architecture
%        information about the muscle. This structure has fields of:
%
%
%     .name    : string name of the muscle. e.g. 'soleus leftN'                      
%
%     .abbr    : abbreviated name e.g. 'solL'
%
%     .fiso    : maximum isometric force (n) of the muscle when the fiber 
%                is at its optimal length, and it is static. e.g. 3800               
%
%     .optimalFiberLength : length (m) of the muscle fiber when it 
%                           generates its maximum isometric force 
%                           e.g. 0.044m
%
%     .maximumNormalizedFiberVelocity : maximum velocity the fiber can
%                              contract or lengthen at in units 
%                              of fiber lengths per second. This 
%                              is normally set to 10 fiber lengths/second.
%
%     .pennationAngle : the angle (radians) between the fiber and the 
%                       tendon when the fiber is at its optimal length.
%
%     .tendonSlackLength : the unloaded length (m) of tendon.
% 
%     .minimumFiberLength : the minimum physical length (m) the fiber is 
%                           allowed to achieve. 
%
%     .minimumFiberLengthAlongTendon : the minimum length (m) of the fiber 		 
%                                      projected along the tendon
%
%     .pennationAngleAtMinumumFiberLength : the pennation angle (radians) 
%                                  of the fiber when it is at its minimum
%                                  length.
%
% @param normMuscleCurves: a structure containing the parameters needed to
%                          evaluate the normalized muscle curves (below) 
%                          using the function calcBezierYFcnXDerivative.
%                          This structure must have the fields
%
%  .activeForceLengthCurve      : gen. by createActiveForceLengthCurve.m
%
%  .activeForceLengthCurveHACK  : gen. by createActiveForceLengthCurve.m 
%                                 s.t. the minimum value is > 0
%
%  .fiberForceVelocityCurve: gen. by createFiberForceVelocityCurve.m                                 
%
%  .fiberForceVelocityCurveHACK: gen. by createFiberForceVelocityCurve.m 
%                                s.t. d(fvN)/d(dlceN) > 0
%
%  .fiberForceVelocityInverseCurveHACK gen. 
%                             by createFiberForceVelocityInverseCurve.m 
%                             s.t. d(fvN)/d(dlceN) > 0
%
%  .fiberForceLengthCurve: gen. by createFiberForceLengthCurve.m
%
%  .tendonForceLengthCurve: gen by createTendonForceLengthCurve.m
%
% @param modelConfig: a structure with the fields
%
%     .useFiberDamping  : 0: no damping element in the fiber
%                         1: Adds a damping element into the fiber. If the
%                            tendon is elastic then using fiber damping
%                            also changes the way the fiber velocity is
%                            calculated to a faster method
%
%     .useElasticTendon = 0: rigid tendon is used, and the model is
%                            stateless. If tendonSlackLength is less
%                            than the optimalFiberLength, this is a good
%                            option; otherwise a compromise is made on
%                            the accuracy of the predicted fiber forces
%                         1: the tendon element is treated as being elastic
%
%     .damping          = 0.1: the amount of damping to apply to the fiber
%                         if useFiberDamping is used. Note that the 
%                         damping force is calculated as damping*dlceN, 
%                         and so this coefficent represents the maximum
%                         damping force ever applied to the fiber.
%
%     .iterMax          = maximum number of Newton iterations permitted
%                         when solving the damped equilibrium equation
%
%     .tol              = normalized force tolerance used when solving the
%                         damped equilibrium equation.
%
%     .minActivation    = minimum activation level permitted. When the
%                         classic elastic model is used (useFiberDamping 
%                         = 0 and useElasticTendon = 1) this value must be 
%                        greater than 0.
%
% @returns mtInfo: A giant structure containing just about everything you
% might want to know about the muscle:
%
%
% mtInfo.muscleLengthInfo.fiberLength                 	  %length            m  
% mtInfo.muscleLengthInfo.fiberLengthAlongTendon      	  %length            m
% mtInfo.muscleLengthInfo.normFiberLength             	  %length/length     m/m        
% mtInfo.muscleLengthInfo.tendonLength                	  %length            m
% mtInfo.muscleLengthInfo.normTendonLength            	  %length/length     m/m        
% mtInfo.muscleLengthInfo.tendonStrain                	  %length/length     m/m        
% mtInfo.muscleLengthInfo.pennationAngle              	  %angle             1/s    
% mtInfo.muscleLengthInfo.cosPennationAngle           	  %NA                NA         
% mtInfo.muscleLengthInfo.sinPennationAngle           	  %NA                NA         	
% mtInfo.muscleLengthInfo.fiberPassiveForceLengthMultiplier %NA                NA
% mtInfo.muscleLengthInfo.fiberActiveForceLengthMultiplier  %NA                NA
%           
% mtInfo.fiberVelocityInfo.fiberVelocity                  %length/time           m/s
% mtInfo.fiberVelocityInfo.fiberVelocityAlongTendon       %length/time           m/s
% mtInfo.fiberVelocityInfo.normFiberVelocity              %(length/time)/length (m/s)/m
% mtInfo.fiberVelocityInfo.pennationAngularVelocity       %angle/time            rad/s
% mtInfo.fiberVelocityInfo.tendonVelocity                 %length/time           m/s
% mtInfo.fiberVelocityInfo.normTendonVelocity             %(length/time)/length  (m/s)/m
% mtInfo.fiberVelocityInfo.fiberForceVelocityMultiplier   %force/force           NA
% 
% mtInfo.muscleDynamicsInfo.activation                	% of muscle active   NA						
% mtInfo.muscleDynamicsInfo.fiberForce                	% force                N
% mtInfo.muscleDynamicsInfo.fiberForceAlongTendon     	% force                N
% mtInfo.muscleDynamicsInfo.normFiberForce            	% force/force          N/N
% mtInfo.muscleDynamicsInfo.activeFiberForce          	% force                N
% mtInfo.muscleDynamicsInfo.passiveFiberForce         	% force                N
% mtInfo.muscleDynamicsInfo.tendonForce               	% force                N
% mtInfo.muscleDynamicsInfo.normTendonForce           	% force/force          N/N
% mtInfo.muscleDynamicsInfo.fiberStiffness            	% force/length         N/m
% mtInfo.muscleDynamicsInfo.fiberStiffnessAlongTendon 	% force/length         N/m
% mtInfo.muscleDynamicsInfo.tendonStiffness           	% force/length         N/m
% mtInfo.muscleDynamicsInfo.muscleStiffness           	% force/length         N/m
% mtInfo.muscleDynamicsInfo.fiberActivePower          	% power   			   W
% mtInfo.muscleDynamicsInfo.fiberPassivePower         	% power   			   W
% mtInfo.muscleDynamicsInfo.tendonPower               	% power   			   W
% mtInfo.muscleDynamicsInfo.musclePower               	% power   			   W
% mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy   %energy   J   
% mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy  %energy    J   
% mtInfo.musclePotentialEnergyInfo.musclePotentialEnergy  %energy    J (Nm)
% 
% %This field does not appear in the OpenSim model
% mtInfo.muscleDynamicsInfo.fiberParallelElementPower 	% power   			   W
% mtInfo.muscleDynamicsInfo.dampingForces             	% force   			   N
% mtInfo.muscleDynamicsInfo.dampingPower              	% power   			   W
% mtInfo.muscleDynamicsInfo.boundaryPower             	% power   			   W
% 
% mtInfo.state.value      
% mtInfo.state.derivative 
% 
% mtInfo.initialization.err                               % force/force           N/N
% mtInfo.initialization.Derr_DlceAT                      % (force/force) / (m/s) (N/N)/(m/s)
%
%%
mtInfo = [];
epsRoot = eps^0.5;

%%
%Get the model configuration
%%

useFiberDamping  = modelConfig.useFiberDamping;
useElasticTendon = modelConfig.useElasticTendon;

if(useFiberDamping == 0 && useElasticTendon == 1)
   assert(modelConfig.minActivation > 0,...
          ['modelConfig.minActivation must be greater than zero',...
           ' when the model is configured as a classic elastic',...
           ' muscle model, otherwise a singularity is possible']);
end

useInitializationCode = 0;
if(isempty(muscleState) == 0)
    if(length(muscleState) == 2)
        useInitializationCode = 1;
    end
end

beta = 0;
if(useFiberDamping == 1)
   beta = modelConfig.damping; 
end

%%
%Break out the input state vectors and model
%parameters into local variables
%%



a     = clampActivation(activationState(1), ...
                        modelConfig.minActivation,...
                        1);

dlp = pathState(1);
lp  = pathState(2);

fiso      =  muscleArchitecture.fiso;
lceOpt    =  muscleArchitecture.optimalFiberLength;
alphaOpt  =  muscleArchitecture.pennationAngle;
ltSlk     =  muscleArchitecture.tendonSlackLength;
dlceMaxN  =  muscleArchitecture.maximumNormalizedFiberVelocity;

lceMin       = muscleArchitecture.minimumFiberLength;
lceATMin     = muscleArchitecture.minimumFiberLengthAlongTendon;
alphaLceMin  = muscleArchitecture.pennationAngleAtMinumumFiberLength;

%Create function handles for the curves to make the code easier to read

calcFalDer = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                       normMuscleCurves.activeForceLengthCurve, ...
                       arg2);
                           
calcFalHACKDer = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                       normMuscleCurves.activeForceLengthCurveHACK, ...
                       arg2);
                           
calcFvDer = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                           normMuscleCurves.fiberForceVelocityCurve, ...
                           arg2);

calcFvHACKDer = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                       normMuscleCurves.fiberForceVelocityCurveHACK, ...
                       arg2);
                           
calcFvInvHACKDer = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
               normMuscleCurves.fiberForceVelocityInverseCurveHACK, ...
               arg2);                           
                           
calcFpeDer = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                           normMuscleCurves.fiberForceLengthCurve, ...
                           arg2);
                           
calcFtDer = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                           normMuscleCurves.tendonForceLengthCurve, ...
                           arg2);  
                           
                           
%%
%0. Check the model inputs
%%
assert( a >= 0       && a <= 1     ,   'Check activation model should be [0, 1]');
assert( lp < 1       && lp >= 0    ,   'Check Units: path length in m!');
assert( alphaOpt < pi/2            ,   'Check Units: Pennation angle must be in radians!');
assert( lceOpt < 0.5 && lceOpt > 0 ,   'Check Units: fiber length in m!');
assert( ltSlk < 0.5  && ltSlk >= 0 ,   'Check Units: tendon slack length in m!');



%%
%1. Calculate the fiber velocity by satisfying
%   the equilibrium equation
%%

%%
%Short hand variable name conventions
%
% Quantity
%   l : length
%   e : strain
%   d : d/dt
%   D : a partial derivative
%   f : force
%   k : stiffness
%   
%   alpha: pennation angle of the fiber w.r.t. tendon
%
% Parts
%   p  : path, as in the path the musuclotendon follows.
%   ce : contractile element
%   pe : parallel element
%   t  : tendon
%
% Modifier
%   N  : normalized. what this means is dependent on the quantity
%        for a precise definition refer to the Millard 2010 paper
%   AT : (a)long (t)endon. The projection of the quantity along the
%        direction of the tendon.
%
%   Slk: slack length 
%%

%Fiber and tendon lengths
lce    = NaN; 
lceAT  = NaN;
alpha  = NaN;
lt     = NaN;
lceN   = NaN;
ltN    = NaN;
etN    = NaN;

sinAlpha = NaN;
cosAlpha = NaN;

%Fiber and tendon velocities
dlce   = NaN;
dlceAT = NaN;
dalpha = NaN;
dlt    = NaN;
dlceN  = NaN;
dltN   = NaN;

%Normalized musculotendon curve values
falN        = NaN;
DfalN_DlceN = NaN;
fvN         = NaN;
fpeN        = NaN;
DfpeN_DlceN = NaN;
ftN         = NaN;
DftN_DltN   = NaN;
  
%Forces, Stiffnesses
ffaN       = NaN; %(f)orce (f)iber (a)ctive normalized 
ffpN       = NaN; %(f)orce (f)iber (p)assive normalized
ffN        = NaN; %fiber force

kt         = NaN; %tendon stiffness in N/m
kf         = NaN; %fiber stiffness in N/m
kfAT       = NaN; %fber stiffness along the tendon
km         = NaN; %stiffness of the entire muscle

%Equations for the initialization error
init              = [];
init.err          = NaN;
init.Derr_DlceAT  = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                       Rigid tendon muscle model
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(useElasticTendon == 0)

   
   lceAT  = lp  - ltSlk;
   lt     = ltSlk;
   
   dlt    = 0;
   dlceAT = dlp - dlt; 

   %Clamp the fiber if its at its lower bound
   fiberState = clampFiberStateAlongTendon(lceAT,dlceAT,lceATMin);
   lceAT      = fiberState.lceAT;
   dlceAT     = fiberState.dlceAT;
   isClamped  = fiberState.isClamped;   
   
   fiberKinematics = calcFixedWidthPennatedFiberKinematics(lceAT,...
                                                           dlceAT,...
                                                           lceOpt,...
                                                           alphaOpt);
   %Fiber & tendon lengths
   lce    = fiberKinematics.fiberLength;
   alpha  = fiberKinematics.pennationAngle;
   sinAlpha= sin(alpha);
   cosAlpha= cos(alpha);
   lt     = ltSlk;       
   lceN   = lce/lceOpt;
   ltN    = 1;
   etN    = 0;

   %Fiber and tendon velocities
     
   dlce   = fiberKinematics.fiberVelocity;
   dalpha = fiberKinematics.pennationAngularVelocity;       
   dlt    = dlp - dlceAT;       
   dlceN  = dlce/(lceOpt*dlceMaxN);
   dltN   = dlt/ltSlk;    
      
   %%
   %3. Evaluate the normalized muscle curves 
   %%
   
   falN        = calcFalDer(lceN,0);
   DfalN_DlceN = calcFalDer(lceN,1);

   fpeN        =  calcFpeDer(lceN,0);
   DfpeN_DlceN =  calcFpeDer(lceN,1);

   fvN         =  calcFvDer(dlceN, 0);


   ffaN       = a*falN*fvN;        %(f)orce (f)iber (a)ctive normalized 
   ffpN       = fpeN + beta*dlceN; %(f)orce (f)iber (p)assive normalized
   ffN        = ffaN + ffpN;       %fiber force
   
   ftN        = ffN*cosAlpha;    
   DftN_DltN  = inf;

   %d/dlce fiso*(a*falN*fvN + fpeN) 
   %           = fiso*(a*dFalNdlceN*fvN + dFpeNdlceN)*(dlceN/dlce)
   %
   kf = fiso*(a*DfalN_DlceN*fvN + DfpeN_DlceN)*(1/lceOpt);
   
   %d/dlce fiso*(a*falN*fvN + fpeN)*cos(alpha)
   %     = kf*cos(alpha) + fiso*ffN*dcosAlphaDlce
   %
   Dalpha_Dlce = calcFixedWidthPennationDalphaDlce(alpha,...
                                                  lce,...
                                                  lceOpt,...
                                                  alphaOpt);
   kfAT      = kf*cos(alpha) - fiso*ffN*sinAlpha*Dalpha_Dlce;                                              
   kt        = inf;
   km        = kfAT; %Since the tendon is infinitely stiff;
   
   if(isClamped == 1)
      ftN = 0; 
      kf = 0;
      kfAT = 0;
      kt   = 0;
      km   = 0;
   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                       Elastic tendon muscle model(s)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(useElasticTendon == 1 )

        
    if(useInitializationCode == 1)
        dlceAT = muscleState(1);
        lceAT  = muscleState(2);
    else
        dlceAT= NaN;
        lceAT = muscleState(1);
    end
    
    
    %Clamp the fiber 
    fiberState = clampFiberStateAlongTendon(lceAT,dlceAT,lceATMin);
    lceAT      = fiberState.lceAT;
    dlceAT     = fiberState.dlceAT;
    isClamped  = fiberState.isClamped;
    
    %Calculate the pennation angle and fiber
    fiberKinematics = calcFixedWidthPennatedFiberKinematics(lceAT,...
                                                            dlceAT,...
                                                            lceOpt,...
                                                            alphaOpt);
    lce     = fiberKinematics.fiberLength;
    dlce    = fiberKinematics.fiberVelocity;
    alpha   = fiberKinematics.pennationAngle;
    dalpha  = fiberKinematics.pennationAngularVelocity;
    
    lceN    = lce/lceOpt;    
    cosAlpha= cos(alpha);
    sinAlpha= sin(alpha);
    
    lt = lp-lceAT;
    ltN= lt/ltSlk;
    etN= ltN-1;

    %Normalized musculotendon curve values that are common to
    %both the damped and inverted musculotendon models.

    fpeN       =  calcFpeDer(lceN,0);
    DfpeN_DlceN =  calcFpeDer(lceN,1);

    ftN        =  calcFtDer(ltN, 0);
    DftN_DltN   =  calcFtDer(ltN, 1);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %
    %                    Classic elastic tendon muscle model
    %    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(useFiberDamping == 0)
        %Compute fiber velocity by inverting the force-velocity curve 
        %How? Here's the force equilibrium equation
        % 
        % (a*falN*fvN + fpeN)*cos(alpha) - ftN = 0;
        %
        %  which we now solve for dlceN
        %
        %  dlceN = fvNInv( [ftN/cos(alpha) - fpeN]/(a*falN) )
        %
        % Which can go singular in 4 different ways
        %  cos(alpha) -> 0
        %  a          -> 0
        %  falN       -> 0
        %  fvNInv     -> +/- inf as dlceN -> +/- 1 
        %
        % First lets evaluate the muscle curves using the HACK versions
        % for the active force length and the force velocity curve. Note
        % that the HACK versions are different than normal:
        %  
        % Hacked falN > 0
        % Hacked fvN is invertible: slope does not go to zero at the ends
        %
        falN       =  calcFalHACKDer(lceN, 0);
        DfalN_DlceN = calcFalHACKDer(lceN, 1);
        
        % Now check for singularities        
        assert(cos(alpha) > eps^0.5,...
          'Cannot invert fvN curve: cos(pennationAngle) <= eps^0.5');
        assert(a     > eps^0.5,...
          'Cannot invert fvN curve: activation <= eps^0.5');
        assert(falN > eps^0.5,...
          ['Cannot invert fvN curve:',...
           ' forceActiveLengthMultiplier <= eps^0.5']);

        %%
        %Evaluate the fiber velocity by solving for fvN and inverting
        %the force-velocity curve
        %%
        if(useInitializationCode == 0)
            fvN    = ((ftN/cos(alpha)) - fpeN) / (a*falN);
            dlceN = calcFvInvHACKDer(fvN,0);       
            dlce  = dlceN*(lceOpt*dlceMaxN);
        end
        
        %Clamp the fiber velocity if necessary
        fiberState = clampFiberState(lce,dlce,lceMin);        
        dlce  = fiberState.dlce;
        dlceN = dlce/(lceOpt*dlceMaxN);
        fvN    = calcFvHACKDer(dlceN,0);
        
        ffaN  = a*falN*fvN;
        ffpN  = fpeN;
        ffN   = ffaN + ffpN;
        
        %d/dlce (a*falN*fvN + fpeN + beta*dlceN)
        %        a*DfalN_DlceN*fvN + DfpeN_DlceN)*dlceNdlce
        kf         = (a*DfalN_DlceN*fvN + DfpeN_DlceN)*(fiso/lceOpt);
        
        Dalpha_Dlce = calcFixedWidthPennationDalphaDlce( alpha,...
                                                        lce,...
                                                        lceOpt,...
                                                        alphaOpt);
        % d/dlce    fiso*ffN*cos(alpha)
        %     =     kf*cos(alpha) - fiso*ffN*sin(alpha)*Dalpha_Dlce 
        %        
        kfAT       = kf*cosAlpha - fiso*ffN*sinAlpha*Dalpha_Dlce;         
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %
    %                    Damped fiber elastic tendon muscle model
    %    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    elseif(useFiberDamping == 1)        
        
        falN       = calcFalDer(lceN,0);
        DfalN_DlceN = calcFalDer(lceN,1);
        
        iterMax = modelConfig.iterMax;
        tol     = modelConfig.tol;

        err     = 10*tol;
        iter    = 1;

        Ddlce_DdlceN = (lceOpt*dlceMaxN);
        
        Dalpha_Dlce = calcFixedWidthPennationDalphaDlce(alpha,...
                                                       lce,...
                                                       lceOpt,...
                                                       alphaOpt);
                                                   
        if(useInitializationCode == 0)
                                                
            fvN    = ((ftN/cos(alpha)) - fpeN) / (a*falN + 1e-8);
            dlceN = calcFvInvHACKDer(fvN,0);       
            if(dlceN > 1)
              dlceN = a*0.99;
            end
            if(dlceN < -1)
              dlceN = -a*0.99;
            end
                     

            %%
            %use Newton's method to polish the root to high precision
            %%
            while abs(err) > tol && iter < iterMax 
          
                fvN         = calcFvDer(dlceN,0);
                DfvN_DdlceN = calcFvDer(dlceN,1);
                
                %Tension is +
                %shortening is -
                % ... thus the sign of beta is +
                ffN         = a*(falN*fvN) + fpeN + beta*dlceN;
                DffN_DdlceN = (a*(falN*DfvN_DdlceN) + beta);

                err = ffN*cosAlpha - ftN;
                Derr_DdlceN  = DffN_DdlceN*cosAlpha;

                if(abs(err) > tol && abs(Derr_DdlceN) > eps*10)               
                     delta = -err/Derr_DdlceN;
                     dlceN = dlceN + delta;
                                          
                elseif( abs(err) > tol && abs(Derr_DdlceN) <= eps*10)
                    %This code should never be entered given that
                    %the smallest gradient is beta*cos(maxPennAngle) 
                    %which is approx 0.0017 when beta is 0.1 and
                    %the max pennation angle is 89 degrees. If the
                    %parameters are chosen to be particularly aggressive
                    %the gradient may hit zero as dlceN -> dlceMaxN. In
                    %this case we set dlceN to the ends of the fvN curve
                    %that have a nonzero slope and hope for the best.                    
                    if(dlceN < dlceNNearMaxShortening)
                        dlceN = dlceNNearMaxShortening;
                    elseif(dlceN > dlceNNearMaxLengthening)
                        dlceN = dlceNNearMaxLengthening; 
                    else
                       msg=['Damped equilibrium Newton method has a\n',...
                           ' zero gradient for an unexpected reason\n',...
                           ' lceNMin %e \n lceN %e ',...
                           ' \n alpha %e  \n dlceN %e\n'];
                       assert(0,sprintf(msg,lceMin,lceN,alpha,dlceN)); 
                    end

                end
                iter = iter+1;
            end

            %Clamp the fiber velocity if necessary
            dlce       = dlceN * lceOpt * dlceMaxN;
            fiberState = clampFiberState(lce,dlce,lceMin); 
            isClamped  = fiberState.isClamped;
            dlce       = fiberState.dlce;
            
            if(isClamped == 0 && abs(err) > tol)
               msg=['Damped equilibrium eqn not satisified\n',...
                       ' lceNMin %e \n lceN %e ',...
                       ' \n alpha %e  \n dlceN %e\n'];
                   
               smsg  = sprintf(msg,lceMin,lceN,alpha,dlceN); 
               
               assert(0, smsg);             
            end

        end
        
        dlceN  = dlce/(lceOpt*dlceMaxN);
        fvN    = calcFvDer(dlceN, 0);
        
        ffaN       = a*(falN*fvN);
        ffpN       = fpeN + beta*dlceN;
        ffN        = ffaN + ffpN;
        
        %d/dlce (a*falN*fvN + fpeN + beta*dlceN)
        %        a*DfalN_DlceN*fvN + DfpeN_DlceN)*dlceNdlce
        kf         = fiso*(a*DfalN_DlceN*fvN + DfpeN_DlceN)*(1/lceOpt);
        
        Dalpha_Dlce = calcFixedWidthPennationDalphaDlce( alpha,...
                                                        lce,...
                                                        lceOpt,...
                                                        alphaOpt);
        % d/dlce    fiso*ffN*cos(alpha)
        %     =     kf*cos(alpha) - fiso*ffN*sin(alpha)*Dalpha_Dlce 
        %        
        kfAT       = kf*cosAlpha - fiso*ffN*sinAlpha*Dalpha_Dlce;         
    end
       
    %Fiber and tendon velocities       
    fiberKinematics = ...
        calcFixedWidthPennatedFiberKinematicsAlongTendon(lce,...
                                                         dlce,...
                                                         lceOpt,...
                                                         alphaOpt);
    dlceAT = fiberKinematics.fiberVelocityAlongTendon;        
    dalpha = fiberKinematics.pennationAngularVelocity;
        
    dlt    = dlp - dlceAT;                              
    dltN   = dlt/ltSlk;
        
    ktN = DftN_DltN;
    kt  = ktN * (fiso / ltSlk);

    %This can blow up: look at the denominator
    %                  then go look at the slope of 
    %                  the active force length curve and the fiber force
    %                  length curve.
    km  = kfAT*kt / (kt + kfAT);    

    %%
    %Initialization Equations
    %  We need:
    %      ffAT - ftN = 0
    %  And the partial derivative        
    %      DffAT_DlceAT - DftN_DlceAT
    %
    %%        
    init.err        = fiso*( (ffaN + ffpN)*cosAlpha - ftN );
    
    Dff_DlceAT       = kfAT;        
    DftN_Dlt          = kt;        
    % lt         = lp - lceAT
    % dlt_dlceAT = -1
    Dlt_DlceAT   = -1;    
    DftN_DlceAT   = DftN_Dlt * Dlt_DlceAT;
    
    init.Derr_DlceAT = Dff_DlceAT - DftN_DlceAT;
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4. Populate the muscle length and velocity information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtInfo.muscleLengthInfo.fiberLength                  = lce;        %length            m  
mtInfo.muscleLengthInfo.fiberLengthAlongTendon       = lceAT;      %length            m
mtInfo.muscleLengthInfo.normFiberLength              = lceN;       %length/length     m/m        
mtInfo.muscleLengthInfo.tendonLength                 = lt;         %length            m
mtInfo.muscleLengthInfo.normTendonLength             = ltN;        %length/length     m/m        
mtInfo.muscleLengthInfo.tendonStrain                 = etN;        %length/length     m/m        
mtInfo.muscleLengthInfo.pennationAngle               = alpha;      %angle             1/s    
mtInfo.muscleLengthInfo.cosPennationAngle            = cosAlpha;   %NA                NA         
mtInfo.muscleLengthInfo.sinPennationAngle            = sinAlpha;   %NA                NA         	
mtInfo.muscleLengthInfo.fiberPassiveForceLengthMultiplier  = fpeN;  %NA             NA
mtInfo.muscleLengthInfo.fiberActiveForceLengthMultiplier   = falN;  %NA             NA
          
mtInfo.fiberVelocityInfo.fiberVelocity                = dlce;    %length/time           m/s
mtInfo.fiberVelocityInfo.fiberVelocityAlongTendon     = dlceAT;  %length/time           m/s
mtInfo.fiberVelocityInfo.normFiberVelocity            = dlceN;   %(length/time)/length  (m/s)/m
mtInfo.fiberVelocityInfo.pennationAngularVelocity     = dalpha;  %angle/time            rad/s
mtInfo.fiberVelocityInfo.tendonVelocity               = dlt;     %length/time           m/s
mtInfo.fiberVelocityInfo.normTendonVelocity           = dltN;    %(length/time)/length  (m/s)/m
mtInfo.fiberVelocityInfo.fiberForceVelocityMultiplier = fvN;     %force/force           NA


%%
%7. Evaluate the potential energy stored in the musculotendon
%%
pt  = 0;
ppe = 0;

if(isempty(normMuscleCurves.tendonForceLengthCurve.integral) == 0)
    if(useElasticTendon == 1)
       ptN = calcFtDer(ltN, -1); 
       
       tendonStrainAtOneNormForce = ...
           normMuscleCurves.tendonForceLengthCurve.integral.xScaling;
       tendonStretchAtOneNormForce = tendonStrainAtOneNormForce*ltSlk;
       
       %To get this scaling we're getting calculating the energy in
       %a square that is 100% strain by 1 maximum isometric force
       %      first in units of Nm 
       %      second in normalized units
       %And then taking the ratio to obtain the scaling
       tendonEnergyScaling = tendonStretchAtOneNormForce*fiso ...
                           / tendonStrainAtOneNormForce*1;

       pt = ptN*tendonEnergyScaling;

    end
end

if(isempty(normMuscleCurves.fiberForceLengthCurve.integral) == 0)
    
    ppeN = calcFpeDer(lceN,-1);
    
    fiberStrainAtOneNormForce = ...
           normMuscleCurves.fiberForceLengthCurve.integral.xScaling;

    fiberStretchAtOneNormForce = fiberStrainAtOneNormForce*lceOpt;

   %To get this scaling we're getting calculating the energy in
   %a square that is 100% strain by 1 maximum isometric force
   %      first in units of Nm 
   %      second in normalized units
   %And then taking the ratio to obtain the scaling    
    fiberEnergyScaling = fiberStretchAtOneNormForce*fiso ...
                        / fiberStrainAtOneNormForce*1;

    ppe = ppeN*fiberEnergyScaling;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%8. Populate the dynamics and energy information structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtInfo.muscleDynamicsInfo.activation                = a;                  % NA                   NA						
mtInfo.muscleDynamicsInfo.fiberForce                = ffN*fiso;           % force                N
mtInfo.muscleDynamicsInfo.fiberForceAlongTendon     = ffN*fiso*cosAlpha;  % force                N
mtInfo.muscleDynamicsInfo.normFiberForce            = ffN;            % force/force          N/N
mtInfo.muscleDynamicsInfo.activeFiberForce          = ffaN*fiso;      % force                N
mtInfo.muscleDynamicsInfo.passiveFiberForce         = ffpN*fiso;      % force                N
mtInfo.muscleDynamicsInfo.tendonForce               = ftN*fiso;        % force                N
mtInfo.muscleDynamicsInfo.normTendonForce           = ftN;             % force/force          N/N
mtInfo.muscleDynamicsInfo.fiberStiffness            = kf;             % force/length         N/m
mtInfo.muscleDynamicsInfo.fiberStiffnessAlongTendon = kfAT;           % force/length         N/m
mtInfo.muscleDynamicsInfo.tendonStiffness           = kt;             % force/length         N/m
mtInfo.muscleDynamicsInfo.muscleStiffness           = km;             % force/length         N/m
mtInfo.muscleDynamicsInfo.fiberActivePower          = -ffaN*fiso*dlce; % force*velocity       W
mtInfo.muscleDynamicsInfo.fiberPassivePower         = -ffpN*fiso*dlce; % force*velocity       W
mtInfo.muscleDynamicsInfo.tendonPower               = -ftN*fiso*dlt;    % force*velocity       W
mtInfo.muscleDynamicsInfo.musclePower               = -ffN*fiso*dlce;  % force*velocity       W

%%
%These next 4 fields do not appear in the OpenSim model
%%
mtInfo.muscleDynamicsInfo.fiberParallelElementPower = -fiso*fpeN*dlce; 
mtInfo.muscleDynamicsInfo.dampingForces             =  fiso*beta*dlceN;
mtInfo.muscleDynamicsInfo.dampingPower              = -(fiso*beta*dlceN)*dlce;
mtInfo.muscleDynamicsInfo.boundaryPower             =  fiso*ftN*dlp;

mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy  = ppe;   %force*distance    J   
mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy = pt;    %force*distance    J   
mtInfo.musclePotentialEnergyInfo.musclePotentialEnergy = ppe+pt;%force*distance    J (Nm)

mtInfo.state.value      = NaN;
mtInfo.state.derivative = NaN;

%This field does not appear in the OpenSim model
mtInfo.initialization = init;

if useElasticTendon == 1
    mtInfo.state.value      = lceAT;
    mtInfo.state.derivative = dlceAT;
end