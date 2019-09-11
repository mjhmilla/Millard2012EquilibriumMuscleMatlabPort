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
function mtInfo = calcElasticTendonDampedEquilibriumMuscleInfo(...
                                            activationState,...                                            
                                            muscleState,...
                                            pathState, ...                                            
                                            muscleArchitecture,...
                                            normMuscleCurves,...
                                            config)


                                        
mtInfo = [];

a = activationState(1);

lceAT = muscleState(1);

dlp = pathState(1);
lp  = pathState(2);

fiso      =  muscleArchitecture.fiso;
lceOpt    =  muscleArchitecture.optimalFiberLength;
alphaOpt  =  muscleArchitecture.pennationAngle;
ltSlk     =  muscleArchitecture.tendonSlackLength;
dlceMaxN  =  muscleArchitecture.maximumNormalizedFiberVelocity;

flNcurve  = normMuscleCurves.activeForceLengthCurve;
fvNcurve  = normMuscleCurves.fiberForceVelocityCurve;
fpeNcurve = normMuscleCurves.fiberForceLengthCurve;
ftNcurve  = normMuscleCurves.tendonForceLengthCurve;

flNcurveHACK    = normMuscleCurves.activeForceLengthCurveHACK;
fvInvNcurveHACK = normMuscleCurves.fiberForceVelocityInverseCurveHACK;
fvNcurveHACK    = normMuscleCurves.fiberForceVelocityCurveHACK;

lceMin       = muscleArchitecture.minimumFiberLength;
lceATMin     = muscleArchitecture.minimumFiberLengthAlongTendon;
alphaLceMin  = muscleArchitecture.pennationAngleAtMinumumFiberLength;

%%
%Set up function handles for the active force length
%and force velocity curves, as these curves are different
%depending on whether fiber damping is used or not.
%%

falNFcn   = [];
fvNFcn    = [];
fvNInvFcn = [];

if(config.damping > 0)
   falNFcn   = @(arg1, arg2)calcBezierYFcnXDerivative(arg1,  flNcurve,  arg2);
   fvNFcn    = @(arg1, arg2)calcBezierYFcnXDerivative(arg1,  fvNcurve,  arg2);
   fvNInvFcn = [];
else
   falNFcn   = @(arg1,arg2)calcBezierYFcnXDerivative(arg1,   flNcurveHACK,  arg2);
   fvNFcn   = @(arg1, arg2)calcBezierYFcnXDerivative(arg1,  fvNcurveHACK,  arg2);
   fvNInvFcn = @(arg1,arg2)calcBezierYFcnXDerivative(arg1,fvInvNcurveHACK,  arg2); 
end
  


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

lce  = NaN;
dlce = NaN;

alpha= NaN;
dalpha=NaN;

lt    =NaN;
dlt   =NaN;

%Clamp the fiber length at its lower bound if necessary
if(lceAT <= lceATMin)
    
    lceAT = lceATMin;
    lce   = lceMin;
    alpha = alphaLceMin;
    lt    = lp - lceAT;  %The tendon will be shorter than ltSlk here
        
    dlceAT  = 0;
    dlce    = 0;
    dalpha = 0;
    dlt     = dlp - dlceAT; %Note that the tendon is buckling here so although its
                         %length isn't changing it is getting shorter    
    flag_fiberClamped = 1;
end
   
%%
%Evaluate the length dependent muscle curves
%%
lceN   = lce/lceOpt;
ltN    = lt/ltSlk;
etN    = ltN-1;

falN = falNFcn(lceN,0);
fpeN = calcBezierYFcnXDerivative(lceN, fpeNcurve, 0);
ftN  = calcBezierYFcnXDerivative(ltN,  ftNcurve,  0);
    
if(flag_fiberClamped == 0 )    
    %calculate the fiber length
    fiberKinematics = calcFixedWidthPennatedFiberKinematics(lceAT,...
                                                            NaN,...
                                                            lceOpt,...
                                                            alphaOpt);
    lce    = fiberKinematics.fiberLength;
    lceN   = lce/lceOpt;
    alpha  = fiberKinematics.pennationAngle;

    lt = lp-lceAT;
    ltN= lt/ltSlk;
    
    %if there is damping iterate over the equilibrium equation
    %to get the fiber state, else invert the force velocity curve.
    
    if(config.damping > 0)
       beta = config.damping; 
        
       
    else
       %Compute fiber velocity by inverting the force-velocity curve 
       %How? Here's the force equilibrium equation
       % 
       % (a*fal*fv + fpe)*cos(alpha) - ft = 0;
       %  dlceN = fvInv( [ft/cos(alpha) - fpe]/(a*fal) )
       %
       % Will go singular as
       %  cos(alpha) -> 0
       %  a -> 0
       %  fal -> 0
       %
       % And so before proceeding we assert that the denominators 
       % cannot be 0.
       assert(cos(alpha) > eps^0.5,...
           'Cannot invert fv curve: cos(pennationAngle) <= eps^0.5');
       assert(a     > eps^0.5,...
           'Cannot invert fv curve: activation <= eps^0.5');
       assert(falN > eps^0.5,...
           'Cannot invert fv curve: forceActiveLengthMultiplier <= eps^0.5');

       fvN   = ((ftN/cos(alpha)) - fpeN) / (a*falN);
       dlceN = fvNInvFcn(fvN,0);       
       dlce  = dlceN*(lceOpt*dlceMaxN);
                    
    end
    
    %%
    %Calculate the pennation angular velocity and the rate of stretch
    %of the tendon
    %%
    fiberKinAlongTendon = ...
        calcFixedWidthPennatedFiberKinematicsAlongTendon(lce,...
                                                     dlce,...
                                                     lceOpt,...
                                                     alpha);
                                                 
    lceAT  = fiberKinAlongTendon.fiberLengthAlongTendon;                                      
    dlceAT = fiberKinAlongTendon.fiberVelocityAlongTendon;
    dalpha = fiberKinAlongTendon.pennationAngularVelocity;
        
    dlt    = dlp-dlceAT;
    
end

%%
%2. Normalize the fiber length and velocity, evaluate velocity dependent
%   curves
%%
dlceN  = dlce/(lceOpt*dlceMaxN);
dltN   = dlt/ltSlk;

fvN    = fvNFcn(dlceN,0); 

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
mtInfo.muscleLengthInfo.cosPennationAngle            = cos(alpha); %NA                NA         
mtInfo.muscleLengthInfo.sinPennationAngle            = sin(alpha); %NA                NA         	
mtInfo.muscleLengthInfo.fiberPassiveForceLengthMultiplier  = fpeN; %NA             NA
mtInfo.muscleLengthInfo.fiberActiveForceLengthMultiplier   = falN; %NA             NA
          
mtInfo.fiberVelocityInfo.fiberVelocity                = dlce;    %length/time           m/s
mtInfo.fiberVelocityInfo.fiberVelocityAlongTendon     = dlceAT;  %length/time           m/s
mtInfo.fiberVelocityInfo.normFiberVelocity            = dlceN;   %(length/time)/length  (m/s)/m
mtInfo.fiberVelocityInfo.pennationAngularVelocity     = dalpha;  %angle/time            rad/s
mtInfo.fiberVelocityInfo.tendonVelocity               = dlt;     %length/time           m/s
mtInfo.fiberVelocityInfo.normTendonVelocity           = dltN;    %(length/time)/length  (m/s)/m
mtInfo.fiberVelocityInfo.fiberForceVelocityMultiplier = fvN;     %force/force           NA



%%
%5. Calculate the forces
%%

fiberForceAlongTendonNorm =  (a*falN*fvN + fpeN)*cos(alpha);
fiberForceNorm            =  (a*falN*fvN + fpeN); %(f)orce (f)iber (n)ormalized
activeFiberForceNorm      =   a*falN*fvN;
passiveFiberForceNorm     =                fpeN;


fiberForceAlongTendon = fiso * fiberForceAlongTendonNorm;
fiberForce            = fiso * fiberForceNorm; %(f)orce (f)iber (n)ormalized
activeFiberForce      = fiso * activeFiberForceNorm;
passiveFiberForce     = fiso * passiveFiberForceNorm;

%Tendon force is equal and opposite unless the tendon is buckling, then
%it is simply zero.
ftN = fiberForceAlongTendonNorm;
ft  = fiberForceAlongTendon;
if(flag_fiberClamped == 1)
    ft = 0;    
end


%%
%6. Calculate the stiffness
%%
dfalDlceN =  calcBezierYFcnXDerivative(lceN,  flNcurve,  1);
dfpeDlceN =  calcBezierYFcnXDerivative(lceN,  fpeNcurve, 1);
dfalDlce  =  dfalDlceN*(1/lceOpt);
dfpeDlce  =  dfpeDlceN*(1/lceOpt);

dAlphaDlce = calcFixedWidthPennationDalphaDlce(alpha,...
                                               lce,...
                                               lceOpt,...
                                               alphaOpt);


kFiber           = fiso*(a*dfalDlce*fvN + dfpeDlce);
kFiberAT         = kFiber*cos(alpha) - fiberForce*sin(alpha)*dAlphaDlce;

kTendon = inf;
if(flag_fiberClamped)
    kTendon = 0;
end

%%
%7. Evaluate the potential energy stored in the musculotendon
%%
pt  = 0; 
ppe = calcBezierYFcnXDerivative(lceN,  fpeNcurve, -1)*(fiso); 
      %Note scaling here is fiso*lceOpt/lceOpt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%8. Populate the dynamics and energy information structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtInfo.muscleDynamicsInfo.activation                = a;                      % NA                   NA						
mtInfo.muscleDynamicsInfo.fiberForce                = fiberForce;             % force                N
mtInfo.muscleDynamicsInfo.fiberForceAlongTendon     = fiberForceAlongTendon;  % force                N
mtInfo.muscleDynamicsInfo.normFiberForce            = fiberForceNorm;         % force/force          N/N
mtInfo.muscleDynamicsInfo.activeFiberForce          = activeFiberForce;       % force                N
mtInfo.muscleDynamicsInfo.passiveFiberForce         = passiveFiberForce;      % force                N
mtInfo.muscleDynamicsInfo.tendonForce               = ft;                     % force                N
mtInfo.muscleDynamicsInfo.normTendonForce           = ftN;                    % force/force          N/N
mtInfo.muscleDynamicsInfo.fiberStiffness            = kFiber;                 % force/length         N/m
mtInfo.muscleDynamicsInfo.fiberStiffnessAlongTendon = kFiberAT;               % force/length         N/m
mtInfo.muscleDynamicsInfo.tendonStiffness           = kTendon;                % force/length         N/m
mtInfo.muscleDynamicsInfo.muscleStiffness           = kFiberAT;               % force/length         N/m
mtInfo.muscleDynamicsInfo.fiberActivePower          = activeFiberForce*dlce;  % force*velocity       W
mtInfo.muscleDynamicsInfo.fiberPassivePower         = passiveFiberForce*dlce; % force*velocity       W
mtInfo.muscleDynamicsInfo.tendonPower               = ft*dlt;                  % force*velocity       W
mtInfo.muscleDynamicsInfo.musclePower               = fiberForce*dlce;        % force*velocity       W

mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy  = ppe;   %force*distance    J   
mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy = pt;    %force*distance    J   
mtInfo.musclePotentialEnergyInfo.musclePotentialEnergy = ppe+pt;%force*distance    J (Nm)

%This muscle is stateless
mtInfo.state.value      = NaN;
mtInfo.state.derivative = NaN;