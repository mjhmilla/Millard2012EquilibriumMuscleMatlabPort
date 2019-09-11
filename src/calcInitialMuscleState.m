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
function initialSolution = calcInitialMuscleState(activationState,...
                                                 pathState,...
                                                 muscleArchitecture,...
                                                 calcMuscleInfo,...
                                                 config)
%%
% This function calculates an initial fiber length for an elastic tendon
% model that satisfies the force equilibrium equation 
%
%   fiberForceAlongTendon - tendonForce = 0
% 
% of the model. The crudest form of initialization assumes that the fiber
% lengthening velocity is zero, that is that the tendon is doing all of the
% lengthening. This assumption is usually poor given that the tendon is
% generally far stiffer than the muscle fiber. That said this crude
% approach is guaranteed to converge.
%
% If the model is configured to do so the initial fiber length and velocity
% are solved for such that the equlibrium equation and an approximation to 
% its first time derivative are satisfied
%
%       fiberForceAlongTendon -      tendonForce = 0  [1]
%  d/dt fiberForceAlongTendon - d/dt tendonForce = 0. [2]
%
% Here the word approximation is used since we often don't know d/dt
% activation, nor do we know d^2 lce/ dt^2 and so these quantites are
% assumed to be zero. Even though this is an approximation the resulting
% initialization produces initial forces and fiber lengths that do not have
% large transients in the beginning of the simulation.
%
%@param activationState : scalar activation of the muscle
%
%@param pathState : 2x1 vector of the path lengthening rate and path length
%                   in m/s and m respectively.
%
%@param muscleArchtecture: a structure containing the architecture
%        information about the muscle. This structure has fields of:
%
%
%     .name    : string name of the muscle. e.g. 'soleus left'                      
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
% @param calcMuscleInfoFcn: a function handle to a muscle function
%                      like (calcMillard2012DampedEquilibriumMuscleInfo.m) 
%                      and takes arguments of activation, pathState, and 
%                      muscle state
%
%@param config a structure that configures the initialization routine
%
%         config.iterMax : maximum number of Newton steps to attempt
%         config.tol     : desired normalized force tolerance to satisfy
%         config.useStaticFiberSolution: 0 : assumes fiber velocity is 0
%                                        1 : estimates fiber velocity using
%                                            by satisfying an approximation
%                                            to the first derivative of the
%                                            equilibrium equation Eqn. 2.
%
%@returns initialSolution a structure with the fields
%
%         .muscleState : starting state of the muscle
%         .converged   : 0 if the initialization attempt failed
%                        1 if the initialization attempt succeeded
%         .isClamped   : 0 if the fiber is in its operational range
%                        1 if the fiber has been clamped to its lower bound
%
%%
 
dlceAT= NaN;
lceAT = NaN;

initialSolution = [];
muscleState = [];
            
epsRoot = eps^0.5;
a   = activationState(1);
dlp = pathState(1);
lp  = pathState(2);

lceATMin = muscleArchitecture.minimumFiberLengthAlongTendon;
ltSlk    = muscleArchitecture.tendonSlackLength;

lceOpt    =  muscleArchitecture.optimalFiberLength;
alphaOpt  =  muscleArchitecture.pennationAngle;
h         = lceOpt*sin(alphaOpt);

isClamped = 0;
converged = 0;

if(lp-ltSlk < lceATMin)
    isClamped = 1;    
end

if(isClamped == 0)
    iterMax                = config.iterMax;
    tol                    = config.tol;
    useStaticFiberSolution = config.useStaticFiberSolution;
    
    %If the path isn't moving to numerical tolerance use the
    %static solution. Otherwise you can get the fiber and 
    %tendon doing opposite things ... because the fiber can
    %have (unphysical) negative stiffness. Time for a new 
    %approach!    
    if(abs(dlp) < epsRoot)
       useStaticFiberSolution = 1; 
    end
    
    lceAT   = lp - 1.01*ltSlk;
    dlceAT  = 0;
        
    iter    = 1;
    err     = 2*tol;
        
    while abs(err) > tol && iter < iterMax
        muscleState = [dlceAT;lceAT];
        
        mtInfo = calcMuscleInfo(activationState,...
                                pathState, ...
                                muscleState);
        
        err         = mtInfo.initialization.err;
        Derr_DlceAT = mtInfo.initialization.Derr_DlceAT;                  
        
        if(abs(err) > tol && abs(Derr_DlceAT) > eps*10)
            delta = -err/Derr_DlceAT;
            lceAT = lceAT + delta;            
        elseif( abs(err) > tol && abs(Derr_DlceAT) <= eps*10)
               %perturb the solution and keep trying
               randStep = rand(1)*lceOpt; 
               randSign = rand(1);
               if(randSign < 0.5)
                  randSign = -1; 
               else
                  randSign =  1;
               end
               lceAT = lceAT + randSign*randStep;
               if(lceAT < lceATMin)
                  lceAT = lceATMin; 
               end          
        end
        
        if(useStaticFiberSolution == 1)
           dlceAT = 0; 
        else
           kt   = mtInfo.muscleDynamicsInfo.tendonStiffness;
           kfAT = mtInfo.muscleDynamicsInfo.fiberStiffnessAlongTendon;
           
           %Here we split the velocity of the path between the
           %fiber and the tendon by taking the derivative of the
           %force eqilibrium equation ()
           %
           % d/dt [(a*fal*fv + fpe + beta*dlceN)cos(alpha) - ft] = 0 
           %
           % assuming that da/dt = 0 and d2lce/dt2 = 0 we get
           %  (a*DfalDlce*fv + DfpeDlce)cos(alpha)*DlceDt
           % -(a*fal*fv + fpe + beta*dlceN)*sin(alpha)*DalphaDlce*DlceDt
           % - DftDlt*DltDt = 0
           %
           % which can be simplified to           
           %  kfAT*DlceDt - kt*DltDt = 0
           %
           % Since
           %   DlceDt + DltDt = DlpDt 
           %   
           % We have
           %   kfAT*DlceDt - kt*(DlpDt-DlceDt) = 0 
           %
           %   DlceDt = kt*DlpDt / (kfAT + kt)
           %
           % And since kfAT can be negative, this can blow up.
           % Fun.
           if(abs(kfAT + kt) > epsRoot)
              dlce = kt*dlp / (kt + kfAT); 
           else
              dlce = 0; 
           end
           
           %Now put this along the tendon
           lce    = sqrt(lceAT*lceAT + h*h);
           fiberKinematics =...
                calcFixedWidthPennatedFiberKinematicsAlongTendon(lce,...
                                                                 dlce,...
                                                                 lceOpt,...
                                                                 alphaOpt);
           dlceAT = fiberKinematics.fiberVelocityAlongTendon;
           
        end
        
        iter = iter + 1;
    end
    
    if abs(err) < tol
       converged = 1; 
    end
        
        
    fiberState = clampFiberStateAlongTendon(lceAT,dlceAT,lceATMin);
    isClamped  = fiberState.isClamped;
end

if(isClamped == 1)
   lceAT  = lceATMin;
   dlceAT = 0;
end

initialSolution.muscleState = [lceAT];
initialSolution.converged   = converged;
initialSolution.isClamped   = isClamped;
