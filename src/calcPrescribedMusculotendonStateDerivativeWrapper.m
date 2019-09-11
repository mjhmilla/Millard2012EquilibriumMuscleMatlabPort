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
function dState = ...
    calcPrescribedMusculotendonStateDerivativeWrapper(t,...
                                                      state,...                                                      
                                                      prescribedPathFcn,...
                                                      prescribedActivationFcn,...
                                                      calcMuscleInfoFcn)
%%
% This is a wrapper ultimately to create a derivative function that takes t
% and muscle state as arguments and returns the muscle state derivative,
% and the various powers of interest, as the muscle undergoes a constant
% activation sinusoidal stretch benchmark simulation
%
% @param t: time in seconds
%
% @param state: the state vector for this benchmark simulation. This vector
%               contains:
%
%               [muscleState;
%                boundaryPower;
%                activeFiberPower;
%                dampingPower]
% 
%                The last three entries are required to numerically
%                evaluate T + V - W = const to ensure that the model is
%                conservative.
%
% @param prescribedPathFcn : a handle to a function that given time t
%                            produces a 2x1 vector containing the velocity 
%                            and length of the path the muscle lies on.
%
% @param prescribedActivationFcn: a handle to a function that given time t
%                            produces a 1x1 scalar of the activation of the
%                            muscle.
%
% @param calcMuscleInfoFcn: a function handle to a muscle function
%                      like (calcMillard2012DampedEquilibriumMuscleInfo.m) 
%                      and takes arguments of activation, pathState, and 
%                      muscle state
%
% @returns dState: the first time derivative of the state vector
%
%               [muscleState;
%                boundaryPower;
%                activeFiberPower;
%                dampingPower]
%
%
%%
                                                  
pathState       = prescribedPathFcn(t);
activationState = prescribedActivationFcn(t);

mtInfo = calcMuscleInfoFcn(activationState,...
                           pathState,...
                           state);      

dlp = pathState(1);                       

boundaryPower       = mtInfo.muscleDynamicsInfo.boundaryPower;
activeFiberPower    = mtInfo.muscleDynamicsInfo.fiberActivePower;
dampingPower        = mtInfo.muscleDynamicsInfo.dampingPower;


dState = [];

if(length(state) > 3)
    dState = [mtInfo.state.derivative; ...
               boundaryPower;...
               activeFiberPower;...
               dampingPower];    
else
    dState = [ boundaryPower;...
               activeFiberPower;...
               dampingPower]; 
end

