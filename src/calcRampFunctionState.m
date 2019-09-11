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
function rampState = calcRampFunctionState(t,initialPauseTime,...
                                            rampStartLength,rampEndLength, ...
                                            rampTime)
%%
% This function evaluates the value and derivative of a ramp function 
% that takes this form:
%
%  
%
%       0    1
%      |*----* 
%      |      \
% y(t) |       \
%      |        \  
%      |_________*2________________________ 
%                      time (s)
%
% @param t: time
% @param initialPauseTime: time between 0 and 1
% @param rampStartLength: length of the function times 0 and 1
% @param rampEndLength  : length of the function at time 2
% @param rampTime       : elapsed time between 1 & 2
%
% @returns rampState 2x1 vector:
%          rampState(1) = dy(t)/dt 
%          rampState(2) =  y(t)
%
%%

rampState = zeros(2,1);
dydt = 0;
y    = rampStartLength;

if(t > initialPauseTime && t < (initialPauseTime+rampTime))
  dydt = (rampEndLength - rampStartLength)/rampTime;
  y    =  rampStartLength ...
        +((t-initialPauseTime)/(rampTime))*(rampEndLength-rampStartLength);
end
if(t>(initialPauseTime+rampTime))
  dydt = 0;
  y    = rampEndLength;
end

rampState(1) = dydt;
rampState(2) = y;
