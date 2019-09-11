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
function fiberState = clampFiberState(lce,dlce,lceMin)
%%
% This function returns the state of the fiber, and if necessary, 
% clamps it to a lower bound defined by lceMin
%
% @param lce: length of the fiber (m)
% @param dlce: fiber lengthening velocity (m/s)
% @param lceMin: the minimum allowable length of the fiber
%
% @return fiberState a structure with fields of
%           .lce : fiber length
%           .dlce: fiber velocity
%           .isClamped: 0 - not clamped
%                       1 - clamped
%%
isClamped = 0;
if(lce < lceMin || (lce == lceMin && dlce < 0))            
    %Clamp the fiber length along the tendon
    lce  = lceMin;
    dlce = 0;    
    isClamped = 1;
end

fiberState.lce     = lce;
fiberState.dlce    = dlce;
fiberState.isClamped = isClamped; 