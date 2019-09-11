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
function fiberState = clampFiberStateAlongTendon(lceAT,dlceAT,lceATMin)
%%
%
% @param lceAT: fiber length along the tendon
% @param dlceAT: d/dt lceAt - fiber velocity along the tendon
% @param lceATMin: minimum fiber length permitted along the tendon
%
% @returns fiberState, a struct with fields of
%
%       .lceAT : length of the fiber along the tendon
%       .dlceAT: fiber velocity along the tendon
%       .isClamped: 0 if the fiber value has not been clamped
%                   1 if the fiber velocity has been clamped
%%
isClamped = 0;
if(lceAT <= lceATMin || (lceAT == lceATMin && dlceAT < 0))            
    %Clamp the fiber length along the tendon
    lceAT  = lceATMin;
    dlceAT = 0;    
    isClamped = 1;
end

fiberState.lceAT     = lceAT;
fiberState.dlceAT    = dlceAT;
fiberState.isClamped = isClamped; 