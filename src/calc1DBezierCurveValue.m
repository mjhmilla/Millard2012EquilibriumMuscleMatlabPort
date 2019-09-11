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
function f = calc1DBezierCurveValue(u, pV)
%%
% This function implements De Casteljau's recursive algorithm for
% evaluating an nth order 1D Bezier curve of the form
%
%  B(u) = sum_{i=0}^n  [(n choose i)(1-u)^{n-1} u^i] p_i
%
% where 
%  u: argument of the curve
%  n: order of the curve
%  p_i: value of the ith control point
%
% For an n th order curve with (n+1) points this algorithm requires n! 
% subtractions, multiplications, and additions to terminate, and a stack 
% that is n! deep. Although this algorithm is very general faster results 
% can be obtained using optimized code if the order of the Bezier curve is 
% known ahead of time using optimized code.
%
% @params u : [0,1] the argument of the Bezier curve
% @params pV: vector of control points
% @returns f: the value of the Bezier curve evaluated at u
%
% Reference
% [1] http://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm
%
%%
bV0 = zeros(size(pV,1)-1,1);


for i=1:1:length(bV0)
    bV0(i) = (pV(i+1)-pV(i))*(u) + pV(i);
end

if(length(bV0) == 1)    
    f  = bV0;
else
    f = calc1DBezierCurveValue(u,bV0);
end


