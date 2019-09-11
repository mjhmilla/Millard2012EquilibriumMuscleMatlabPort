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
function curveValues = calcBezierYFcnXCurveSampleVector(curveParams, npts)
%%
% This function evaluates a Bezier spline curve (x(u), y(u)) across its
% entire domain, and and slightly beyond to get values in the extrapolated
% region. The spline is evaluated at its value, and first 3 derivatives. In
% addition if the curve has its integral function defined, then the integral
% is also evaluated.
%
% @param curveParams : a structure with the fields
%           .xpts
%           .ypts
%           .integral
%
%   xpts: n x m matrix of Bezier control points defining x(t). 
%
%              n x m matrix 
%              n: Bezier curve order + 1 
%              m: # spline sections
%                
%              Each column defines the control points used for a Bezier 
%              section. Thus a quintic Bezier spline with 3 sections will 
%              have a 6 x 3 matrix for xpts
%              
%   ypts: n x m matrix of Bezier control points defining y(t)
%
%   integral: if the integral curve has not been numerically computed, 
%             then this field will be empty. Otherwise it will have fields
%             of
%
%                      .xptsN : points the integral is evaluated at  
%                      .yptsN : numerical value of the integral
%                      .y1ptsN: '' integral's first derivative.
%                      .y2ptsN: '' integral's 2nd derivative.
%
% @param npts: the number of samples to use across the curve domain
%              
% @return curveValues, a struct with the fields
%
%	.x      : vector of x values of the Bezier spline
%	.y      : vector of y values of the Bezier spline
%	.dydx   : first derivative        dy/dx
%	.d2ydx2 : second derivative       d^2y/dx^2
%	.d3ydx3 : third derivative        d^3y/dx^3
%	.intYdx : integral of y w.r.t. x  dy/dx
%%

xmin  = min(min(curveParams.xpts));
xmax  = max(max(curveParams.xpts));
delta = xmax-xmin;
xmin  = xmin-delta/5;
xmax  = xmax+delta/5;

ymin  = min(min(curveParams.ypts));
ymax  = max(max(curveParams.ypts));
delta = ymax-ymin;
ymin  = ymin-delta/5;
ymax  = ymax+delta/5;

x = [xmin:((xmax-xmin)/npts):xmax]';  
    
y      = zeros(size(x));
dydx   = zeros(size(x));
d2ydx2 = zeros(size(x));
d3ydx3 = zeros(size(x));
intYdx = [];


if(isempty(curveParams.integral) == 0)
   intYdx = zeros(size(x)); 
end

intYdx = zeros(size(x));

for k=1:1:length(x)
   if(isempty(curveParams.integral) == 0)
     intYdx(k) = calcBezierYFcnXDerivative(x(k), curveParams, -1); 
   end 
    
   y(k)      = calcBezierYFcnXDerivative(x(k), curveParams, 0);
   dydx(k)   = calcBezierYFcnXDerivative(x(k), curveParams, 1);
   d2ydx2(k) = calcBezierYFcnXDerivative(x(k), curveParams, 2);
   d3ydx3(k) = calcBezierYFcnXDerivative(x(k), curveParams, 3);
end

curveValues.x      = x;
curveValues.y      = y;
curveValues.dydx   = dydx;
curveValues.d2ydx2 = d2ydx2;
curveValues.d3ydx3 = d3ydx3;
curveValues.intYdx = intYdx;
