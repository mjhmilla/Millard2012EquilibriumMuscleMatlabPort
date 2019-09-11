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
function integralStruct = createCurveIntegralStructure(curveParams, npts, tol, xScaling)
%%
% Numerically evaluates the integral at npts over the domain of the 2D 
% Bezier curve defined by the matrix of x control points and y control.
% The first and second derivative of the curve integral are also evaluated
% at these points so that it is possible to interpolate the curve using a
% quintic hermine spline.
%
% @param x: value to 
% @param xpts: n x m matrix of Bezier control points defining x(t). 
%
%              n x m matrix 
%              n: Bezier curve order + 1 
%              m: # spline sections
%                
%              Each column defines the control points used for a Bezier 
%              section. Thus a quintic Bezier spline with 3 sections will 
%              have a 6 x 3 matrix for xpts
%              
% @param ypts: n x m matrix of Bezier control points defining y(t)
%
% @param npts: number of intermediate points between xmin and xmax to
%              evaluate the integeral.
% @param tol: relative and absolute tolerance on the integral.
%
% @return integralStruct: A structure containing the numerically calculated
%                         integral, its first and second derivative so that 
%                         the integral curve can be interpolated using a 
%                         quintic Hermite spline.
%
%           Has fields of
%                      .xptsN : points the integral is evaluated at  
%                      .yptsN : numerical value of the integral
%                      .y1ptsN: '' integral's first derivative.
%                      .y2ptsN: '' integral's 2nd derivative.
%
%%

fcn = @(arg,arg1)calcBezierYFcnXDerivative(arg, curveParams, 0);
xmin = min(min(curveParams.xpts));
xmax = max(max(curveParams.xpts));
xv   = [xmin:((xmax-xmin)/(npts-1)):xmax];


options = odeset('RelTol',tol,'AbsTol',tol);
[xe ye] = ode45(fcn,xv,0,options);

ye1 = zeros(size(ye));
ye2 = zeros(size(ye));

for i=1:1:length(ye)
   ye1(i) =  calcBezierYFcnXDerivative(xe(i), curveParams, 0);
   ye2(i) =  calcBezierYFcnXDerivative(xe(i), curveParams, 1);
end

integralStruct.xptsN  = xe;
integralStruct.yptsN  = ye;
integralStruct.y1ptsN = ye1;
integralStruct.y2ptsN = ye2;
integralStruct.xScaling = xScaling;

