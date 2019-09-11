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
function curveInv = createInverseBezierCurve(curve)

curveInv.xpts = curve.ypts;
curveInv.ypts = curve.xpts;

curveInv.xEnd = curve.yEnd;
curveInv.yEnd = curve.xEnd;

col  = size(curve.xpts,2);
row  = size(curve.xpts,1);
n    = row-1; %polynomial order

dx1 = n.*(curveInv.xpts(2:row,:)-curveInv.xpts(1:(row-1),:));
dy1 = n.*(curveInv.ypts(2:row,:)-curveInv.ypts(1:(row-1),:));

dx2 = (n-1).*(dx1(2:(row-1),:)-dx1(1:(row-2),:));
dy2 = (n-1).*(dy1(2:(row-1),:)-dy1(1:(row-2),:));



dydx0 = dy1(1,1)/dx1(1,1);
dydx1 = dy1(row-1,col)/dx1(row-1,col);

d2ydx20 = 0;
d2ydx21 = 0;

if(curve.d2ydx2End(1) ~= 0)
    d2ydx20 = (dy2(1,1)*dx1(1,1)...
             - dy1(1,1)*dx2(1,1))...
              /( dx1(1,1)^2 );
end

if(curve.d2ydx2End(2) ~= 0)
    d2ydx21 = (dy2(row-2,col)*dx1(row-1,col)...
             - dy1(row-1,col)*dx2(row-2,col))...
               /( dx1(row-1,col)^2 );
end

curveInv.dydxEnd   = [  dydx0,  dydx1];
curveInv.d2ydx2End = [d2ydx20, d2ydx21];

%I'm leaving this empty for now, as the scaling information
%in the y direction is not embedded in the structure, and
%don't yet have a use case for the integral of the inverse
%function
curveInv.integral = [];

curveInv.name = [curve.name,'INV'];