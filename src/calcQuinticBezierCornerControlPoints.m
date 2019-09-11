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
function xyPts = calcQuinticBezierCornerControlPoints(x0,y0,dydx0,d2ydx20,...
                                              x1,y1,dydx1,d2ydx21,...
                                              curviness)
%%
% Note: This is an improved version of the code that is in OpenSim: this
%       function allows you to set the 2nd derivative at the end points of
%       the Bezier curve section. The code in OpenSim automatically sets
%       these endpoints to 0.
%
%        Calculates the location of quintic Bezier curve control points to 
%        create a C shaped curve.
%
%        @param x0       First intercept x location
%        @param y0       First intercept y location
%        @param dydx0    First intercept slope
%        @param x1       Second intercept x location
%        @param y1       Second intercept y location
%        @param dydx1    Second intercept slope
%        @param curviness A parameter that ranges between 0 and 1 to denote a 
%                         straight line or a curve
%        @throws OpenSim::Exception 
%         -If the curviness parameter is less than 0, or greater than 1;
%         -If the points and slopes are chosen so that an "S" shaped curve would 
%          be produced. This is tested by examining the points (x0,y0) and 
%          (x1,y1) together with the intersection (xC,yC) of the lines beginning 
%          at these points with slopes of dydx0 and dydx1 form a triangle. If the 
%          line segment from (x0,y0) to (x1,y1) is not the longest line segment, 
%          an exception is thrown. This is an overly conservative test as it 
%          prevents very deep 'V' shapes from being respresented.
%
%        @return a SimTK::Matrix of 6 points Matrix(6,2) that correspond to the 
%                         X, and Y control points for a quintic Bezier curve that
%                         has the above properties
%
%
%        Calculates the location of quintic Bezier curve control points to 
%        create a C shaped curve that intersects points 0 (x0, y0) and point 1
%        (x1, y1) with slopes dydx0 and dydx1 respectively, and a second 
%        derivative of 0. The curve that results can approximate a line 
%        (curviness = 0), or in a smooth C shaped curve (curviniess = 1)
%
%        The current implementation of this function is not optimized in anyway
%        and has the following costs:
%
%        <B>Computational Costs</B>
%        \verbatim
%            ~55 flops
%        \endverbatim
%
%
%        <B>Example:</B>
%            @code
%            double x0 = 1;
%            double y0 = 0;
%            double dydx0 = 0;
%            double x1 = 1.04;
%            double y1 = 1;
%            double dydx1 = 43;
%            double c = 0.75;
%
%            SimTK::Matrix p0 = SegmentedQuinticBezierToolkit::
%               calcQuinticBezierCornerControlPoints(x0, y0, dydx0,x1,y1,dydx01,
%                                                                     c);
%            @endcode
%
%        */
%%
                                              
xyPts = zeros(6,2);
assert( (curviness >= 0 && curviness <= 1),...
        'Error: curviness must be [0 1]');       

 
    
xScl = 1/(x1-x0);    
%1. Calculate the location where the two lines intersect
% (x-x0)*dydx0 + y0 = (x-x1)*dydx1 + y1
%   x*(dydx0-dydx1) = y1-y0-x1*dydx1+x0*dydx0
%                 x = (y1-y0-x1*dydx1+x0*dydx0)/(dydx0-dydx1);

xC = 0;
yC = 0;
rootEPS = eps^0.5;
if(abs(dydx0-dydx1) > rootEPS)
    xC = (y1-y0-x1*dydx1+x0*dydx0)/(dydx0-dydx1);    
else
    xC = (x1+x0)/2;
end

yC = (xC-x1)*dydx1 + y1;
%Check to make sure that the inputs are consistent with a corner, and will
%not produce an 's' shaped section. To check this we compute the sides of
%a triangle that is formed by the two points that the user entered, and 
%also the intersection of the 2 lines the user entered. If the distance
%between the two points the user entered is larger than the distance from
%either point to the intersection loctaion, this function will generate a
%'C' shaped curve. If this is not true, an 'S' shaped curve will result, 
%and this function should not be used.

xCx0 = (xC-x0);
yCy0 = (yC-y0);
xCx1 = (xC-x1);
yCy1 = (yC-y1);
x0x1 = (x1-x0);
y0y1 = (y1-y0);

a = xCx0*xCx0 + yCy0*yCy0;
b = xCx1*xCx1 + yCy1*yCy1;
c = x0x1*x0x1 + y0y1*y0y1;

%This error message needs to be better.
%assert( ((c > a) && (c > b)), ...
%    'The endpoints and slopes do not intersect.');

%Start point
xyPts(1,1) = x0;
xyPts(1,2) = y0;
%End point
xyPts(6,1) = x1;
xyPts(6,2) = y1;


%Original code - leads to 2 localized corners
xyPts(2,1) = x0 + curviness*(xC-xyPts(1,1));
xyPts(2,2) = y0 + curviness*(yC-xyPts(1,2));



%d2ydx20 = 0; %later this will be an input

dxdu0   = 5*(xyPts(2,1)-xyPts(1,1));
dydu0   = 5*(xyPts(2,2)-xyPts(1,2));

xyPts(3,1) = xyPts(2,1) + 0.5*(xC-xyPts(2,1));
d2xdu20 = 20*(xyPts(3,1) - 2*xyPts(2,1) + xyPts(1,1));



d2ydu20 = (dxdu0*dxdu0*(d2ydx20) + d2xdu20*(dydx0));

xyPts(3,2) = d2ydu20*(1/20) + 2*xyPts(2,2) - xyPts(1,2) ;
%chkX = 20 *( xyPts(3,1)-2*xyPts(2,1)+xyPts(1,1));
%chkY = 20 *( xyPts(3,2)-2*xyPts(2,2)+xyPts(1,2));
%chkdydx2 = (d2ydu20*(1/dxdu0) - dydu0*d2xdu20/(dxdu0*dxdu0))/dxdu0;


xyPts(5,1) = xyPts(6,1) + curviness*(xC-xyPts(6,1));
xyPts(5,2) = xyPts(6,2) + curviness*(yC-xyPts(6,2));

%d2ydx21 = 0; %later this will be an input
dxdu1 = 5*(xyPts(6,1)-xyPts(5,1));

xyPts(4,1) = xyPts(5,1) + 0.5*(xC-xyPts(5,1));
d2xdu21 = 20*(xyPts(4,1) - 2*xyPts(5,1) + xyPts(6,1) );

d2ydu21 = (dxdu1*dxdu1*(d2ydx21) + d2xdu21*(dydx1));


xyPts(4,2) = d2ydu21*(1/20) + 2*xyPts(5,2) - xyPts(6,2) ;

    
        
                                              

                                              