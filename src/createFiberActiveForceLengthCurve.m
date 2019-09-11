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
function activeForceLengthCurve = createFiberActiveForceLengthCurve(...
            lce0, lce1, lce2, lce3, ...
            minActiveForceLengthValue, plateauSlope, ...
            curviness, computeIntegral, ...
            muscleName)
%%
%       This is a function that will produce a C2 (continuous to the second
%       derivative) active force length curve.
%
%
%       @param lce0   Normalized fiber length at the left-most shoulder of the 
%                     active force-length curve. The value of the active force
%                     length curve for lce < lce0 will be equal to the value
%                     set in shoulderVal. Normally lce0 is approximately 0.5
%       
%       @param lce1   Normalized fiber length at the transition point between 
%                     the ascending limb and the plateau region of the active 
%                     force length curve.
%       
%       @param lce2   Normalized fiber length at the maximum active force length
%                     curve value of 1. Normally lce2 is by definition 1.
%       
%       @param lce3   Normalized fiber length of the at the right most shoulder
%                     of the active-force length curve. The value of the active
%                     force length curve for lce > lce2 will be equal to the 
%                     value of shoulderVal. Normally lce3 is approximately 1.5
%
%       @param minActiveForceLengthValue
%                             The minimum value of the active force length 
%                             curve. A physiological non-equibrium muscle model
%                             would have this value set to 0. An equilibrium 
%                             muscle model would have a non-zero lower bound on 
%                             this value of 0.1 typically. shoulderVal must be 
%                             greater than, or equal to 0.
%                           
%       @param plateauSlope   The slope of the plateau of the active force
%                             length curve between lce1 and lce2. This parameter
%                             can vary depending on the muscle model, but a 
%                             value of 0.8616 is a good place to start.
%
%       @param curviness  The dimensionless 'curviness' parameter that 
%                         can vary between 0 (a line) to 1 (a smooth, but 
%                         sharply bent elbow). A value of 0 will yield an active 
%                         force length curve that is composed of slightly curved 
%                         line segments. A value of 1 will yield an active force
%                         length curve that is smoothly rounded.
%
%       @param computeIntegral If this is true, the integral for this curve
%                              is numerically calculated and splined. If false, 
%                              this integral is not computed, and a call to 
%                              .calcIntegral will throw an exception
%
%       @param muscleName The name of the muscle this curve applies to. This 
%                         curve name should have the name of the muscle and the
%                         curve in it (e.g. "bicep_fiberActiveForceLengthCurve") 
%                         sothat if this curve ever causes an exception, a 
%                         userfriendly error message can be displayed to the
%                         end user to help them debug their model.
%
%       @throws SimTK::Exception if these conditions aren't met
%           -0 < lce0 < lce1 < lce2 < lce3 
%           -shoulderVal >= 0
%           -0 <= plateauSlope < (1/(lce3-lce2))
%           -0 <= curviness <= 1
%
%       @return activeForceLengthCurve 
%               A structure that the function calcNormalizedMuscleCurveDerivative can use to 
%               can use to evaluate the active force length curve value
%               or up to the 3rd derivative.
%
%
%      
%       <B>Conditions:</B>
%
%       <B>Computational Costs</B>
%       \verbatim 
%           Without Integral :   ~20,500 flops
%           With Integral    :  ~870,500 flops
%       \endverbatim
%
%       <B>Example:</B>
%       @code
%           lce0 = 0.5;
%           lce1 = 0.75;
%           lce2 = 1;
%           lce3 = 1.5;
%           shoulderVal  = 0.1;
%           plateauSlope = 0.75;
%           curviness    = 0.9;
%
%           fiberfalCurve = 
%               createFiberActiveForceLengthCurve(lce0, lce1, lce2, lce3, 
%                             shoulderVal, plateauSlope, curviness,false,"test");           
%       @endcode       
%%



activeForceLengthCurve = [];
activeForceLengthCurve.name = sprintf('%s.%s',muscleName,'activeForceLengthCurve');


x0 = lce0;
x1 = lce1;
x2 = lce2;
x3 = lce3;
ylow =  minActiveForceLengthValue;
dydx = plateauSlope;

%%
%Check inputs
%%
rootEPS = eps^0.5;
assert( (x0>=0 && x1>x0+rootEPS  && x2>x1+rootEPS && x3>x2+rootEPS),...    
    sprintf('%s: This must be true: 0 < lce0 < lce1 < lce2 < lce3',...
    activeForceLengthCurve.name));

assert( ylow >= 0,...
    sprintf('%s: shoulderVal must be greater than, or equal to 0',...
    activeForceLengthCurve.name));

dydxUpperBound = (1-ylow)/(x2-x1);
assert(dydx >= 0 && dydx < dydxUpperBound,...    
    sprintf('%s: plateauSlope must be greater than 0 and less than %f',...
    activeForceLengthCurve.name,dydxUpperBound));

assert( (curviness >= 0 && curviness <= 1),...
    sprintf('%s: curviness must be between 0 and 1',...
    activeForceLengthCurve.name));

%%
%Create the curve
%%
%Translate the users parameters into Bezier curves 
c = scaleCurviness(curviness);

%The active force length curve is made up of 5 elbow shaped sections. 
%Compute the locations of the joining point of each elbow section.

%Calculate the location of the shoulder
xDelta = 0.05*x2; %half the width of the sarcomere 0.0259, 
                       %but TM.Winter's data has a wider shoulder than
                       %this

xs    = (x2-xDelta);%x1 + 0.75*(x2-x1);

%Calculate the intermediate points located on the ascending limb
y0    = 0;   
dydx0 = 0;

y1    = 1 - dydx*(xs-x1);
dydx01= 1.25*(y1-y0)/(x1-x0);%(y1-y0)/(x1-(x0+xDelta));

x01   = x0 + 0.5*(x1-x0); %x0 + xDelta + 0.5*(x1-(x0+xDelta));
y01   = y0 + 0.5*(y1-y0);

%Calculate the intermediate points of the shallow ascending plateau
x1s   = x1 + 0.5*(xs-x1);
y1s   = y1 + 0.5*(1-y1);
dydx1s= dydx;

%dydx01c0 = 0.5*(y1s-y01)/(x1s-x01) + 0.5*(y01-y0)/(x01-x0);
%dydx01c1 = 2*( (y1-y0)/(x1-x0));
%dydx01(1-c)*dydx01c0 + c*dydx01c1; 

%x2 entered
y2 = 1;
dydx2 = 0;

%Descending limb
%x3 entered
y3 = 0;
dydx3 = 0;

x23 = (x2+xDelta) + 0.5*(x3-(x2+xDelta)); %x2 + 0.5*(x3-x2);
y23 = y2 + 0.5*(y3-y2);

%dydx23c0 = 0.5*((y23-y2)/(x23-x2)) + 0.5*((y3-y23)/(x3-x23));
%dydx23c1 = 2*(y3-y2)/(x3-x2);
dydx23   = (y3-y2)/((x3-xDelta)-(x2+xDelta)); 
%(1-c)*dydx23c0 + c*dydx23c1; 

%Compute the locations of the control points
   p0 = calcQuinticBezierCornerControlPoints( x0, ylow,   dydx0, 0, ...
                                             x01,  y01,  dydx01, 0, c);
                                         
   p1 = calcQuinticBezierCornerControlPoints(x01,  y01,  dydx01, 0, ...
                                             x1s,  y1s,  dydx1s, 0, c);
                                         
   p2 = calcQuinticBezierCornerControlPoints(x1s,  y1s,  dydx1s, 0,...
                                              x2,   y2,   dydx2, 0, c);
                                          
   p3 = calcQuinticBezierCornerControlPoints( x2,   y2,   dydx2, 0,...
                                             x23,  y23,  dydx23, 0, c);
                                         
   p4 = calcQuinticBezierCornerControlPoints(x23,  y23,  dydx23, 0,...
                                             x3,  ylow,   dydx3, 0, c);

   
   xpts = [p0(:,1) p1(:,1) p2(:,1) p3(:,1) p4(:,1)];
   ypts = [p0(:,2) p1(:,2) p2(:,2) p3(:,2) p4(:,2)];

%Create the curve structure
activeForceLengthCurve.xpts    = xpts;
activeForceLengthCurve.ypts    = ypts;

activeForceLengthCurve.xEnd         = [x0, x3];
activeForceLengthCurve.yEnd         = [ylow, ylow];
activeForceLengthCurve.dydxEnd      = [dydx0, dydx3];
activeForceLengthCurve.d2ydx2End    = [0, 0];



activeForceLengthCurve.integral = [];

if(computeIntegral == 1)
    activeForceLengthCurve.integral = ...
        createCurveIntegralStructure(activeForceLengthCurve, 1000, 1e-12,NaN);    
end

        