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
function tendonForceLengthCurve = ...
          createTendonForceLengthCurve( eIso, kIso, ...
                                        fToe, curviness, ...
                                        computeIntegral, ...
                                        muscleName)
%%

%Will generate a C2 continous (continuous to the second derivative) 
%curve in a MuscleFunctionObject object that fits a tendon's tensile 
%force length curve. 
%
%
%
%@param eIso   The tendon strain at which the tendon develops 1 unit
%			of normalized force (1 maximum isometric force). Note that 
%			the'1' is left off. Thus an e0 of 0.04 means that the tendon 
%			will develop an 1 normalized force unit when it is strained 
%			by 4% of its resting length, at a normalized length of 
%			1.04
%
%@param kIso    The normalized stiffness (or slope) of the tendon
%				curve when the tendon is strained by e0 
%				(or has a length of 1+e0) under a load of 1 maximum
%				isometric unit of force.        
%
%@param fToe    The normalized force at which the tendon smoothly
%			   transitions from the curved low stiffness region to 
%			   the linear stiffness region.
%
%@param curviness    The dimensionless 'curviness' parameter that 
%					can vary between 0 (a line) to 1 (a smooth, but 
%					sharply bent elbow)
%
%@param computeIntegral  If this is true, the integral for this curve
%						is numerically calculated and splined. If false, 
%						this integral is not computed, and a call to 
%						.calcIntegral will throw an exception
%
% @param curveName The name of the muscle this curve applies to. This 
%				  curve name should have the name of the muscle and the
%				  curve in it (e.g. 'bicep_tendonForceLengthCurve') 
%				  sothat if this curve ever causes an exception, a 
%				  userfriendly error message can be displayed to the
%				  end user to help them debug their model.
%
%@throws SimTK::Exception unless the following conditions are met:
%	-0 < fToe < 1
%	-e0 > 0
%	-kiso > 1/e0
%	-0 <= curviness <= 1
%
%@return SmoothSegmentedFunction*
%
%\image html fig_SmoothSegmentedFunctionFactory_fseCurve.png
%
%
%<B>Computational Costs</B>
%\verbatim 
%	Without Integral :   ~4,100 flops
%	With Integral    : ~174,100 flops
%\endverbatim
%
%<B>Example:</B>
%@code
%	e0   = 0.04;
%	kiso = 42.79679348815859;
%	fToe = 1.0/3.0
%	c    = 0.75;
%
%	SmoothSegmentedFunction* tendonCurve = SmoothSegmentedFunctionFactory::
%										createTendonForceLengthCurve(
%										  e0,kiso,fToe,c,true,'test');
%	tendonCurve.printMuscleCurveToFile();  
%@endcode
%%
tendonForceLengthCurve = [];
tendonForceLengthCurve.name = sprintf('%s.%s',muscleName,'tendonForceLengthCurve');

%Check the input arguments
%eIso>0 
assert( eIso>0 , ... 
  sprintf('%s: eIso must be greater than 0, but %f was entered',... 
  tendonForceLengthCurve.name,eIso));

assert( (fToe>0 && fToe < 1) , ...
  sprintf('%s: fToe must be greater than 0 and less than 1, but %f was entered',... 
  tendonForceLengthCurve.name,fToe));

assert( kIso > (1/eIso) , ... 
  sprintf('%s : kIso must be greater than 1/eIso, (%f), but kIso (%f) was entered',...
  tendonForceLengthCurve.name, (1/eIso),kIso));

assert( (curviness>=0 && curviness <= 1) , ... 
  sprintf('%s : curviness must be between 0.0 and 1.0, but %f was entered',...
  tendonForceLengthCurve.name,curviness));



%Translate the user parameters to quintic Bezier points
c = scaleCurviness(curviness);
x0 = 1.0;
y0 = 0;
dydx0 = 0;

xIso = 1.0 + eIso;
yIso = 1;
dydxIso = kIso;

%Location where the curved section becomes linear
yToe = fToe;
xToe = (yToe-1)/kIso + xIso;


%To limit the 2nd derivative of the toe region the line it tends to
%has to intersect the x axis to the right of the origin
xFoot = 1.0+(xToe-1.0)/10.0;
yFoot = 0;
dydxToe = (yToe-yFoot)/(xToe-xFoot);

%Compute the location of the corner formed by the average slope of the
%toe and the slope of the linear section
yToeMid = yToe*0.5;
xToeMid = (yToeMid-yIso)/kIso + xIso;
dydxToeMid = (yToeMid-yFoot)/(xToeMid-xFoot);

%Compute the location of the control point to the left of the corner
xToeCtrl = xFoot + 0.5*(xToeMid-xFoot); 
yToeCtrl = yFoot + dydxToeMid*(xToeCtrl-xFoot);



%Compute the Quintic Bezier control points
p0 = calcQuinticBezierCornerControlPoints(x0,      y0,     dydx0, 0,...
                                    xToeCtrl,yToeCtrl,dydxToeMid, 0, c);
p1 = calcQuinticBezierCornerControlPoints(xToeCtrl, yToeCtrl, dydxToeMid, 0,...
                                              xToe,     yToe,    dydxIso, 0, c);
xpts = [p0(:,1),p1(:,1)];
ypts = [p0(:,2),p1(:,2)];

tendonForceLengthCurve.xpts    = xpts;
tendonForceLengthCurve.ypts    = ypts;

tendonForceLengthCurve.xEnd         = [x0, xToe];
tendonForceLengthCurve.yEnd         = [y0, yToe];
tendonForceLengthCurve.dydxEnd      = [dydx0, dydxIso];
tendonForceLengthCurve.d2ydx2End    = [0, 0];

tendonForceLengthCurve.integral = [];

if(computeIntegral == 1)
    xScaling = eIso;
    tendonForceLengthCurve.integral = ...
        createCurveIntegralStructure(tendonForceLengthCurve, ...
                                     1000,...
                                     1e-12,...
                                     xScaling);    
end


