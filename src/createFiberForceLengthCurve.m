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
function fiberForceLengthCurve = createFiberForceLengthCurve(...
                                       eZero, eIso, ...
                                       kLow, kIso, curviness,...
                                       computeIntegral, muscleName)
%%

%This function will generate a C2 continuous curve that fits a fiber's 
%tensile force length curve.
%
%@param eZero The fiber strain at which the fiber begins to develop force.
%			 Thus an e0 of 0.0 means that the fiber will start to develop
%			 passive force when it has a normalized length of 1.0. Note
%			 that e0 can be postive or negative.
%
%@param eIso The fiber strain at which the fiber develops 1 unit of 
%			normalized force (1 maximum isometric force). Note that the 
%			'1' is left off. Thus an e0 of 0.6 means that the fiber 
%			will develop an 1 normalized force unit when it is strained 
%			by 60% of its resting length, or to a normalized length of 
%			1.6
%
%@param kLow   The normalized stiffness (or slope) of the fiber curve 
%			  close to the location where the force-length curve 
%			  approaches a normalized force of 0. This is usually 
%			  chosen to be a small, but non-zero fraction of kIso 
%			  (kLow = 0.025 kIso is typical).
%
%@param kIso   The normalized stiffness (or slope) of the fiber curve 
%			  when the fiber is strained by eIso (or has a length of 
%			  1+eIso) under a load of 1 maximum isometric unit of force.
%
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
%				  curve in it (e.g. "bicep_fiberForceLengthCurve") 
%				  so that if this curve ever causes an exception, a 
%				  userfriendly error message can be displayed to the
%				  end user to help them debug their model.
%
%@throws exception unless the following conditions are met
%	-eIso > eZero            
%	-kIso > 1/(eIso-eZero)
%	-0 < kLow < kIso
%	-0 <= curviness <= 1
%
%@return fiberForceLengthCurve 
%       A structure that the function calcNormalizedMuscleCurveDerivative 
%       can use to evaluate the active force length curve value
%       or up to the 3rd derivative.
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
%	 eIso      = 0.6;
%	 eZero     = 0.0;
%	 kIso      = 4.0/(eIso-eZero);
%	 kNearZero = 0.025*kIso
%	 c         = 0.5;
%
%	 fiberFLCurve = createFiberForceLengthCurve(eZero, eIso,
%								  kLow, kIso, c, true,"test");
%@endcode
%%
fiberForceLengthCurve = [];
fiberForceLengthCurve.name = sprintf('%s.%s',muscleName,'fiberForceLengthCurve');

%%
%Check the input arguments
%%
assert( eIso > eZero , ...
    sprintf('%s: The following must hold: eIso  > eZero',fiberForceLengthCurve.name));

assert( kIso > (1.0/(eIso-eZero)) , ...
    sprintf('%s: kiso must be greater than 1/(eIso-eZero) (%f)',...
    fiberForceLengthCurve.name, (1.0/(eIso-eZero))));

assert(kLow > 0.0 && kLow < 1/(eIso-eZero),...
    sprintf('%s: kLow must be greater than 0 and less than or equal to 1',...
    fiberForceLengthCurve.name));

assert( (curviness>=0 && curviness <= 1),...      
    sprintf('%s: curviness must be between 0.0 and 1.0',...
            fiberForceLengthCurve.name));


%%
%Translate the user parameters to quintic Bezier points
%%
c = scaleCurviness(curviness);
xZero = 1+eZero;
yZero = 0;

xIso = 1 + eIso;
yIso = 1;

deltaX = min(0.1*(1.0/kIso), 0.1*(xIso-xZero));

xLow     = xZero + deltaX;
xfoot    = xZero + 0.5*(xLow-xZero);
yfoot    = 0;
yLow     = yfoot + kLow*(xLow-xfoot);

%Compute the Quintic Bezier control points
p0 = calcQuinticBezierCornerControlPoints(xZero, yZero,   0, 0, ...
                                           xLow, yLow, kLow, 0,c);

p1 =  calcQuinticBezierCornerControlPoints(xLow, yLow, kLow, 0, ...
                                           xIso, yIso, kIso, 0, c);
xpts = [p0(:,1) p1(:,1)];
ypts = [p0(:,2) p1(:,2)];



%Create the curve structure
fiberForceLengthCurve.xpts    = xpts;
fiberForceLengthCurve.ypts    = ypts;

fiberForceLengthCurve.xEnd         = [xZero, xIso];
fiberForceLengthCurve.yEnd         = [yZero, yIso];
fiberForceLengthCurve.dydxEnd      = [0, kIso];
fiberForceLengthCurve.d2ydx2End    = [0, 0];

fiberForceLengthCurve.integral = [];

if(computeIntegral == 1)
    xScaling = eIso;
    fiberForceLengthCurve.integral = ...
        createCurveIntegralStructure(fiberForceLengthCurve, ...
                                     1000,...
                                     1e-12,...
                                     xScaling);    
end
                                   