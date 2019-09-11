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
function fiberForceVelocityCurve =  createFiberForceVelocityCurve2012(...
                                        fmaxE, ...
                                        dydxC,dydxNearC, ...
                                        dydxIso, ...
                                        dydxE, dydxNearE,...
                                        concCurviness, eccCurviness,...
                                        computeIntegral, muscleName)
%%
%This function will generate a C2 continous (continuous to the second 
%derivative) force velocity curve of a single muscle fiber. The main 
%function of this element is to model the amount the force enhancement or 
%attenuation that is associated with contracting at a particular velocity.
%
%@param fmaxE  The normalized maximum force the fiber can generate when 
%			  is being stretched. This value is reported to range 
%			  between 1.1 and 1.8 in the literature, though all values
%			  are above 1.
%
%@param dydxC  The slope of the fv(dlce(t)/dt) curve at the maximum 
%			  normalized concentric contraction velocity. Although 
%			  physiologically the value of dydxC at the maximum 
%			  concentric contracton velocity is by definition 0, a value
%			  of 0 is often used. If you are using an equilbrium type 
%			  model this term must be positive and greater than zero so
%			  that the fv curve can be inverted.
%			  <br /><br />
%			  Minimum Value: 0
%			  Maximum Value: dydxC < 1 
%			  <br /><br />
%
%@param dydxNearC The slope of the force velocity curve as it approaches
%				 the maximum concentric (shortening) contraction velocity.
%				 <br /><br />
%				  Minimum Value: > dydxC
%				  Maximum Value: dydxNearC < 1 
%				  <br /><br />
%
%
%@param dydxIso  The slope of the fv curve when dlce(t)/dt = 0. 
%				<br /><br />
%				Minimim Value: dydxIso > 1.0
%				Maximum Value: dydxIso < Inf
%				
%@param dydxE    The analogous term of dydxC parameter but for the 
%				eccentric portion of the force-velocity curve. As with
%				the dydxC term, the physiologically accurate value for
%				this parameter is 0, though a value of 0 is rarely used
%				in muscle models.  If you are using an equilbrium type 
%				model this term must be positive and greater than zero 
%				so that the fv curve can be inverted. 
%				<br /><br />
%				Minimum Value: 0
%				Maximum Value: dydxC < (fmaxE-1).
%				<br /><br />
%				As with the dydxC term, 
%				the size of this term also affects the stiffness of the 
%				integration problem for equilibrium-type muscle models: 
%				the closer to zero this term is, the stiffer the model 
%				will be (but only when (dlce(t)/dt)/vmax approaches 1.
%
%@param dydxNearE The slope of the force velocity curve as it approaches
%				 the maximum eccentric (lengthening) contraction velocity.
%				 <br /><br />
%				  Minimum Value: > dydxE
%				  Maximum Value: dydxNearE < (fmaxE-1)
%				  <br /><br />
%
%
%@param concCurviness    The dimensionless 'curviness' parameter that 
%						can vary between 0 (a line) to 1 (a smooth, but 
%						sharply bent elbow). This parameter affects only
%						the concentric side of the fv curve.
%
%@param eccCurviness     The dimensionless 'curviness' parameter that 
%						can vary between 0 (a line) to 1 (a smooth, but 
%						sharply bent elbow). This parameter affects only 
%						the eccentric side of the fv curve.
%
%@param computeIntegral  If this is true, the integral for this curve
%						is numerically calculated and splined. If false, 
%						this integral is not computed, and a call to 
%						SmoothSegmentedFunction::calcIntegral() will throw 
%						an exception
%
%@param curveName The name of the muscle this curve applies to. This 
%				  curve name should have the name of the muscle and the
%				  curve in it (e.g. 'bicep_fiberForceVelocityCurve') 
%				  sothat if this curve ever causes an exception, a 
%				  userfriendly error message can be displayed to the
%				  end user to help them debug their model.
%
%@throws exception unless these conditions are met
%
%	-0 <= dydxC < 1
%	-dydxC < dydxNearC < 1
%	-1 < dydxIso
%	-dydxE < (fmaxE-1) 
%	-dydxE < dydxNearC < (fmaxE-1)
%	-0<= concCurviness <=0
%	-0 <= eccCurviness <= 0
%
%@return fiberForceVelocityCurve 
%       A structure that the function calcNormalizedMuscleCurveDerivative 
%       can use to evaluate the active force length curve value
%       or up to the 3rd derivative.
%
%
%<B>Computational Costs</B>
%\verbatim 
%	Without Integral :   ~8,200 flops
%	With Integral    : ~348,200 flops
%\endverbatim
%
%<B>Example:</B>
%@code
%	fmaxE = 1.8;
%	dydxC = 0.1;
%	dydxNearC = 0.25;
%	dydxE = 0.1;
%	dydxNearE = 0.15;
%	dydxIso= 5;
%	concCurviness = 0.1;
%	eccCurviness = 0.75;
%
%	SmoothSegmentedFunction fiberFVCurve = SmoothSegmentedFunctionFactory::
%		createFiberForceVelocityCurve(fmaxE, 
%			dydxC, dydxNearC, dydxIso, dydxE, dydxNearE,
%			concCurviness,  eccCurviness,false,'test');
%	fiberFVCurve.printMuscleCurveToFile();
%@endcode             
%%
fiberForceVelocityCurve = [];
fiberForceVelocityCurve.name = sprintf('%s.%s',muscleName,...
                                       'fiberForceVelocityCurve');

%Ensure that the inputs are within a valid range
assert( fmaxE > 1.0, ...
  sprintf('%s: fmaxE must be greater than 1',fiberForceVelocityCurve.name));

assert( (dydxC >= 0.0 && dydxC < 1), ...
  sprintf('%s: dydxC must be greater than or equal to 0 and less than 1',...
  fiberForceVelocityCurve.name));

assert( (dydxNearC >= dydxC && dydxNearC <= 1), ...
  sprintf('%s: dydxNearC must be greater than or equal to 0 and less than 1',...
  fiberForceVelocityCurve.name));

assert( dydxIso > 1, ...    
  sprintf('%s: dydxIso must be greater than (fmaxE-1)/1 (%f)',...
  fiberForceVelocityCurve.name,((fmaxE-1.0)/1.0)));

assert( (dydxE >= 0.0 && dydxE < (fmaxE-1)), ...    
  sprintf('%s: dydxE must be greater than or equal to 0 and less than fmaxE-1 (%f)',...
  fiberForceVelocityCurve.name,(fmaxE-1)));

assert( (dydxNearE >= dydxE && dydxNearE < (fmaxE-1)), ...    
  sprintf('%s: dydxNearE must be greater than or equal to dydxE and less than fmaxE-1 (%f)',...
  fiberForceVelocityCurve.name,(fmaxE-1)));
    
assert( (concCurviness <= 1.0 && concCurviness >= 0), ...    
  sprintf('%s: concCurviness must be between 0 and 1', ...
  fiberForceVelocityCurve.name));

assert( (eccCurviness <= 1.0 && eccCurviness >= 0), ...    
  sprintf('%s: eccCurviness must be between 0 and 1',...
           fiberForceVelocityCurve.name));



%Translate the users parameters into Bezier point locations
cC = scaleCurviness(concCurviness);
cE = scaleCurviness(eccCurviness);

%Compute the concentric control point locations
xC   = -1;
yC   = 0;

xNearC = -0.9;
yNearC = yC + 0.5*dydxNearC*(xNearC-xC) + 0.5*dydxC*(xNearC-xC);

xIso = 0;
yIso = 1;

xE   = 1;
yE   = fmaxE;

xNearE = 0.9;
yNearE = yE + 0.5*dydxNearE*(xNearE-xE) + 0.5*dydxE*(xNearE-xE);


concPts1 =...
    calcQuinticBezierCornerControlPoints(    xC,     yC,    dydxC,0,... 
                                         xNearC, yNearC,dydxNearC,0,cC);
                                     
concPts2 =...
    calcQuinticBezierCornerControlPoints(xNearC,yNearC,dydxNearC, 0,...
                                           xIso,  yIso,  dydxIso, 0, cC);
eccPts1  = ...
    calcQuinticBezierCornerControlPoints(xIso,     yIso,   dydxIso, 0,... 
                                         xNearE, yNearE, dydxNearE, 0, cE);

eccPts2  = ...
    calcQuinticBezierCornerControlPoints(xNearE, yNearE, dydxNearE, 0,... 
                                             xE,     yE,     dydxE, 0, cE);


xpts = [concPts1(:,1) concPts2(:,1) eccPts1(:,1) eccPts2(:,1)];
ypts = [concPts1(:,2) concPts2(:,2) eccPts1(:,2) eccPts2(:,2)];

fiberForceVelocityCurve.xpts = xpts;
fiberForceVelocityCurve.ypts = ypts;

fiberForceVelocityCurve.xEnd = [xC xE];
fiberForceVelocityCurve.yEnd = [yC yE];
fiberForceVelocityCurve.dydxEnd= [dydxC dydxE];
fiberForceVelocityCurve.d2ydx2End=[0,0];

fiberForceVelocityCurve.integral = [];

if(computeIntegral == 1)
    fiberForceVelocityCurve.integral = ...
        createCurveIntegralStructure(fiberForceVelocityCurve, 1000, 1e-12,NaN);    
end
