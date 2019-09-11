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
%
% If you use this code in your work please cite this paper
%
%  Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). 
%    Flexing computational muscle: modeling and simulation of 
%    musculotendon dynamics. Journal of biomechanical engineering, 
%    135(2), 021005.
%
% Updates   : This function will fit the Bezier curve so that the 
%             concentric side closely matches a Hill hyperbola that
%             passes through these three points:
%
%             (-vmax        , 0           )
%             (-0.5*vmax    , fvAtHalfVmax) 
%             (1            , 1           )
%
%             This makes it far easier to create curves that are consistent
%             with slow-twitch (fvAtHalfVmax ~= 0.15) or fast-twitch muscles
%             (fvAtHalfVmax ~= 0.22) measured by Ranatunga 1984. See Fig.3
%             of Ranatunga to get the mentioned parameters
%
%   Ranatunga, K. W. (1984). The force-velocity relation of rat
%   fast-and slow-twitch muscles examined at different temperatures.
%   The Journal of Physiology, 351, 517.
%


%> This function will generate a C2 continous (continuous to the second 
%> derivative) force velocity curve of a single muscle fiber. The main 
%> function of this element is to model the amount the force enhancement or 
%> attenuation that is associated with contracting at a particular velocity.
%>
%>
%>
%>
%>
%> 
%> @param fmaxE  The normalized maximum force the fiber can generate when 
%>        is being stretched. This value is reported to range 
%>        between 1.1 and 1.8 in the literature, though all values
%>        are above 1.
%>        
%> @param dydxE    The analogous term of dydxC parameter but for the 
%>        eccentric portion of the force-velocity curve. As with
%>        the dydxC term, the physiologically accurate value for
%>        this parameter is 0, though a value of 0 is rarely used
%>        in muscle models.  If you are using an equilbrium type 
%>        model this term must be positive and greater than zero 
%>        so that the fv curve can be inverted. 
%>        <br /><br />
%>        Minimum Value: 0
%>        Maximum Value: dydxC < (fmaxE-1).
%>        <br /><br />
%>        As with the dydxC term, 
%>        the size of this term also affects the stiffness of the 
%>        integration problem for equilibrium-type muscle models: 
%>        the closer to zero this term is, the stiffer the model 
%>        will be (but only when (dlce(t)/dt)/vmax approaches 1.
%> 
%> @param dydxC  The slope of the fv(dlce(t)/dt) curve at the maximum 
%>        normalized concentric contraction velocity. Although 
%>        physiologically the value of dydxC at the maximum 
%>        concentric contracton velocity is by definition 0, a value
%>        of 0 is often used. If you are using an equilbrium type 
%>        model this term must be positive and greater than zero so
%>        that the fv curve can be inverted.
%>        <br /><br />
%>        Minimum Value: 0
%>        Maximum Value: dydxC < 1 
%>        <br /><br />
%>
%> @param flag_smoothenNonZeroDyDxC Setting this flag to 1 will linearly
%>        extrapolate the force velocity curve to match the terminal
%>        derivative of the Hill force-velocity curve for shortening 
%>        contraction velocities faster than vmax. 
%>
%> @param dydxNearE The slope of the force velocity curve as it approaches
%>         the maximum eccentric (lengthening) contraction velocity.
%>         <br /><br />
%>          Minimum Value: > dydxE
%>          Maximum Value: dydxNearE < (fmaxE-1)
%>          <br /><br />
%> 
%> @param fvAtHalfVmax: the value of the force velocity curve at half of
%>        the maximum shortening velocity. Note that this value must be 
%>        within the limits listed below.
%>         <br /><br />
%>          Minimum Value: > 0.05
%>          Maximum Value: <= 0.45
%>          <br /><br />
%>
%> @param eccCurviness     The dimensionless 'curviness' parameter that 
%>            can vary between 0 (a line) to 1 (a smooth, but 
%>            sharply bent elbow). This parameter affects only 
%>            the eccentric side of the fv curve.
%> 
%> 
%> @param curveName The name of the muscle this curve applies to. This 
%>          curve name should have the name of the muscle and the
%>          curve in it (e.g. 'bicep_fiberForceVelocityCurve') 
%>          sothat if this curve ever causes an exception, a 
%>          userfriendly error message can be displayed to the
%>          end user to help them debug their model.
%> 
%> @throws exception unless these conditions are met
%> 
%>  -0 <= dydxC < 1
%>  -dydxC < dydxNearC < 1
%>  -1 < dydxIso
%>  -dydxE < (fmaxE-1) 
%>  -dydxE < dydxNearC < (fmaxE-1)
%>  -0<= concCurviness <=0
%>  -0 <= eccCurviness <= 0
%> 
%> @return fiberForceVelocityCurve 
%>        A structure that the function calcNormalizedMuscleCurveDerivative 
%>        can use to evaluate the active force length curve value
%>        or up to the 3rd derivative.
%> 
%> 
%> Computational Costs
%> \verbatim 
%>  Without Integral :   ~8,200 flops
%> \endverbatim
%> 
%> Example:
%> @code
%>  fmaxE = 1.8;
%>  dydxC = 0.1;
%>  dydxNearC = 0.25;
%>  dydxE = 0.1;
%>  dydxNearE = 0.15;
%>  dydxIso= 5;
%>  concCurviness = 0.1;
%>  eccCurviness = 0.75;
%> 
%>  SmoothSegmentedFunction fiberFVCurve = SmoothSegmentedFunctionFactory::
%>    createFiberForceVelocityCurve(fmaxE, 
%>      dydxC, dydxNearC, dydxIso, dydxE, dydxNearE,
%>      concCurviness,  eccCurviness,false,'test');
%>  fiberFVCurve.printMuscleCurveToFile();
%> @endcode    
function fiberForceVelocityCurve =  createFiberForceVelocityCurve2018(...
                                        fmaxE,...
                                        dydxE,...
                                        dydxC,...
                                        flag_smoothenNonZeroDyDxC,...
                                        dydxNearE,...
                                        fvAtHalfVmax,...
                                        eccCurviness,...
                                        muscleName)

         


fiberForceVelocityCurve = [];
fiberForceVelocityCurve.name = sprintf('%s.%s',muscleName,...
                                       'fiberForceVelocityCurve');

vMax     = 1;                                   
vMaxC    = -vMax;   %since internally a concentric velocity is negative ...
                        %because the fiber is getting shorter.
vMaxE    =  vMax;
                  
assert( fvAtHalfVmax <= 0.45,...
        sprintf('%s: fvAtHalfVmax',...
        ' must be < 0.45',fiberForceVelocityCurve.name));
                                   
%assert( dydxIso >  1.0, ...
%  sprintf('%s: dydxIso must be < -1.1',fiberForceVelocityCurve.name));

%Ensure that the inputs are within a valid range
assert( fmaxE > 1.0, ...
  sprintf('%s: fmaxE must be greater than 1',fiberForceVelocityCurve.name));

assert( (dydxC >= 0.0 && dydxC < 1), ...
  sprintf('%s: dydxC must be greater than or equal to 0 and less than 1',...
  fiberForceVelocityCurve.name));

%assert( dydxIso > 1, ...    
%  sprintf('%s: dydxIso must be greater than (fmaxE-1)/1 (%f)',...
%  fiberForceVelocityCurve.name,((fmaxE-1.0)/1.0)));

assert( (dydxE >= 0.0 && dydxE < (fmaxE-1)), ...    
  sprintf('%s: dydxE must be greater than or equal',...
          ' to 0 and less than fmaxE-1 (%f)',...
  fiberForceVelocityCurve.name,(fmaxE-1)));

assert( (dydxNearE >= dydxE && dydxNearE < (fmaxE-1)), ...    
  sprintf('%s: dydxNearE must be greater than or',...
          ' equal to dydxE and less than fmaxE-1 (%f)',...
          fiberForceVelocityCurve.name,(fmaxE-1)));
    
%assert( (concCurviness <= 1.0 && concCurviness >= 0), ...    
%  sprintf('%s: concCurviness must be between 0 and 1', ...
%  fiberForceVelocityCurve.name));

assert( (eccCurviness <= 1.0 && eccCurviness >= 0), ...    
  sprintf('%s: eccCurviness must be between 0 and 1',...
           fiberForceVelocityCurve.name));



%We are going to set these parameters so that they fit a Hill-type hyperbolic
%fv curve that has the same slope as dydxIso.
dydxNearC     = 0;
concCurviness = 0;

%First we compute the terms that are consistent with a Hill-type concentric 
%contraction. Starting from Hill's hyperbolic equation
%
% f(w)     = (fiso*b - a*w) / (b+w)
% df(w)/dw = [-(a)*(b+w) - (b*fiso-a*w)] / (b+w)^2
%
% since this is a normalized curve
%
% at w = vMaxC the numerator goes to 0
%
% (fiso*b - a*vMaxC) / (b + vMaxC) = 0
% (fiso*b - a*vMaxC) = 0;
%  b =  a*vMaxC/fiso;
%
% Subtituting this expression for b into the expression when
%  f(wHalf) = fvAtHalfVmax yields this expression for parameter a
%
%  a = fvAtHalfVmax*w*fiso ...
%      / (vMaxC*fvAtHalfVmax + fiso*vMaxC - fiso*w);
%
%%
fiso = 1;
w = 0.5*vMaxC;
a = -fvAtHalfVmax*w*fiso ...
    / (vMaxC*fvAtHalfVmax - fiso*vMaxC + fiso*w);
b =  a*vMaxC/fiso;

yCheck  = (b*fiso-a*w)/(b+w);
assert(abs(yCheck-fvAtHalfVmax) < sqrt(eps));

w = 0*vMaxC;
dydxIso =(-(a)*(b+w) - (b*fiso-a*w)) / ((b+w)*(b+w));

w         = 0.9*vMaxC;
dydxNearC = (-(a)*(b+w) - (b*fiso-a*w)) / ((b+w)*(b+w));

assert( (abs(dydxNearC) > abs(dydxC) || abs(dydxNearC) < abs(1/vMaxC)), ...
  sprintf('%s: dydxNearC must be greater than or equal to 0 and less than 1',...
  fiberForceVelocityCurve.name));


%Solve for the concCurviness that results in the Bezier curve that minimizes
%the error between the Bezier curve and Hill's concentric contraction equations
%at 10 points between omega = 0, and omega = 0.9*vMax 
iter    = 1;
iterMax = 100;
tol     = sqrt(eps);

xNearC =  0.9*vMaxC;
yNearC = (b*fiso-a*w)/(b+w);
xIso = 0;
yIso = 1.0;

cC = 0.5;
pts = calcQuinticBezierCornerControlPoints(xIso,  yIso,   dydxIso, 0,...
                                         xNearC,yNearC, dydxNearC, 0,...
                                                                  cC);
curve.xpts = pts(:,1);
curve.ypts = pts(:,2);
curve.xEnd = [0,vMaxC];
curve.yEnd = [1,0];
curve.dydxEnd = [dydxIso, dydxNearC];

%Get the initial error between Hill and the Bezier curve for cC = 0.5.
nSample = 10;
f = 0;
for j=1:1:nSample
    w    = (j-1)*vMaxC/nSample;
    yHill = (b*fiso-a*w)/(b+w);        
    f = f + abs( calcBezierYFcnXDerivative(w, curve, 0) - yHill);
end

fBest = f;
cCBest = cC;

h = 0.25;
curveLeft.xpts = [];
curveLeft.ypts = [];
curveLeft.xEnd = [0,1];
curveLeft.yEnd = [1,0];
curveLeft.dydxEnd = [dydxIso, dydxNearC];

curveRight.xpts = [];
curveRight.ypts = [];
curveRight.xEnd = [0,1];
curveRight.yEnd = [1,0];
curveRight.dydxEnd = [dydxIso,dydxNearC];

%Use the bisection method to find the best curviness value
for i=1:1:10


  cCLeft = cC-h;
  ptsLeft =...
    calcQuinticBezierCornerControlPoints( xIso,  yIso,  dydxIso, 0,...
                                        xNearC,yNearC,dydxNearC, 0,...
                                        cCLeft);
  curveLeft.xpts = ptsLeft(:,1);
  curveLeft.ypts = ptsLeft(:,2);
  
  cCRight = cC+h;
  ptsRight =...
    calcQuinticBezierCornerControlPoints(xIso,  yIso,  dydxIso, 0,...
                                       xNearC,yNearC,dydxNearC, 0,...
                                       cCRight);
  curveRight.xpts = ptsRight(:,1);
  curveRight.ypts = ptsRight(:,2);

  %Compute the error at 10 points between -vMax and 0. 
  fLeft = 0;
  for j=1:1:nSample
    w    =  (j-1)*vMaxC/(nSample-1);
    yHill = (b*fiso-a*w)/(b+w);      
    fLeft = fLeft + abs( calcBezierYFcnXDerivative(w, curveLeft, 0) - yHill);
  end

  fRight = 0;
  for j=1:1:nSample
    w    = (j-1)*vMaxC/(nSample-1);
    yHill = (b*fiso-a*w)/(b+w);     
    fRight = fRight + abs( calcBezierYFcnXDerivative(w, curveRight, 0) - yHill);
  end
         
  %disp(sprintf('f: %e, fl: %e, fr: %e, cC: %e',f,fLeft,fRight,cC));  
  %Update the current solution
  if(abs(fLeft) < abs(f))
      f = fLeft;
      cC = cCLeft;
  end
  
  if(abs(fRight) < abs(f))
     f = fRight;
     cC = cCRight; 
  end
  
  h = h/2;
  iter = iter+1;
end

%assert( (abs(f) <= tol), ...
%  sprintf('%s: Fitting the curve to the concentric force-velocity curve of ',...
%    'A.V. Hill failed.',...
%  fiberForceVelocityCurve.name));


%Translate the users parameters into Bezier point locations
%cC = scaleCurviness(concCurviness);
cE = scaleCurviness(eccCurviness);

%Compute the concentric control point locations
xC   = vMaxC;
yC   = 0;
if(flag_smoothenNonZeroDyDxC==1)
  dydxC = 0.5*dydxNearC;
end

xNearC = 0.9*vMaxC;
yNearC = yC + 0.5*dydxNearC*(xNearC-xC) + 0.5*dydxC*(xNearC-xC);



xIso = 0;
yIso = 1;

concPts1 =...
    calcQuinticBezierCornerControlPoints(    xC,     yC,    dydxC,0,...
                                         xNearC, yNearC,dydxNearC,0, cC);
concPts2 =...
    calcQuinticBezierCornerControlPoints(xNearC,  yNearC,dydxNearC, 0, ...
                                           xIso,    yIso,  dydxIso, 0,  cC);


yIsoH    = fmaxE - (0.1*vMaxE*dydxE) - (dydxNearE*0.9*vMaxE);
dydxIsoH = max(5.0*dydxNearE, 5.0*dydxIso);
xIsoH    = (yIsoH-1.0)/dydxIsoH;

xIsoHm = 2*0.5*xIsoH;
yIsoHm  = 1 + (yIsoH-1)*0.5;

eccPts1  = ...
    calcQuinticBezierCornerControlPoints(  xIso,   yIso,  dydxIso, 0,...
                                         xIsoHm, yIsoHm, dydxIsoH, 0, cE);

xE   = vMaxE;
yE   = fmaxE;

xNearE    = 0.9*vMaxE;
dydxNearE = dydxNearE/vMaxE;
yNearE    = yE + 0.5*dydxNearE*(xNearE-xE) + 0.5*dydxE*(xNearE-xE);


eccPts2  = ...
    calcQuinticBezierCornerControlPoints(xIsoHm, yIsoHm,  dydxIsoH, 0,...
                                         xNearE, yNearE, dydxNearE, 0,...
                                          cE);
                                       



eccPts3  = ...
    calcQuinticBezierCornerControlPoints(xNearE, yNearE, dydxNearE, 0, ...
                                             xE,     yE,     dydxE, 0, cE);

                                                                              

xpts = [concPts1(:,1) concPts2(:,1) eccPts1(:,1) eccPts2(:,1) eccPts3(:,1)];
ypts = [concPts1(:,2) concPts2(:,2) eccPts1(:,2) eccPts2(:,2) eccPts3(:,2)];

fiberForceVelocityCurve.xpts = xpts;
fiberForceVelocityCurve.ypts = ypts;

fiberForceVelocityCurve.xEnd = [xC xE];
fiberForceVelocityCurve.yEnd = [yC yE];
fiberForceVelocityCurve.dydxEnd  = [dydxC dydxE];
fiberForceVelocityCurve.d2ydx2End= [0,0];
fiberForceVelocityCurve.integral = [];
