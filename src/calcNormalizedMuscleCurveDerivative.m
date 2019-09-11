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
function val = calcNormalizedMuscleCurveDerivative(x,der,curveIndex,muscleCurveParams)
%%
% This function evaluates a musculotendon curve at x. The order of the
% derivative to evaluate is set using der, which can range from 0 to 3.
% All curves are represented using quintic Bezier curves except for the
% tendon strain energy curve, which is just an interpolation of a numerical
% integral within the defined range, and a quadratic extrapolation outside
% of the curve interval.
% 
% @param x: the argument of the curve function
%
% @param der: the desired derivative of the curve. 
%             values of der must be [0,3], where 0 corresponds to the value
%             of the curve
%
% @param curveIndex: an integer that accesses the correct curve. These
%                    indices are stored in muscleCurveParams(1).index
%                    Currently these are the indices:
%
%  index    curve
%  1        fiber active force length curve
%  2        fiber passive force length curve
%  3        fiber force velocity curve
%  4        fiber force velocity inverse curve
%  5        tendon force length curve
%  6        tendon force length strain energy curve
%
% @param muscleCurveParams: the structure that contains, for each curve,
%                           the control points required to evaluate the 
%                           splines that define each curve and its first 3 
%                           derivatives.
%
% @returns val: the value of the curve derivative at x
%%
val = NaN;

if(curveIndex ~= muscleCurveParams(1).index.fte)
    i=curveIndex;
    xpts = muscleCurveParams(i).xpts;
    ypts = muscleCurveParams(i).ypts;
    val  = calcBezierYFcnXDerivative(x, xpts, ypts, der);
else
    %we're evaluating the tendon strain energy curve.
    idxE=curveIndex;
    xmin = min(muscleCurveParams(idxE).xptsN);
    xmax = max(muscleCurveParams(idxE).xptsN);

    switch der
        case 0
            if(x < xmin)
                %tendon is slack, no stored energy
                val = 0;
            elseif(x > xmax)
                %energy stored is a quadratic function of the
                %terminal slope of the tendon curve
                x0 = xmax;
                y0 = max(muscleCurveParams(idxE).y1pts);                
                
                y1pts = muscleCurveParams(idxE).y1pts;
                nrowT = size(y1pts,1);
                ncolT = size(y1pts,2);
                
                y1 = y1pts(nrowT,ncolT);
                y2pts = diff(y1pts).*(nrowT-1);
                y2    = y2pts(nrowT-1,ncolT);
                
                val = y0 + y1*(x-xmax) + (0.5*y2)*(x-xmax)^2;
                
                
            else
                [y0 y1] = calc5thOrderInterp(x, ...
                                    muscleCurveParams(idxE).xptsN,...
                                    muscleCurveParams(idxE).yptsN,...
                                    muscleCurveParams(idxE).y1ptsN,...
                                    muscleCurveParams(idxE).y2ptsN);
                 val = y0;               
            end
        case 1
            xpts = muscleCurveParams(idxE).xpts;
            ypts = muscleCurveParams(idxE).y1pts;
            val  = calcBezierYFcnXDerivative(x, xpts, ypts, 0);
            
        case 2
            xpts = muscleCurveParams(idxE).xpts;
            ypts = muscleCurveParams(idxE).y1pts;
            val  = calcBezierYFcnXDerivative(x, xpts, ypts, 1);
            
        case 3
            xpts = muscleCurveParams(idxE).xpts;
            ypts = muscleCurveParams(idxE).y1pts;
            val  = calcBezierYFcnXDerivative(x, xpts, ypts, 2);
    end
end