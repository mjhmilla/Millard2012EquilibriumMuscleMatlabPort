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
function val = calcBezierYFcnXDerivative(x, curveParams, der)
%%
% This function takes the control points for the Bezier curves
% x(t) and f(t) and treats it as a function:
%
%        val = y(x)
%
% To do this we assume that
% 1. x(t) is a monotonic increasing function
% 2. d/dt x(t) > 0
%
% In addition, to maintain continuity we assume that the second derivative
% of the curve goes to 0 at the end points so that we can linearly
% extrapolate the curve outside of [xmin, xmax]
%
% @param x: value to 
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
% @param der: integer defining the order of the desired derivative. 
%             The value of der is restricted to:
%             
%             0: y = f(x)
%             1: dy/dx  
%             2: d^2y/dx^2 
%             3: d^3y/dx^3
%
% @returns y : d^n/dx^n y(x), the n^th derivative of the function y(x)
%%

val = NaN;

assert((der >= -1 && der <= 3),'der must be within [0,3]');

xpts = curveParams.xpts;
ypts = curveParams.ypts;

nrow = size(xpts,1);
ncol = size(xpts,2);
xmin = min(min(xpts));
xmax = max(max(xpts));

if( (x < xmin) || (x > xmax) )
    idxEnd = NaN;
    if(x <= xmin)
       idxEnd = 1; 
    else
       idxEnd = 2;
    end
    
    %evaluate the desired derivative of y at x
    switch der
        case -1
            assert(isempty(curveParams.integral)==0,...
               'Integral function has not been computed for this curve');
            if(x <= xmin)
               val = 0; 
            else
               len = length(curveParams.integral.xptsN);
               y0  = curveParams.integral.yptsN(len);
               f1  = curveParams.integral.y1ptsN(len);
               f2  = curveParams.integral.y2ptsN(len); 
               x0  = curveParams.integral.xptsN(len);
               
               val = y0 + f1*(x-x0) + (0.5*f2)*(x-x0)^2;
            end            
        case 0
            x0   = curveParams.xEnd(idxEnd);
            y0   = curveParams.yEnd(idxEnd);
            dydx = curveParams.dydxEnd(idxEnd);
                        
            val = dydx*(x-x0) + y0;
             
        case 1
            dydx = curveParams.dydxEnd(idxEnd);
            val  = dydx;
        otherwise
            val = 0;
    end        
else
    %Find the spline section that x is in
    options.moduloRange = [];
    options.tol    = eps; 
    xInterval      = [xpts(   1,:); ...
                      xpts(nrow,:)];
    col =  calcIndex(x, xInterval, options);
                               
    %Extract the vector of control points, and calculate
    %the control points that define the derivative curve
    xV  = xpts(:,col);
    x1V = diff(xV) .*(nrow-1);
    
    
    %Find the value of u that corresponds to the desired value of x
    u      = (x-xpts(1, col)) / (xpts(nrow, col)-xpts(1, col));
    iter    = 1;
    iterMax = 100;    
    tol     = eps*10;
    err     = tol*10;
    
    while iter < iterMax && abs(err) > tol
       err   = calc1DBezierCurveValue(u, xV) - x;
       derr = calc1DBezierCurveValue(u, x1V);
       
       if(abs(err) > tol && abs(derr) > eps)          
           du = -err/derr;
           u  = u + du;
           
           %For very nonlinear curves Newton's method can
           %become pathological. If u is outside [0,1] we
           %kick it back into the interval by some random small
           %amount.
           if(u < 0 || u > 1.0)
                if(u < 0)
                   u = clampU(u);
                   u = u + rand(1)*0.1;
                else
                   u = clampU(u);
                   u = u - rand(1)*0.1;                    
                end
           end
       end
     
       iter = iter+1;
    end
        
    %Evaluate the desired derivative of y   
    switch der
        case -1
             assert(isempty(curveParams.integral)==0,...
                  'Integral function has not been computed for this curve');
           
             [y0 y1] = calc5thOrderInterp(x, ...
                                    curveParams.integral.xptsN,...
                                    curveParams.integral.yptsN,...
                                    curveParams.integral.y1ptsN,...
                                    curveParams.integral.y2ptsN);
                 val = y0;   
              
        case 0            
            yV  =  ypts(:,col);
            val = calc1DBezierCurveValue(u, yV);
        case 1
            yV  =  ypts(:,col);
            y1V =  diff(yV) .*(nrow-1);
            x1  = calc1DBezierCurveValue(u, x1V);
            y1  = calc1DBezierCurveValue(u, y1V);
            
            val = y1/x1;
        case 2
            x2V =  diff(x1V).*(nrow-2); 
            
            yV  =  ypts(:,col);
            y1V =  diff(yV) .*(nrow-1);
            y2V =  diff(y1V).*(nrow-2);
    
            x1  = calc1DBezierCurveValue(u, x1V);
            y1  = calc1DBezierCurveValue(u, y1V);
            x2  = calc1DBezierCurveValue(u, x2V);
            y2  = calc1DBezierCurveValue(u, y2V);
            
            t1 = 1/x1;
            t3 = x1*x1;
            
            %val = ((x1*y2-x2*y1)/(x1*x1))/x1;
            val = (y2 * t1 - y1 / t3 * x2) * t1;
        case 3
            x2V = diff(x1V).*(nrow-2);
            x3V = diff(x2V).*(nrow-3);
            
            yV  =  ypts(:,col);
            y1V =  diff(yV) .*(nrow-1);
            y2V =  diff(y1V).*(nrow-2);
            y3V =  diff(y2V).*(nrow-3);
            
            x1  = calc1DBezierCurveValue(u, x1V);
            y1  = calc1DBezierCurveValue(u, y1V);
            
            x2  = calc1DBezierCurveValue(u, x2V);
            y2  = calc1DBezierCurveValue(u, y2V);
            
            x3  = calc1DBezierCurveValue(u, x3V);
            y3  = calc1DBezierCurveValue(u, y3V);
            
            t1 = 1 / x1;
            t3 = x1*x1;
            t4 = 1 / t3;
            t11 = x2*x2;
            t14 = y1 * t4;
            
            %val = (  (x1*x1*y3-2*x1*x2*y2+(2*x2*x2-x1*x2)*y1)...
            %        /(x1*x1*x1);
            %val =(    (x2*y2+x1*y3-x3*y1-x2*y2)/(x1*x1) ...
            %              - 2*(x1*y2-x2*y1)*x2/(x1*x1*x1)  )/x1;
            val = ((y3*t1 - 2*y2*t4*x2 ...
                  + 2*y1/t3/x1 * t11 - t14 * x3) * t1 ...
                - (y2*t1 - t14*x2)*t4*x2) * t1;              
    end
            
    
end
    