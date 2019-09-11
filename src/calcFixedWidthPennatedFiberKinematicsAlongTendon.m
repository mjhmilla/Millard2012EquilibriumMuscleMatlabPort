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
function fiberKinematics =...
    calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                        fiberLength,...
                                        fiberVelocity,...
                                        optimalFiberLength,...
                                        pennationAngleAtOptimalFiberLength)
%%
% This function calculates the length and velocity of the pennated 
% fiber in the direction of the tendon using the fixed width pennation
% model.
%
% This model very simply assumes that the fibers remain parallel that their  
% length and pennation angle vary such that the width and height of the
% parallelogram remain constant. This is a simple approximation to the
% constant volume property of muscle: the area of this parallelogram will
% be constant, and so if you give this a constant depth the volume of the
% parallelpiped remains constant.
%
%                       
%                          /----------/=====  -
%                         /          /        | h
%                   =====/----------/         -
%                       |<--  w -- >|
% 
%  = : tendon
%  / : fiber
%
% Note that this the terms 'height' and 'thickness' are used interchangably
% with 'width' when referring to this model. 
%
% @param fiberLengthAlongTendon  (m)
% @param fiberVelocityAlongTendon (m/s)
% @param optimalFiberLength (m)
% @param pennationAngleAtOptimalFiberLength (radians)
%
% @returns fiberKinematics, a structure with the fields:
%
%         .fiberLengthAlongTendon    (m)
%         .fiberVelocityAlongTendon  (m)
%         .pennationAngle            (radians)
%         .pennationAngularVelocity  (radians/second)
%
%%

lce       = fiberLength;    % lce*cos(aPen)
dlce      = fiberVelocity;



lceAT  = [];
dlceAT = [];
alpha  = [];
dalpha = [];

epsRoot = eps^0.5;
if(pennationAngleAtOptimalFiberLength > epsRoot)
    %%
    %Length information
    %%
    
    lopt      = optimalFiberLength;         
    alphaOpt  = pennationAngleAtOptimalFiberLength;
    
    h     = lopt*sin(alphaOpt); %the height (aka width or thickness) of the 
                                 %pennated fiber, which is constant
    lceAT = sqrt(lce*lce - h*h);
    alpha = atan2(h,lceAT);

    %%
    %Velocity information: obtained by solving:
    % [sin(alpha) ,  lce*cos(alpha)] (dlce/dt  ) = 0        [1]
    % [cos(alpha) , -lce*sin(alpha)] (dalpha/dt) = vceAT    [2]
    %
    % Eqn 1 is just dh/dt = d/dt (lce*sin(alpha)) = 0
    % Eqn 2 is just d/dt (lce*cos(alpha)) = fiberVelocityAlongTendon
    %
    % dlce/dt*sin(alpha) + lce*cos(alpha)*dalpha/dt = 0
    % dlce/dt*cos(alpha) - lce*sin(alpha)dalpha/dt = vceAT
    %

    assert(lce > epsRoot, ...
          'Impedending singularity: fiberLength -> 0!');
    assert(alpha < pi/2-epsRoot,...
          'Impedending singularity: alpha -> pi/2!');
      
    dalpha = -(dlce/lce)*tan(alpha); 
    dlceAT  = dlce*cos(alpha) - lce*sin(alpha)*dalpha; 
else
    lceAT  = lce;
    dlceAT = dlce;
    alpha  = 0;
    dalpha = 0;
end

fiberKinematics.fiberLengthAlongTendon    = lceAT;
fiberKinematics.fiberVelocityAlongTendon  = dlceAT;
fiberKinematics.pennationAngle            = alpha;
fiberKinematics.pennationAngularVelocity  = dalpha;



