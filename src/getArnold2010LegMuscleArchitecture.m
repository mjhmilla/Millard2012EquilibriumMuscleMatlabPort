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
function arnold2010LegMuscleArchitecture = getArnold2010LegMuscleArchitecture(unitsMKSN)
%%
% This function reads in Table 1 of
%
% Arnold, E. M., Ward, S. R., Lieber, R. L., & Delp, S. L. (2010). 
% A model of the lower limb for analysis of human movement. 
% Annals of biomedical engineering, 38(2), 269-279.
%
% into a structure. 
%
% @param unitsType: 0 : N, cm, deg as in the table
%                   1 : N, m, rad
%
% @returns arnold2010LegMuscleArchitecture a structure with the fields
%          that correspond directly to the columns of Table 1. Note too  
%          each row corresponds directly to the rows of Table 1, and thus
%          some rows refer to groups of muscles (like row 3) and do not
%          contain all of the fields.
%
%      field                 units        derived how?
%      .abbrevation          n/a          n/a
%      .names                n/a          n/a
%      .optimalFiberLength   cm           cadaver study + heuristic for groups 
%      .PSCA                 cm^2         cadaver study + heuristic for groups 
%      .peakForce            N            .PCSA * 2 * 30.5 N/cm^2
%      .pennationAngle       deg          cadaver study
%      .tendonSlackLength    cm           from kinematic model
%
%        Since some of these parameters are likely to vary a lot from
%        person to person you should read the paper to determine if these
%        values are appropriate for your application.
%%

arnold2010LegMuscleArchitectureFiles = { ...
        'arnold2010LegMuscleArchitectureAbbreviation.txt',       ... 
        'arnold2010LegMuscleArchitectureNames.txt',              ...
        'arnold2010LegMuscleArchitectureOptimalFiberLength.txt', ...
        'arnold2010LegMuscleArchitecturePCSA.txt',               ...
        'arnold2010LegMuscleArchitecturePeakForce.txt',          ...
        'arnold2010LegMuscleArchitecturePennationAngle.txt',     ... 
        'arnold2010LegMuscleArchitectureTendonSlackLength.txt' };


typeString = 1;
typeNumber = 0;

arnold2010LegMuscleArchitecture = [];

arnold2010LegMuscleArchitecture.abbrevation = readSingleColumnTextData(...
                                     arnold2010LegMuscleArchitectureFiles{1},...
                                        typeString);
                                    
arnold2010LegMuscleArchitecture.names       = readSingleColumnTextData(...
                                     arnold2010LegMuscleArchitectureFiles{2},...
                                        typeString);

arnold2010LegMuscleArchitecture.optimalFiberLength = readSingleColumnTextData(...
                                     arnold2010LegMuscleArchitectureFiles{3},...
                                        typeNumber);

arnold2010LegMuscleArchitecture.PSCA = readSingleColumnTextData(...
                                     arnold2010LegMuscleArchitectureFiles{4},...
                                        typeNumber);
                                    
arnold2010LegMuscleArchitecture.peakForce = readSingleColumnTextData(...
                                     arnold2010LegMuscleArchitectureFiles{5},...
                                        typeNumber);                                    
                                    
arnold2010LegMuscleArchitecture.pennationAngle = readSingleColumnTextData(...
                                     arnold2010LegMuscleArchitectureFiles{6},...
                                        typeNumber);
                                    
arnold2010LegMuscleArchitecture.tendonSlackLength = readSingleColumnTextData(...
                                     arnold2010LegMuscleArchitectureFiles{7},...
                                        typeNumber); 
if( unitsMKSN == 1)
    arnold2010LegMuscleArchitecture.optimalFiberLength = ...
        arnold2010LegMuscleArchitecture.optimalFiberLength./100;
    
    arnold2010LegMuscleArchitecture.PSCA = ...
        arnold2010LegMuscleArchitecture.PSCA ./ (100*100);
    
    arnold2010LegMuscleArchitecture.pennationAngle = ...
        arnold2010LegMuscleArchitecture.pennationAngle .* (pi/180);
    
    arnold2010LegMuscleArchitecture.tendonSlackLength = ...
        arnold2010LegMuscleArchitecture.tendonSlackLength ./ 100;
end
                                    
                                    