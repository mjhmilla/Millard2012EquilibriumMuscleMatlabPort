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

clc;
close all;
clear all;

flag_useArnold2010SoleusArchitecture = 1;
muscleAbbrArnold2010                 = 'soleus';

flag_plotNormMuscleCurves            = 1;
flag_updateCurves                    = 1;

flag_runRigidBench               = 0;
flag_runClassicElasticBench      = 0;
flag_runDampedFiberElasticBench  = 0;


normPathStretch = 0.5;
cycleTime       = 1;

%%
%Parameters that apply to all muscles
%%
maximumNormalizedFiberVelocity = 10; %in units of norm fiber lengths/second
maximumPennationAngle          = 89*(pi/180); %if we go to 90 the 
                                              %classic formulation goes
                                              %singular.

%%
% Muscle color
%%
              
rigidTendonPlotColor              = [0,0,0];
classicElasticTendonPlotColor     = [1,0,0];
dampedFiberElasticTendonPlotColor = [0,0,1];

%%
%Get a copy of the default muscle curves 
%%
muscleAbbr  = [];
if(flag_useArnold2010SoleusArchitecture == 1)
    muscleAbbr = muscleAbbrArnold2010;
else
    muscleAbbr = 'compBench';
end

tendonStrainAtOneNormForce = -1; %If this is greater than 0 this value will
                                 %be used to make the tendon-force-length
                                 %curve. Otherwise the default of 0.049 is
                                 %taken.

normMuscleCurves = ...
    createDefaultNormalizedMuscleCurves(muscleAbbr,...
                                        tendonStrainAtOneNormForce,...
                                        flag_updateCurves,...
                                        flag_plotNormMuscleCurves);
                                    
                                    
%%
%Get a muscle and extract out its architecture information
%%

muscleName  = [];
fiso        = [];
lceOpt      = [];
alphaOpt    = [];
ltSlk       = [];


if(flag_useArnold2010SoleusArchitecture ==1)
    unitsMKSN = 1;
    arnold2010LegArch = getArnold2010LegMuscleArchitecture(unitsMKSN);

    idx =  getArnold2010MuscleIndex(muscleAbbrArnold2010,...
                          arnold2010LegArch.abbrevation);
    
    muscleName  = arnold2010LegArch.names{idx};
    fiso        = arnold2010LegArch.peakForce(idx);
    lceOpt      = arnold2010LegArch.optimalFiberLength(idx);
    alphaOpt    = arnold2010LegArch.pennationAngle(idx);
    ltSlk       = arnold2010LegArch.tendonSlackLength(idx);
else 
    muscleName                     = 'compBenchMillard2010';
    fiso    = 1;
    lceOpt  = 0.02;
    alphaOpt= 30*(pi/180);
    ltSlk   = 0.20;
end


muscleArch = [];
    muscleArch.name                = muscleName;
    muscleArch.abbr                = muscleAbbr;
    muscleArch.fiso                = fiso;
    muscleArch.optimalFiberLength  = lceOpt;
    muscleArch.maximumNormalizedFiberVelocity = ...
                maximumNormalizedFiberVelocity;
    muscleArch.pennationAngle      = alphaOpt;
    muscleArch.tendonSlackLength   = ltSlk;
        
    minimumActiveFiberNormalizedLength = ...
        normMuscleCurves.activeForceLengthCurve.xEnd(1);
        
    minFiberKinematics = calcFixedWidthPennatedFiberMinimumLength(...
                minimumActiveFiberNormalizedLength,...
                maximumPennationAngle,...
                muscleArch.optimalFiberLength,...
                muscleArch.pennationAngle);
    
    muscleArch.minimumFiberLength = ...
               minFiberKinematics.minimumFiberLength;
                                    
    muscleArch.minimumFiberLengthAlongTendon =...
               minFiberKinematics.minimumFiberLengthAlongTendon;
                         
    muscleArch.pennationAngleAtMinumumFiberLength = ...
               minFiberKinematics.pennationAngleAtMinimumFiberLength;



%%
%Run the computational benchmark described in:
%  Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). 
%    Flexing computational muscle: modeling and simulation of 
%    musculotendon dynamics. Journal of biomechanical engineering, 
%    135(2), 021005.
%%
benchConfig.npts   = 100;
benchConfig.relTol = 1e-6;
benchConfig.absTol = 1e-6;
benchConfig.activationVector = [0:0.1:1]';
benchConfig.tspan  = [0,1];

%%
%Create the path function
%%
lceOpt   = muscleArch.optimalFiberLength;
alphaOpt = muscleArch.pennationAngle;
ltSlk    = muscleArch.tendonSlackLength;
fiso     = muscleArch.fiso;

lp0     = lceOpt*cos(alphaOpt) + ltSlk;
lpDelta = lceOpt*normPathStretch;
omega   = 2*pi*(1/cycleTime);



pathFcn = @(t)calcSinusoidState(t,lp0,lpDelta, omega);                            

benchConfig.pathFcn = pathFcn;
benchConfig.tspan   = [0, cycleTime];

%%=========================================================================
%Rigid Tendon Model Benchmark
%%=========================================================================

if(flag_runRigidBench == 1)
    disp('Rigid-tendon model: sinusoidal stretch, constant activation');
    figRigidTendonBasic  = figure;
    figRigidTendonEnergy = figure;
    figRigidTendonPower  = figure;

    rigidConfig.useFiberDamping  = 0;
    rigidConfig.useElasticTendon = 0;
    rigidConfig.damping          = 0;
    rigidConfig.iterMax          = 100;
    rigidConfig.tol              = 1e-12;
    rigidConfig.minActivation    = 0.0;

    calcRigidTendonMuscleInfoFcn =...
        @(actState1,pathState2,mclState3)...
        calcMillard2012DampedEquilibriumMuscleInfo(  ...
                                    actState1,...
                                    pathState2, ... 
                                    mclState3,...                                                                           
                                    muscleArch,...
                                    normMuscleCurves,...
                                    rigidConfig);   



    benchConfig.numberOfMuscleStates = 0;
    benchConfig.minimumActivation = rigidConfig.minActivation;
    benchConfig.name = 'RT';
    benchConfig.color = rigidTendonPlotColor;
    
    calcInitalRigidMuscleState = [];

    benchRecordRigid = ...
        runMillard2012ComputationalBenchmark(calcRigidTendonMuscleInfoFcn,... 
                                             calcInitalRigidMuscleState ,...
                                             benchConfig,...
                                             figRigidTendonBasic,...
                                             figRigidTendonEnergy,...
                                             figRigidTendonPower);

    save('benchRecordRigid.mat','benchRecordRigid');
    saveas(figRigidTendonBasic,  'figRigidTendonBasic.fig','fig');
    saveas(figRigidTendonEnergy, 'figRigidTendonEnergy.fig','fig');
    saveas(figRigidTendonPower,  'figRigidTendonPower.fig','fig');   
end




%%=========================================================================
%Classic Elastic Tendon Model Benchmark
%%=========================================================================

if flag_runClassicElasticBench == 1
    disp('Classic elastic-tendon model: sinusoidal stretch, constant activation');
    figClassicElasticBasic  = figure;
    figClassicElasticEnergy = figure;
    figClassicElasticPower  = figure;

    classicElasticTendonConfig.useFiberDamping  = 0;
    classicElasticTendonConfig.useElasticTendon = 1;
    classicElasticTendonConfig.damping          = 0;
    classicElasticTendonConfig.iterMax          = 100;
    classicElasticTendonConfig.tol              = 1e-6;
    classicElasticTendonConfig.minActivation    = 0.05;

    calcClassicElasticTendonMuscleInfoFcn =...
        @(actState1,pathState2,mclState3)...
        calcMillard2012DampedEquilibriumMuscleInfo(  ...
                                    actState1,...
                                    pathState2, ... 
                                    mclState3,...                                                                           
                                    muscleArch,...
                                    normMuscleCurves,...
                                    classicElasticTendonConfig);   

    calcClassicElasticTendonInitialMuscleStateFcn = ...
        @(actState1,pathState2,calcMuscleInfo3, initConfig4) ...
            calcInitialMuscleState(actState1,...
                                   pathState2,...
                                   muscleArch,...
                                   calcMuscleInfo3,...
                                   initConfig4);


    benchConfig.numberOfMuscleStates = 1;
    benchConfig.minimumActivation    = ...
        classicElasticTendonConfig.minActivation;
    benchConfig.name = 'CE';
    benchConfig.color = classicElasticTendonPlotColor;


    benchRecordClassicElastic = ...
        runMillard2012ComputationalBenchmark(...
             calcClassicElasticTendonMuscleInfoFcn,...                                         
             calcClassicElasticTendonInitialMuscleStateFcn,...
             benchConfig,...
             figClassicElasticBasic,...
             figClassicElasticEnergy,...
             figClassicElasticPower);    


    save('benchRecordClassicElastic.mat','benchRecordClassicElastic');
    saveas(figClassicElasticBasic,'figClassicElasticBasic.fig','fig');
    saveas(figClassicElasticEnergy,'figClassicElasticEnergy.fig','fig');
    saveas(figClassicElasticPower,'figClassicElasticPower.fig','fig');
    
end
%%=========================================================================
%Damped Equilibrum Elastic Model Benchmark
%%=========================================================================

if flag_runDampedFiberElasticBench == 1
    disp('Damped-fiber elastic-tendon model: sinusoidal stretch, constant activation');
    disp('(Default elastic tendon model formulation in OpenSim)');
    figDampedFiberElasticBasic  = figure;
    figDampedFiberElasticEnergy = figure;
    figDampedFiberElasticPower  = figure;

    dampedFiberElasticTendonConfig.useFiberDamping  = 1;
    dampedFiberElasticTendonConfig.useElasticTendon = 1;
    dampedFiberElasticTendonConfig.damping          = 0.1;
    dampedFiberElasticTendonConfig.iterMax          = 100;
    dampedFiberElasticTendonConfig.tol              = 1e-6;
    dampedFiberElasticTendonConfig.minActivation    = 0.0;

    
    
    calcDampedFiberElasticTendonMuscleInfoFcn =...
        @(actState1,pathState2,mclState3)...
        calcMillard2012DampedEquilibriumMuscleInfo(  ...
                                    actState1,...
                                    pathState2, ... 
                                    mclState3,...                                                                           
                                    muscleArch,...
                                    normMuscleCurves,...
                                    dampedFiberElasticTendonConfig);   

    calcDampedFiberElasticTendonInitialMuscleStateFcn = ...
        @(actState1,pathState2,calcMuscleInfo3, initConfig4) ...
            calcInitialMuscleState(actState1,...
                                   pathState2,...
                                   muscleArch,...
                                   calcMuscleInfo3,...
                                   initConfig4);


    benchConfig.numberOfMuscleStates = 1;
    benchConfig.minimumActivation    = ...
        dampedFiberElasticTendonConfig.minActivation;
    benchConfig.name = 'DFE';
    benchConfig.color = dampedFiberElasticTendonPlotColor;

    benchRecordDampedFiberElasticTendon = ...
        runMillard2012ComputationalBenchmark(...
             calcDampedFiberElasticTendonMuscleInfoFcn,...                                         
             calcDampedFiberElasticTendonInitialMuscleStateFcn,...
             benchConfig,...
             figDampedFiberElasticBasic,...
             figDampedFiberElasticEnergy,...
             figDampedFiberElasticPower);  

    save('benchRecordDampedFiberElasticTendon.mat',...
        'benchRecordDampedFiberElasticTendon');
    saveas(figDampedFiberElasticBasic,...
        'figDampedFiberElasticBasic.fig','fig');
    saveas(figDampedFiberElasticEnergy,...
        'figDampedFiberElasticEnergy.fig','fig');
    saveas(figDampedFiberElasticPower,...
        'figDampedFiberElasticPower.fig','fig');
    
end
