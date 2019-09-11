%/* -------------------------------------------------------------------------- *
% *                           singleMuscleBench.cpp                            *
% * -------------------------------------------------------------------------- *
% * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
% * See http://opensim.stanford.edu and the NOTICE file for more information.  *
% * OpenSim is developed at Stanford University and supported by the US        *
% * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
% * through the Warrior Web program.                                           *
% *                                                                            *
% * Copyright (c) 2005-2012 Stanford University and the Authors                *
% * Author(s): Matthew Millard                                                 *
% *                                                                            *
% * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
% * not use this file except in compliance with the License. You may obtain a  *
% * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
% *                                                                            *
% * Unless required by applicable law or agreed to in writing, software        *
% * distributed under the License is distributed on an "AS IS" BASIS,          *
% * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
% * See the License for the specific language governing permissions and        *
% * limitations under the License.                                             *
% * -------------------------------------------------------------------------- */
%
% Derivative work
% Date      : March 2015
% Authors(s): Millard
% Updates   : Ported to code to Matlab which included some major rework
%
% If you use this code in your work please cite this paper
%
%  Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). 
%    Flexing computational muscle: modeling and simulation of 
%    musculotendon dynamics. Journal of biomechanical engineering, 
%    135(2), 021005.
%%
function benchRecord = ...
    runMillard2012ComputationalBenchmark(calcMuscleInfoFcn,...   
                                         calcInitialMuscleStateFcn,...
                                         benchConfig,...
                                         figBasicInfo,...
                                         figEnergyInfo,...
                                         figPowerInfo)
%%
% This function runs the computational benchmark that is described in this
% paper. 
%
%    Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). 
%    Flexing computational muscle: modeling and simulation of 
%    musculotendon dynamics. Journal of biomechanical engineering, 
%    135(2), 021005.
%
% @param calcMuscleInfoFcn: a function handle to a muscle function
%                      like (calcMillard2012DampedEquilibriumMuscleInfo.m) 
%                      and takes arguments of activation, pathState, and 
%                      muscle state
%
% @param calcInitialMuscleStateFcn: a function handle to a function that
%                                 takes arguments of activation, path
%                                 state, calcMuscleInfoFcn, and
%                                 a struct initConfig (see
%                                 calcInitialMuscleState.m for details)
%                                 and returns a muscle state that satisfies
%                                 the equilibrium equation
%
% @param benchConfig: a structure that configures this constant activation
% sinusoidal stretch benchmark simulation:
%
%   benchConfig.npts                 : number of points to evaluate the
%                                      model at during a 1 second constant 
%                                      activation sinusoidal stretch 
%                                      simulation
%
%   benchConfig.numberOfMuscleStates : 0 for rigid tendon, 
%                                      1 for elastic tendon
%
%   benchConfig.minimumActivation    : 0 unless the classic elastic tendon
%                                      model is being used, then it should 
%                                      be > 0, and probably 0.05
%
%   benchConfig.name                 : short musclemodel name e.g.
%                                      'soleusRigidTendon'   
%
% @param figBasicInfo   : handle to an empty figure
%
% @param figEnergyInfo  : handle to an empty figure
%
% @param figPowerInfo   : handle to an empty figure
%
% @return benchRecord: A structure containing a set of matricies, where
%                      each column represents the data of 1 
%                      constant-activation sinusoidal-stretch simulation.
% 
%     benchRecord.activation                  
%     benchRecord.cpuTime                     
%     benchRecord.normfiberForceAlongTendon   
%     benchRecord.normFiberLength              
%     benchRecord.pennationAngle              
%     benchRecord.normFiberVelocity           
%     benchRecord.pennationAngVelocity        
%     benchRecord.fiberStiffnessAlongTendon   
%     benchRecord.tendonStiffnessAlongTendon  
%     benchRecord.muscleStiffness             
%     benchRecord.pathLength                  
%     benchRecord.pathVelocity                
% 
%     benchRecord.dSystemEnergyLessWork       
%     benchRecord.systemEnergyLessWork        
%     benchRecord.tendonPotentialEnergy       
%     benchRecord.fiberPotentialEnergy        
%     benchRecord.fiberActiveWork             
%     benchRecord.dampingWork                 
%     benchRecord.boundaryWork                
% 
%     benchRecord.tendonPower                 
%     benchRecord.fiberParallelElementPower   
%     benchRecord.fiberActivePower            
%     benchRecord.dampingPower                
%     benchRecord.boundaryPower               

%
%%
                                     


npts    = benchConfig.npts;

tmin    = benchConfig.tspan(1);
tmax    = benchConfig.tspan(2);

tV      = [tmin:((tmax-tmin)/(npts-1)):tmax];

aV      = benchConfig.activationVector;
if(min(aV) < benchConfig.minimumActivation)
  [val,idx] = min(aV);
  aV(idx) = benchConfig.minimumActivation;
end


pathFcn = benchConfig.pathFcn;
        
benchRecord = [];
benchRecord.activation                  = aV;
benchRecord.cpuTime                     = zeros(size(aV));
benchRecord.normfiberForceAlongTendon   = zeros(npts,length(aV));
benchRecord.normFiberLength             = zeros(npts,length(aV));    
benchRecord.pennationAngle              = zeros(npts,length(aV));
benchRecord.normFiberVelocity           = zeros(npts,length(aV));
benchRecord.pennationAngVelocity        = zeros(npts,length(aV));
benchRecord.fiberStiffnessAlongTendon   = zeros(npts,length(aV));
benchRecord.tendonStiffnessAlongTendon  = zeros(npts,length(aV));
benchRecord.muscleStiffness             = zeros(npts,length(aV));
benchRecord.fiberVelocity               = zeros(npts,length(aV));
benchRecord.fiberVelocityAlongTendon    = zeros(npts,length(aV));
benchRecord.tendonVelocity              = zeros(npts,length(aV));
benchRecord.pathLength                  = zeros(npts,length(aV));
benchRecord.pathVelocity                = zeros(npts,length(aV));

benchRecord.dSystemEnergyLessWork       = zeros(npts,length(aV));
benchRecord.systemEnergyLessWork        = zeros(npts,length(aV));
benchRecord.tendonPotentialEnergy       = zeros(npts,length(aV));
benchRecord.fiberPotentialEnergy        = zeros(npts,length(aV));
benchRecord.fiberActiveWork             = zeros(npts,length(aV));
benchRecord.dampingWork                 = zeros(npts,length(aV));
benchRecord.boundaryWork                = zeros(npts,length(aV));

benchRecord.tendonPower                 = zeros(npts,length(aV));
benchRecord.fiberParallelElementPower   = zeros(npts,length(aV));
benchRecord.fiberActivePower            = zeros(npts,length(aV));
benchRecord.dampingPower                = zeros(npts,length(aV));
benchRecord.boundaryPower               = zeros(npts,length(aV));


for i = 1:1:length(aV)

    
    activationFcn = @(t)calcConstantState(t,aV(i));

    muscleStateV = [];
    
    t0 = 0;
    cpuTime = 0;
    
       
               
    dfcn = @(argt,argState)...
            calcPrescribedMusculotendonStateDerivativeWrapper(...
                              argt,...
                              argState,...                                                      
                              pathFcn,...
                              activationFcn,...
                              calcMuscleInfoFcn);

    muscleState0 = [];
    if(benchConfig.numberOfMuscleStates ~= 0)     
        initConfig.iterMax = 100;
        initConfig.tol     = 1e-8;
        initConfig.useStaticFiberSolution = 0;
        initSoln = calcInitialMuscleStateFcn(aV(i),...
                                          pathFcn(0),...                                          
                                          calcMuscleInfoFcn,...
                                          initConfig);

        if(initSoln.converged == 0 && initSoln.isClamped == 0)
            %%
            %If we're here then we've been unlucky enough to be in the
            %negative stiffness region of the active force length curve and
            %we will have to settle for an initial condition where the
            %velocity of the fiber is 0 ... which often leads to force
            %transients at the beginning of the simulation.
            %%
            initConfig.useStaticFiberSolution = 1;
            initSoln = calcInitialMuscleStateFcn(aV(i),...
                                          pathFcn(0),...
                                          calcMuscleInfoFcn,...
                                          initConfig);

            assert(initSoln.converged == 1 || initSoln.isClamped == 1,...
                   'Failed to bring the muscle to a valid initial solution');
        end

        muscleState0 = initSoln.muscleState(:);

    end

    options = odeset('RelTol',benchConfig.relTol,...
             'AbsTol',benchConfig.absTol,...
             'Stats','off');
    t0 = tic;
    disp([' ',num2str(i),' of ',num2str(length(aV))]);
    [xe ye] = ode15s(dfcn,tV,[muscleState0(:);0;0;0],options);                                                   
    cpuTime = toc(t0);
    workOfBoundary     = zeros(size(xe));
    workOfActiveFiber = zeros(size(xe));
    workOfDamping     = zeros(size(xe));

    n = benchConfig.numberOfMuscleStates;
    if(n >= 1)            
        muscleStateV = ye(:,1:1:n);
    end

    workOfBoundary    = ye(:,n+1);
    workOfActiveFiber = ye(:,n+2);
    workOfDamping     = ye(:,n+3);
                        
    if(benchConfig.numberOfMuscleStates == 0) 
        t0 = tic;
    end

    T0V0 = 0;
    for j=1:1:length(tV)
        
        muscleState     = NaN;
        
        if(benchConfig.numberOfMuscleStates ~= 0)
           muscleState = muscleStateV(j,:)'; 
        end
        
        
        activationState = activationFcn(tV(j));
        pathState       = pathFcn(tV(j));

        dlp = pathState(1);
        lp  = pathState(2);
        
        mtInfo = calcMuscleInfoFcn(activationState,...                                   
                                   pathState,...
                                   muscleState);

        tendonForce = mtInfo.muscleDynamicsInfo.tendonForce;
               
        lceN    = mtInfo.muscleLengthInfo.normFiberLength;
        alpha   = mtInfo.muscleLengthInfo.pennationAngle;
        dlceN   = mtInfo.fiberVelocityInfo.normFiberVelocity;
        dalpha  = mtInfo.fiberVelocityInfo.pennationAngularVelocity;
                
        fNAT    = mtInfo.muscleDynamicsInfo.normFiberForce*cos(alpha);
        kFiberAT= mtInfo.muscleDynamicsInfo.fiberStiffnessAlongTendon;
        kTendon = mtInfo.muscleDynamicsInfo.tendonStiffness;
        kMuscle = mtInfo.muscleDynamicsInfo.muscleStiffness;

        %%
        %force and kinematic information
        %%
        
        benchRecord.normfiberForceAlongTendon(j,i)   = fNAT;
        benchRecord.normFiberLength(j,i)             = lceN;    
        benchRecord.pennationAngle(j,i)              = alpha;
        benchRecord.normFiberVelocity(j,i)           = dlceN;
        benchRecord.pennationAngVelocity(j,i)        = dalpha;
        benchRecord.fiberStiffnessAlongTendon(j,i)   = kFiberAT;
        benchRecord.tendonStiffnessAlongTendon(j,i)  = kTendon;
        benchRecord.muscleStiffness(j,i)             = kMuscle;
        benchRecord.pathLength(j,i)                  = lp;
        benchRecord.pathVelocity(j,i)                = dlp;
        
        benchRecord.fiberVelocity(j,i) = mtInfo.fiberVelocityInfo.fiberVelocity;
        benchRecord.fiberVelocityAlongTendon(j,i) = mtInfo.fiberVelocityInfo.fiberVelocityAlongTendon;
        benchRecord.tendonVelocity(j,i) = mtInfo.fiberVelocityInfo.tendonVelocity;

        %%
        %Energetics & power information
        %%
        
        fiberVelocity    = mtInfo.fiberVelocityInfo.fiberVelocity;
        activeFiberForce = mtInfo.muscleDynamicsInfo.activeFiberForce;
        tendonForce      = mtInfo.muscleDynamicsInfo.tendonForce;
        dampingForce     = mtInfo.muscleDynamicsInfo.dampingForces;

        benchRecord.tendonPower(j,i)               = mtInfo.muscleDynamicsInfo.tendonPower;
        benchRecord.fiberParallelElementPower(j,i) = mtInfo.muscleDynamicsInfo.fiberParallelElementPower;
        benchRecord.fiberActivePower(j,i)          = mtInfo.muscleDynamicsInfo.fiberActivePower;
        benchRecord.dampingPower(j,i)              = mtInfo.muscleDynamicsInfo.dampingPower;
        benchRecord.boundaryPower(j,i)             = mtInfo.muscleDynamicsInfo.boundaryPower;

        boundaryPower       = mtInfo.muscleDynamicsInfo.boundaryPower;
        activeFiberPower    = mtInfo.muscleDynamicsInfo.fiberActivePower;
        dampingPower        = mtInfo.muscleDynamicsInfo.dampingPower;

        tendonPower               = mtInfo.muscleDynamicsInfo.tendonPower;
        fiberParallelElementPower = mtInfo.muscleDynamicsInfo.fiberParallelElementPower;
        
        dTpVmW = - tendonPower...
                 - fiberParallelElementPower ...
               - activeFiberPower ...
               - dampingPower ...
               - boundaryPower;
        
        if(j==1)
            T0V0 = mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy ...
                +  mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy;
        end
        
        TpVmW   =  mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy ...
                 + mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy ...
                 - workOfActiveFiber(j) ...
                 - workOfDamping(j) ...
                 - workOfBoundary(j) ...
                 - T0V0;
                
        benchRecord.systemEnergyLessWork(j,i)        = TpVmW;
        benchRecord.dSystemEnergyLessWork(j,i)       = dTpVmW;
        benchRecord.tendonPotentialEnergy(j,i)       = ...
            mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy;
        benchRecord.fiberPotentialEnergy(j,i)        = ...
            mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy;
        benchRecord.fiberActiveWork(j,i)             = workOfActiveFiber(j);
        benchRecord.dampingWork(j,i)                 = workOfDamping(j);
        benchRecord.boundaryWork(j,i)                = workOfBoundary(j);
        
        
        
    end       

    if(benchConfig.numberOfMuscleStates == 0)
        cpuTime = toc(t0);
    end
    
    benchRecord.cpuTime(i) = cpuTime;
    %clr01 = aV(i).*[1,0,0] + (1-aV(i)).*[0,0,1];
    clr01 = benchConfig.color.*aV(i) + [0.5,0.5,0.5].*(1-aV(i));
    
    if(isempty(figBasicInfo) == 0) 
        figure(figBasicInfo);


        subplot(2,2,1)
            plot(tV, benchRecord.normFiberLength(:,i),...
                 'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Norm. Length (m/lceOpt)');
            title([benchConfig.name, ' Norm. Fiber Length']);                
            hold on;   


        subplot(2,2,2)
            plot(tV, benchRecord.normfiberForceAlongTendon(:,i),...
                 'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('norm. Force (N/fiso)');
            title([benchConfig.name, ' Norm. Fiber Force Along Tendon']);                     
            hold on;

            set(gca,'XTick',[0,0.25,0.5,0.75,1].*tmax);
            set(gca,'YTick',[0,0.25,0.5,0.75,1,1.25]);
            
        subplot(2,2,3)
            plot(tV, benchRecord.fiberVelocityAlongTendon(:,i).*1000,...
                 '--','Color',clr01.*0.5+[0.5,0.5,0.5]);               
            hold on;
            plot(tV, benchRecord.fiberVelocity(:,i).*1000,...
                 'Color',clr01);
            hold on;
            legend('Along tendon','Along fiber','Location','SouthEast');
            xlabel('Time (s)'); 
            ylabel('Velocity (mm/s)');
            title([benchConfig.name, 'Fiber Velocity']);                      
            hold on;
            
            %set(gca,'YTick',[-1,0,1]);
            
        subplot(2,2,4)
            %plot(tV,...
            %     benchRecord.muscleStiffness(:,i),...
            %     'Color',clr01);
            %xlabel('Time (s)'); 
            %ylabel('Stiffness (N/m)');
            %title([benchConfig.name, ' Musculotendon Stiffness']);                
            plot(tV, benchRecord.tendonVelocity(:,i).*1000,...
                 'Color',clr01);
            hold on;
            xlabel('Time (s)'); 
            ylabel('Velocity (mm/s)');
            title([benchConfig.name, 'Tendon Velocity']);                      
            hold on;
            
            
            hold on;   


        
    end
    
    if(isempty(figEnergyInfo) == 0)
        figure(figEnergyInfo);
        subplot(2,3,1);
            plot(tV, benchRecord.tendonPotentialEnergy(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Energy (J)');
            title([benchConfig.name, ' Tendon Potential Energy']);                      
            hold on;
            
        subplot(2,3,2);
            plot(tV, benchRecord.fiberPotentialEnergy(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Energy (J)');
            title([benchConfig.name, ' Fiber Potential Energy']);                      
            hold on;
                
            
        subplot(2,3,3);
            plot(tV, benchRecord.systemEnergyLessWork(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Energy (J)');
            title([benchConfig.name, ' T+V-W-(V0+T0): should be 0']);                      
            hold on;
            
        subplot(2,3,4);
            plot(tV, benchRecord.fiberActiveWork(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Energy (J)');
            title([benchConfig.name, ' Active Fiber Work']);                      
            hold on;
            
        subplot(2,3,5);
            plot(tV, benchRecord.dampingWork(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Energy (J)');
            title([benchConfig.name, ' Damping Work']);                      
            hold on;

        subplot(2,3,6);
            plot(tV, benchRecord.boundaryWork(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Energy (J)');
            title([benchConfig.name, ' Boundary Work']);                      
            hold on;
    end
    
    if(isempty(figPowerInfo) == 0)
        figure(figPowerInfo);
        subplot(2,3,1);
            plot(tV, benchRecord.tendonPower(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Power (W)');
            title([benchConfig.name, ' Tendon Potential Power']);                      
            hold on;
            
        subplot(2,3,2);
            plot(tV, benchRecord.fiberParallelElementPower(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Power (W)');
            title([benchConfig.name, ' Fiber Potential Power']);                      
            hold on;
                
        subplot(2,3,3);
            plot(tV, benchRecord.dSystemEnergyLessWork(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Power (W)');
            title([benchConfig.name, ' d/dt T+V-W: should be 0']);                      
            hold on;
            
        subplot(2,3,4);
            plot(tV, benchRecord.fiberActivePower(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Power (W)');
            title([benchConfig.name, ' Active Fiber Power']);                      
            hold on;
            
        subplot(2,3,5);
            plot(tV, benchRecord.dampingPower(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Power (W)');
            title([benchConfig.name, ' Damping Power']);                      
            hold on;

        subplot(2,3,6);
            plot(tV, benchRecord.boundaryPower(:,i),...
                'Color',clr01);
            xlabel('Time (s)'); 
            ylabel('Power (W)');
            title([benchConfig.name, ' Boundary Power']);                      
            hold on;
    end
           
end



                                         
                                     