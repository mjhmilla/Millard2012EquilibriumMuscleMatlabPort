# Millard2012EquilibriumMuscleMatlabPort
A Matlab port of the Millard2012EquilibriumMuscle that is a part of OpenSim

Introduction
============

This repository contains a Matlab implementation of the muscle model described
in this paper:
```
Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). Flexing computational 
muscle: modeling and simulation of musculotendon dynamics. Journal of 
biomechanical engineering, 135(2), 021005.
```
If you make use of this Matlab implementation please cite this paper in your work. 

This is an implementation that has been made for teaching, thus it is been coded so that it is relatively easy to read, but is slow. If you want to dig into the details of how this model works, and are unfamiliar with C++, then this implementation is a good place to start. If you need something faster please use the C++ implementation provided by OpenSim: it is around 100x faster than this Matlab implementation. Furthermore, you can access this implementation using the python/matlab interfaces that OpenSim provides: not knowing C++ is not an excuse that should stop you from using the official implementation. 

Finally, this script can be made to run in octave. To do so, simply change the call to ode45/ode15s on line 224 of runMillard2012ComputationalBenchmark.m  to use octaves lsode integrator.


Matlab vs. OpenSim's C++ implementation
============

The Matlab implementation contains at least one difference with the OpenSim
model: the force velocity curve in the Matlab version is fitted to Hill's 
hyperbola. In the OpenSim implementation the parameters concentric and eccentric
curviness are exposed to the user and it is up to them to choose these values
appropriately. There may be other minor differences between the implementations
as it has been a while since I've done a side-by-side numerical comparison of 
the two implementations. These differences will be restricted to parameter 
values: the underlying mathematics are identical.


Getting Started
============
1. Run main.m from Matlab: this should produce plots that show all of the normalized muscle curve values, and first 2 derivatives.
   - Now in main.m set
```
flag_plotNormMuscleCurves = 0;
```
unless you want these plots popping up all the time.

2. Run a sinusoidal-stretch constant activation simulation of a muscle: set one or more of these flags to 1 in main.m. For example:
```
flag_runRigidBench               = 0;
flag_runClassicElasticBench      = 0;
flag_runDampedFiberElasticBench  = 1;
```
These simulations may take a few minutes, so be patient. Perhaps read the rest of the text below while its running.

Basic Code Exploration
============
1. In main.m near the top have a look at
```
flag_useArnold2010SoleusArchitecture = 0;
muscleAbbrArnold2010                 = 'soleus';
```
Setting the flag to 1 will mean the muscle architecture will match the soleus muscle that appears in Arnold et al. You can also choose from a variety of other architectures using an abbreviation from the file 'arnold2010LegMuscleArchitectureAbbreviation.txt'.
```
Arnold, E. M., Ward, S. R., Lieber, R. L., & Delp, S. L. (2010). A model of 
the lower limb for analysis of human movement. Annals of biomedical 
engineering, 38(2), 269-279.
```
2. Also in main.m near the top 
```
normPathStretch = 0.5; 
cycleTime       = 1; 
```
These two variables set the length change and period of the sinusoidal stretch.

Advanced Code Exploration
============

1. Examine 
```
calcMillard2012DampedEquilibriumMuscleInfo.m
```
This function evaluates the state derivative and a large number of variables of interest.

2. Examine
```
calcInitialMuscleState.m
```
which is a function that initializes the state of the muscle.

3. The code that creates the Bezier curves for each muscle are found in
```
createCurveIntegralStructure.m
createDefaultNormalizedMuscleCurves.m
createFiberActiveForceLengthCurve.m
createFiberForceLengthCurve.m
createFiberForceVelocityCurve2018.m
createTendonForceLengthCurve.m
```
  4. Code that evaluates the Bezier curves
```
calc1DBezierCurveValue.m
calc5thOrderInterp.m
calcBezierYFcnXCurveSampleVector.m
calcBezierYFcnXDerivative.m
```
  5. Code that evaluates the pennation model
```
calcFixedWidthPennatedFiberKinematicsAlongTendon.m
calcFixedWidthPennatedFiberKinematics.m
calcFixedWidthPennatedFiberMinimumLength.m
calcFixedWidthPennationDalphaDlce.m
```
  6. Code that evaluates the state derivative of activation
```
calcFirstOrderActivationDerivative.m
```
  8. The script that runs the sinusoidal-stretch constant activation simulations.
```
runMillard2012ComputationalBenchmark.m
```
Details
=========
<b>Authors</b> M.Millard

<b>Revision Notes</b>  

- 2011-2012: Partial implementation
- March 2015: Giant overhaul in preparation for a mini-course on muscle modelling
- May 2018: Small improvements to the force-velocity curve
- September 2019: further sprucing up

<b>License</b>
```
Licensed under the Apache License, Version 2.0 (the "License"); you may    
not use this file except in compliance with the License. You may obtain a  
copy of the License at http://www.apache.org/licenses/LICENSE-2.0. 
```

