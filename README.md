# Embodying Control in Soft Multistable Robots from Morphofunctional Co-design

Soft robots are distinguished by their flexibility and adaptability, allowing them to perform nearly impossible tasks for rigid robots. However, controlling their behavior is challenging due to their nonlinear material response and infinite degrees of freedom. A potential solution to these challenges is to discretize the infinite-dimensional configuration space into a finite but sufficiently large number of functional modes with programmed dynamics. We present a strategy for co-designing the desired tasks and morphology of pneumatically actuated soft robots with multiple encoded stable states and dynamic responses. Our approach introduces a general method to capture the soft robots’ response using an energy-based analytical model, the parameters of which are obtained using Recursive Feature Elimination. The resulting lumped-parameter model facilitates inverse co-design of the robot’s morphology and planned tasks by embodying specific dynamics upon actuation. We illustrate our approach’s ability to explore the configuration space by co-designing kinematics with optimized stiffnesses and time responses to obtain robots capable of classifying the size and weight of objects and displaying adaptable locomotion with minimal feedback control. This strategy offers a framework for simplifying the control of soft robots by exploiting the nonlinear mechanics of multistable structures and embodying mechanical intelligence into soft material systems.

## Key Features:

+ CAD_Files -> STL Files of all the optimized geometries for the Dome Phalanx Gripper and Dome Phalanx Walker.
+ DPG_Static -> MATLAB Code to simulate the static behavior of the optimized Dome Phalanx Gripper. 
+ Dynamic_Model_DPF ->  Model to simulate the dynamic behavior of the Dome Phalanx Finger.
+ Dynamic_DPW   -> Model to simulate the dynamic behavior of the Dome Phalanx Walker leg.

## Getting Started:
+ DPG_Static -> To get the different shapes of the robot (for different stable configurations), run file DPG_Static_Response.m
+ Dynamic_Model_DPF ->  To get the energy vs. time plot and an animation of the dynamic behavior of the DPF,  run the file DPF_Dyamic.m
+ Dynamic_DPW   -> To obtain the leg tip displacement under different dome arrangements. Run file Dynamic_DPW.m.

For details, see the paper:
Embodying Control in Soft Multistable Robots from Morphofunctional Co-design. Juan C. Osorio, Jhonatan S. Rincon, Harith Morgan, Andres F. Arrieta