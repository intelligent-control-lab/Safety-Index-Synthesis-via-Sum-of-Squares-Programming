# Safety-Index-Synthesis-via-Sum-of-Squares-Programming

By Weiye Zhao, Tairan He, Tianhao Wei, Simin Liu and Changliu Liu
(the code of ``Safety Index Synthesis via Sum of Squares Programming'', American Control Conference, 2023)

### Introduction
**Safety Index Synthesis Using Sum of Squares Programming** is currently implemented under two experimental settings, including (i) varying degrees-of-freedom planar manipulator collision avoidance with limited acceleration, and (ii) unicycle collision avoidance with limited acceleration and angular velocity

### Requirements: software

0.	MATLAB 2016a or later.

### Requirements: hardware

CPU, Windows 7 or later, MAC OS.

### Demo
0.	Run `experiments/unicycle_SI_optimization.m` to synthesis an optimal safety index for unicycle collision avoidance. the acceleration limits are [-1,1], and the angular velocity limits are [-0.1,0.1]. velocity limits are [0,1], and heading angle range is [0, pi/2].
0.	Run `experiments/arm_ndof_SI_optimization.m` to synthesis an optimal safety index for **n** degrees-of-freedom planar manipulator collision avoidance. (The safety index is searched within a limited range, and feasible safety index will be collected and optimal safety index can then be chosen) The acceleration limits for all joints are [-1,1]. The state space limits are specified in the code.


### Note
0. For **Planar Manipulator Experiments**, it is recommend to first try `experiments/arm_1dof_SI_optimization.m`, which is a simplest version where only 1 link planar robot is considered. The explicit SOSP derivative is elaborated in the [paper section VII B](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=10156463), the only difference is that code in this repo considers a stricter condtion *non-empty set of safe control for all possible states*, meaning only f1,f2,f3,f4 are included in (18). 
0. For **Unicycle Experiment**, it is recommend to read [this paper](https://openreview.net/forum?id=UGp6FDaxB0f), which analytically derived the safety index synthesis rule for unicycle. According to the (11) of Appendex C.1, when $\cos(\alpha) > 0$, safety index should satisfy the following condition:
$$
\forall (a,v), \exists (a,w), \text{ s.t. } -\frac{a}{v} + \tan(\alpha) w \geq \frac{n d^{n-1}}{k},
$$
which indicates setting $n=1$, and take the maximizer of acceleration and angular velocity, the following condition should hold
$$
\forall (a,v), k + \tan(\alpha)vk - v \geq 0.
$$

Therefore, we can then following the paper and construct the refute problem and associated SOSP as documented in `experiments/unicycle_SI_optimization.m`.

