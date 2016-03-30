# Reentry_Apollo_11
The Reentry problem for Apollo 11 type reentry vehicles.

Based on the German version written by Prof. Folkmar Bornemann, TUM, WiSe 2002/2003
for the lectures on Numerical Mathematics III.

Reentry angle <=> Flight path angle 

Questions:
1. Change the critical point in the Maple script which restricts 'u' unnecessarily to the interval [-pi/2, pi/2]
   What is the condition of the minimum principle now?
2. What will be the maneuver time and the value of J to get the
   reentry angle of -6.5 degrees?
3. Which reentry angle optimizes J? How high is the gain over that studied in the lecture situation of -8.1 Degrees?

If you think hard, so you need to solve (1) - (3)
no more than a handful of changes to the Maple script and
make Matlab files.

Answers:
To (1)
	The limitation of the range of u is avoided
    by
    
    U [solve]: = arctan (+ lambda [2] * c [3], + lambda [1] * x [1] * c [2]):

    or

    U [solve]: = arctan (-Lambda [2] * c [3], - lambda [1] * x [1] * c [2]):

    The Pontryagin's Mimimum principle is always satisfied for the second option as compared to never for the first.
    So we choose the second option.

    Therefore, it also eliminates the condition lambda [1] < 0
    How to choose now the useful starting values ​​for the adjoint variables?

To (2)
    The confirmation of the solution of the auxiliary problem to the control parameter u with values ​​in the range [-pi/2, pi/2] should be removed.
    The auxiliary problem should be solved with backward shooting method.

To (3)
    We need for gamma at t = 0, the natural boundary condition.
    We must choose, the associated adjoint variable to be zero at t = 0.
    This provides the following changes:

    in reentry_bc.jl

        r[2] = xa[5];

    in reentry.jl

        before the solution of the optimal control problem,
        force natural boundary condition

        x[5,1] = 0;

        A few more multiple shooting nodes will not hurt.

That's it.