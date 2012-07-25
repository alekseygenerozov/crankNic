Liu & Schapiro 2010 Solution
============================

This recreates the steady state solution from L&S 2010, and is
a good way to debug the torque profile and our solver's correct
handling of it

Optimal Resoltion
-----------------

The following table gives optimal lambda values, for a
certain # of grid cells. They assume lMin is at the ISCO
and lOut is at 300.0. They optimize at l_a = 10.0

N			lambda		dl @ l_a
------------------------
300		1.01709		0.16
400		1.01280		0.12
500		1.01023		0.10
600		1.00852		0.08
700		1.00730		0.07
800		1.00639		0.06 
900		1.00568		0.05 
1000	1.00511		0.05 
