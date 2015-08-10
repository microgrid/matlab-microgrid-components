# matlab-microgrid-components
A simulation to find the optimized sizes of microgrid components (PV and battery) constrained by a certain acceptable loss of load percentage and by budget. This simulation is written by Stefano Mandelli and expanded by Håkon Duus.


## todo list (possibilities)
- [ ] optimize the minimal State of Charge of the battery, SoC_min, w.r.t. cost since SoC_min influences the battery lifetime and depends on the battery type. Then SoC_min will be an output value instead of input. 
N.B. CyclesToFailure() also depends on the battery and must be changed then.
- [ ] Play with the value of LLP_var. Is now very small. 
- [ ] Could include X days of autonomy in LoadCurve (randomly placed Y times)
- [ ] A. Algorithm e.g. start with graphical point of view
- [ ] B. include error bar (and shifts) in load curve
- [ ] C. Rewrite in Python
- [ ] D. spikes. Play with stepsize and make a test.
- [ ] E. visual interface of graphs (Honorat)
- [ ] F. Think about dependance on square meters. x, t, m2.

## todo list (tests)
- [ ] give error if SoC < 20% (It may never be zero)
- [ ] check that range of PV and battery is divisible by PV_step and Batt_step
- [ ] evt. check that size of dimensions is larger than one s.t. taken over right dimension?

## todo list (assumptions)
- [ ] coeff_cost_BoSeI is now 20% of total cost. Is this reliable or can it vary a lot? May have large impact on costs.
