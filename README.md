# matlab-microgrid-components
A simulation to find the optimized sizes of microgrid components (PV and battery) constrained by a certain acceptable loss of load percentage and by budget. This simulation is written by Stefano Mandelli and expanded by Håkon Duus.


## todo list (possibilities)
- [ ] optimize the minimal State of Charge of the battery, SoC_min, w.r.t. cost since SoC_min influences the battery lifetime and depends on the battery type. Then SoC_min will be an output value instead of input. 
N.B. CyclesToFailure() also depends on the battery and must be changed then.

## todo list (tests)
- [ ] give error if SoC < 20% (It may never be zero)

## todo list (assumptions)
- [ ] 
