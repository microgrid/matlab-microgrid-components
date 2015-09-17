# Matlab Microgrid Components
A simulation to find the optimized sizes of microgrid components (PV and battery) constrained by a certain acceptable loss of load percentage and by budget. This simulation is written by Stefano Mandelli and expanded by Håkon Duus.


## Todo list (possibilities)
- [ ] optimize the minimal State of Charge of the battery, SoC_min, w.r.t. cost since SoC_min influences the battery lifetime and depends on the battery type. Then SoC_min will be an output value instead of input. 
N.B. CyclesToFailure() also depends on the battery and must be changed then.
- [ ] Play with the value of LLP_var. Is now very small. 
- [ ] Could include X days of autonomy in LoadCurve (randomly placed Y times)
- [ ] A. Algorithm e.g. start with graphical point of view
- [ ] B. include error bar (and shifts) in load curve
- [ ] C. Rewrite in Python
- [ ] D. spikes. Play with stepsize and make a test.
- [ ] E. visual interface of graphs (Honorat)
- [ ] F. Think about dependance on square meters. x, t, m2. (m2 can be obtained from the rated power PV_size and the type of PV). Alternative is to write this out (and ask for these inputs explicitly).
- [ ] genetic algorithm? (Honorat)
- [ ] check that fig(6) is correct i.e. that always 8 sets of batteries. It is not! Compare with output of original script (stepsizes 2 and 5).
- [ ] determine search range from integration (needed since if no value is found at all for LLP=10 then script does not know which way to go)
- [ ] parallellize

## Todo list (tests)
- [ ] give error if SoC < 20% (It may never be zero)
- [x] give error if NPC_opt empty and define PV_opt, Batt_opt for that case 
- [x] check that range of PV and battery is divisible by PV_step and Batt_step
- [ ] evt. check that size of dimensions is larger than one s.t. taken over right dimension?
- [ ] give error if squeeze() squeezes more than time due to weird initial values

## Todo list (assumptions)
- [ ] coeff_cost_BoSeI is now 20% of total cost. Is this reliable or can it vary a lot? May have large impact on costs.

## Todo rewrite (would be nice; not needed)
- [ ] The loop 'for i = 1 : size(posPV, 1)' s.t. loop not needed
- [x] Combine all ifs to ifelse 
- [x] change batt_balance(1,t) to vector i.e. without 1 everywhere
