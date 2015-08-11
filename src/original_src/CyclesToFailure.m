function [ c_f ] = CyclesToFailure( DOD )

%  Compute the battery lifetime consumption with reference to the battery cycles exploited
c_f = 15790*exp(-11.96*DOD)+2633*exp(-1.699*DOD);

end
