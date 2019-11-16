function [o] = GenerateRandomModel (recharge)

% addpath('../MonteCarlo/NaiveMCRuns/');

o = tdcFV('setdef');

Flds = fieldnames(o);
mpInd = GetMPInds(Flds);
[mbnds, logParams] = GetBounds(mpInd, recharge);

Vals = mbnds(:,1) + rand(length(mbnds),1).*diff(mbnds,[],2);

if recharge
    Vals(mpInd.Omega) = 10^Vals(mpInd.Omega);
end

o = FillFields(Flds, Vals, logParams);


end