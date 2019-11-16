function [o] = FillFields (ParamNames, ParamVals, logParams)
% fills the structure of model parameters for tdcFV,
% o.ParamNames = ParamVals

if nargin == 2
    mpInd = GetMPInds(ParamNames);
    [~, logParams] = GetBounds(mpInd, ParamVals(strcmp(ParamNames, 'Omega')));
end

for ifld = 1:length(ParamNames)
    if (logParams(ifld))
        o.(ParamNames{ifld}) = 10.^ParamVals(ifld);
    else
        o.(ParamNames{ifld}) = ParamVals(ifld);
    end
end


end