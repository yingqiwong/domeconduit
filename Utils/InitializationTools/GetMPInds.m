function [mpInd] = GetMPInds (Flds)
% produces a structure mpInd which contains the indices of each field in
% the model parameter vector. 
% i.e. mpInd.V0 = 1; mpInd.Omega = 2; etc.
% useful for indexing the model parameter vector without knowing the order
% of parameters in the vector.

for ifld = 1:length(Flds)
    mpInd.(Flds{ifld}) = ifld;
end

end