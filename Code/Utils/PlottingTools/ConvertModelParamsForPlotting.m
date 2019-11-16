function [FldNames, FldVals] = ConvertModelParamsForPlotting (params, mvec, units)
% Converts the code version of the model parameter names and values to
% text version of field names and intuitive units
%
% Inputs
% params    EITHER o structure (containing variables as defined in tdcFV)
%           OR names of parameters cell array (Nparams x 1)
% mvec      matrix containing parameter values (Nparams x Niter)

if isstruct(params) % if input is structure o
    Flds = fieldnames(params);
else                % if input is parameter name cell array
    Flds = params;
end

if nargin == 1
    mvec = struct2array(params)';
    units = 1;
elseif nargin == 2 
    units = 1;
end

NFlds = length(Flds);

% check that size of mvec is Nparams x Niter
if size(mvec,1) ~= NFlds, mvec = mvec'; end

for ifld = 1:NFlds
    mpInd.(Flds{ifld}) = ifld;
end

if (units)
    FldNames{mpInd.V0}          = 'V_0 (km^3)';
    FldNames{mpInd.Omega}       = '\Omega (m^3/day/MPa)';
    FldNames{mpInd.AR}          = '\alpha';
    FldNames{mpInd.L}           = 'L (km)';
    FldNames{mpInd.op}          = '\Deltap_0 (MPa)';
    FldNames{mpInd.R}           = 'R (m)';
    FldNames{mpInd.total_h2o}   = '\chi_h^{ch} (wt%)';
    FldNames{mpInd.total_co2}   = '\chi_c^{ch} (ppm)';
    FldNames{mpInd.k_lat}       = 'k_c (m^2)';
    FldNames{mpInd.phi_gc}      = '\phi_{gc}';
    FldNames{mpInd.f0}          = 'f_0';
    FldNames{mpInd.a}           = 'a';
else
    FldNames{mpInd.V0}          = 'V_0';
    FldNames{mpInd.Omega}       = '\Omega';
    FldNames{mpInd.AR}          = '\alpha';
    FldNames{mpInd.L}           = 'L';
    FldNames{mpInd.op}          = '\Deltap_0';
    FldNames{mpInd.R}           = 'R';
    FldNames{mpInd.total_h2o}   = '\chi_h^{ch}';
    FldNames{mpInd.total_co2}   = '\chi_c^{ch}';
    FldNames{mpInd.k_lat}       = 'k_c';
    FldNames{mpInd.phi_gc}      = '\phi_{gc}';
    FldNames{mpInd.f0}          = 'f_0';
    FldNames{mpInd.a}           = 'a';
end

FldVals = mvec;
FldVals(mpInd.V0,:)     = 1e-9*mvec(mpInd.V0,:);
FldVals(mpInd.Omega,:)  = ConvertOmega(mvec(mpInd.Omega,:), 'm3/s/Pa');
FldVals(mpInd.L,:)      = 1e-3*mvec(mpInd.L,:);
FldVals(mpInd.op,:)     = 1e-6*mvec(mpInd.op,:);

end