
function [mbnds, LogParams] = GetBounds (mpInd, recharge)


% bounds (all units in SI unless otherwise stated)
mbnds(mpInd.AR,:)         = [0.5,1];
mbnds(mpInd.L,:)          = [2,8]*1e3;
mbnds(mpInd.V0,:)         = [10,200]*1e9;
mbnds(mpInd.Omega,:)      = [0,0];
mbnds(mpInd.op,:)         = [10,30]*1e6;
mbnds(mpInd.total_h2o,:)  = [3,7];        % wt%
mbnds(mpInd.total_co2,:)  = [100,8000];   % ppm
mbnds(mpInd.R,:)          = [30,100];
mbnds(mpInd.f0,:)         = [-2,log10(0.6)]; % in log scale
mbnds(mpInd.a,:)          = [-3,-1];     % in log scale
mbnds(mpInd.k_lat,:)      = [-20,-10];       % in log scale
mbnds(mpInd.phi_gc,:)     = [0.2,0.4];

if (recharge)
    mbnds(mpInd.Omega,:)  = [0,6] + log10(ConvertOmega(1, 'm3/day/MPa'));       % in log space (1/1e6/24/3600 ~ 1e-11)
end

LogParams = zeros(size(mbnds,1),1);
LogParams(mpInd.f0)    = 1;
LogParams(mpInd.a)     = 1;
LogParams(mpInd.k_lat) = 1;


end