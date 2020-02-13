
function [Err] = PlotConvergence (td, m, NzVec)

Err = CalcConvergence(td, m, NzVec);
Vars = {'pch', 'vExit', 'phigExit', 'Vol', 'Def', 'CO2'};

figure; 
[Nrow, Ncol] = GetSubplotRowCol(size(Err,2));
for pi = 1:size(Err,2)
    subplot(Nrow, Ncol, pi);
    loglog(NzVec, Err(:,pi), '+-');
    title(Vars{pi});
end


end

function Err = CalcConvergence (td, m, NzVec)

Nsols = length(td);

Err  = nan(Nsols, 6);
% pch, vExit, phigExit, Vol, Def, CO2 at end time

[~, ind] = max(NzVec);
tdTrue = td(ind);
mTrue  = m(ind);
tdvTrue = extract_y(tdTrue, mTrue);
[VolTrue, DefTrue, CO2True] = CompareTDBE_Plotting('CalcData', tdTrue, mTrue);

for si = 1:Nsols
    if ind == si, continue; end
    
    tdv = extract_y(td(si),m(si));
    Err(si,1) = (tdv.p(1,end) - tdvTrue.p(1,end))/tdvTrue.p(1,end);
    Err(si,2) = (tdv.v(end,end) - tdvTrue.v(end,end))/tdvTrue.v(end,end);
    Err(si,3) = (tdv.phi_g(end,end) - tdvTrue.phi_g(end,end))/tdvTrue.phi_g(end,end);
    
    [vol, def, co2] = CompareTDBE_Plotting('CalcData', td(si), m(si));
    Err(si,4) = (vol(end) - VolTrue(end))/VolTrue(end);
    Err(si,5) = (def(end) - DefTrue(end))/DefTrue(end);
    Err(si,6) = (co2(1) - CO2True(1))/CO2True(1);
end
Err = 100*abs(Err);
end







