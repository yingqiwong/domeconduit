clear all;

Scales = [1e5,1e6,1e7,1e8];

o = tdcFV('setdef');
opts = tdcFV('ss_init',o);
[ss, ssflag] = smf_rad_dz('solve', opts);
[m, y0, z] = tdcFV('td_init', ss.m, ss, 0, 0);

RunTime = zeros(length(Scales),1);

for si = 1:length(Scales)
    m.slv.syVol = Scales(si);
    
    tic;
    [td, m, flag] = tdcFV('run_tdc', y0, z, m);
    RunTime(si) = toc;
    
    sol(si).m  = m;
    sol(si).td = td;
end

%%

figure;
for si = 1:length(Scales)
    tyr = ConvertSecToYear(sol(si).td.x);
    plot(tyr, sol(si).td.y(end,:)*sol(si).m.slv.syVol*1e-6); hold on;
    pause
end
