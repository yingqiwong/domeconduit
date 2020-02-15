function [td, m, tdRunTime] = RunConvergence (o, NzVec, FileName)

% initialize output variables
td = struct('x', [], 'y', [], 'z', [], 'yp', [], 'xe', [], 'ye', [], 'ie', []);
m  = [];      

% initialize ss options and allow structure alteration
opts = tdcFV('ss_init',o);

tdRunTime = nan(length(NzVec),1);
PlugDepth = nan(length(NzVec),1);
zphigc    = nan(length(NzVec),1);

for zi = 1:length(NzVec)
    
    opts.Nz = NzVec(zi);
    
    fprintf('Nz for ss = %d.\n', opts.Nz);
    
    % run steady state model
    [ss, opts, ssflag] = tdcFV('run_ssc_opts', opts);
    if ssflag~=1, return; end
    
    [m, y0, z]    = tdcFV('td_init', ss.m, ss, 0, 0);
    m.slv.dtFixed = 1; % do this to make the file sizes manageable.
    m.Nt          = 101;
    m.Tspan       = ConvertYearToSec([0,5]);
    PlugDepth(zi) = m.PlugDepth;
    zphigc(zi)    = m.zphigc;
    
    fprintf('Starting td solution...\n');
    try
        tic;
        [tdout, m, flag] = tdcFV('run_tdc', y0, z, m);
        tdRunTime(zi) = toc;
        td(zi) = tdout;
        fprintf('Finished td solution.\n');
    catch
        fprintf('No td solution.\n');
        continue;
    end
    
end

if nargin == 3
    save(FileName, 'NzVec', 'tdRunTime', 'td', 'm', 'PlugDepth', 'zphigc');
    fprintf('Finished solution, saved in file %s.\n', FileName);
end



end