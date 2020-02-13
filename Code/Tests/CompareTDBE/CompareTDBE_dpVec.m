function varargout = CompareTDBE_dpVec (varargin)
[varargout{1:nargout}] = feval(varargin{:});
end

function [td, m, tdRunTime, be, mbe, beRunTime] = main (o, dpVec, varargin)

% initialize output variables
td = []; m = [];      tdRunTime = [];
be = struct('x', [], 'y', [], 'z', [], 'yp', [], 'xe', [], 'ye', [], 'ie', []);
beRunTime = nan(length(dpVec),1);
mbe = [];

% initialize ss options and allow structure alteration
opts = tdcFV('ss_init',o);
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opts.(args{1,ia}) = args{2,ia};
end

% run steady state model
[ss, opts, ssflag] = tdcFV('run_ssc_opts', opts);
if ssflag~=1, return; end

[m, y0, z] = tdcFV('td_init', ss.m, ss, 0, 0);
% m.PlugDepth = 0;

fprintf('Starting td solution...\n');
try
    tic;
    [td, m, flag] = tdcFV('run_tdc', y0, z, m);
    tdRunTime = toc;
catch
    td = [];
end
fprintf('Finished td solution.\n');

% if flag ~= 1, return; end

for idp = 1:length(dpVec)
    fprintf('Starting BE solution, dp = %.2f...\n', 1e-6*dpVec(idp));
    [mbe, y0, z] = tdcFV('td_init', ss.m, ss, 1, 0);
    mbe.be.dp = dpVec(idp);
%     mbe.PlugDepth = 0;
    tic;
    [be(idp), flag] = RunBE(y0, z, mbe, 1);
    beRunTime(idp) = toc;
    fprintf('Finished BE solution with dp = %.2f.\n', 1e-6*dpVec(idp));
end

end
