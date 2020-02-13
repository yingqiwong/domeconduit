
clear all;
dbstop if error

AddPaths;
addpath('../CompareTDBE/');
BEFolder = '../../../../ModelTestFiles/TestBE/';
addpath(BEFolder);

%% load BE file info

load TestBE_manydp_Info.mat
LogParam = zeros(1,12);

VelErr = (Vel(:,1:end-1) - Vel(:,end))./Vel(:,end);
VolErr = (Vol(:,1:end-1) - Vol(:,end))./Vol(:,end);
DefErr = (Def(:,1:end-1) - Def(:,end))./Def(:,end);
CO2Err = (CO2(:,1:end-1) - CO2(:,end))./CO2(:,end);

%% choose a few files to evaluate

[~, tmpinds] = maxk(abs(VolErr),10);
inds = unique(tmpinds(:));
inds(inds==length(FileNames)) = [];

%% define some variables

Nf = length(FileNames);
NzVec = [201, 301, 401, 601, 801]

o = tdcFV('setdef');
mNames = fieldnames(o);

%% run test

NumWorkers = 10;
parpool(NumWorkers);

parfor (i = 1:length(inds), NumWorkers)
    fi = inds(i);
    
    fprintf('Running solution on file #%d...\n', fi);
    FileName = [BEFolder FileNames{fi}(1:end-4) '_dzbetd.mat'];
    RunModels(FileName, mNames, Model(fi,:), LogParam, NzVec, dpVec)
    
end
delete(gcp('nocreate'))
 
%% functions

function [] = RunModels (FileName, mNames, model, LogParam, NzVec, dpVec)

Nz = length(NzVec);
Np = length(dpVec);

td(Nz,1)  = struct('x', [], 'y', [], 'z', [], 'yp', [], 'xe', [], 'ye', [], 'ie', []);
be(Nz,Np) = struct('x', [], 'y', [], 'z', [], 'yp', [], 'xe', [], 'ye', [], 'ie', []);

tdRunTime = nan(Nz,1);
beRunTime = nan(Nz,Np);
PlugDepth = nan(Nz,1);
zphigc    = nan(Nz,1);

m = [];
mbe = [];

o = FillFields(mNames, model, LogParam);

for zi = 1:Nz
    fprintf('Nz for ss = %d.\n', NzVec(zi));
    
    try
        [tdtmp, m, tdRTtmp, betmp, mbe, beRTtmp] = RunTDBE(o, dpVec, NzVec(zi));
    catch
        fprintf('Something wrong with model with filename %s...\n',FileName);
        continue;
    end
    
    if isempty(tdtmp), continue; end
    
    td(zi)   = tdtmp;
    be(zi,:) = betmp;
    
    tdRunTime(zi)   = tdRTtmp;
    beRunTime(zi,:) = beRTtmp;
    PlugDepth(zi)   = m.PlugDepth;
    zphigc(zi)      = m.zphigc;
end

save(FileName, 'NzVec', 'dpVec', 'tdRunTime', 'td', 'be', ...
    'tdRunTime', 'beRunTime', 'm', 'mbe', 'PlugDepth', 'zphigc');

fprintf('Finished solution on file #%d.\nSaved in file %s.\n', fi, FileName);

end



function [td, m, tdRunTime, be, mbe, beRunTime] = RunTDBE (o, dpVec, Nz)

% initialize output variables
td = []; m = [];      tdRunTime = [];
be = struct('x', [], 'y', [], 'z', [], 'yp', [], 'xe', [], 'ye', [], 'ie', []);
beRunTime = nan(length(dpVec),1);
mbe = [];

% initialize ss options and allow structure alteration
opts = tdcFV('ss_init',o);
opts.Nz = Nz;

% run steady state model
[ss, opts, ssflag] = tdcFV('run_ssc_opts', opts);
if ssflag~=1, return; end

[m, y0, z] = tdcFV('td_init', ss.m, ss, 0, 0);
m.slv.dtFixed = 1; % do this to make the file sizes manageable.
m.Nt = 101;
m.Tspan = ConvertYearToSec([0,5]);

fprintf('Starting td solution...\n');
try
    tic;
    [td, m, flag] = tdcFV('run_tdc', y0, z, m);
    tdRunTime = toc;
catch
    td = [];
end
fprintf('Finished td solution.\n');


for idp = 1:length(dpVec)
    fprintf('Starting BE solution, dp = %.2f...\n', 1e-6*dpVec(idp));
    [mbe, y0, z] = tdcFV('td_init', ss.m, ss, 1, 0);
    mbe.be.dp = dpVec(idp);
    tic;
    [be(idp), ~] = RunBE(y0, z, mbe, 1);
    beRunTime(idp) = toc;
    fprintf('Finished BE solution with dp = %.2f.\n', 1e-6*dpVec(idp));
end

end