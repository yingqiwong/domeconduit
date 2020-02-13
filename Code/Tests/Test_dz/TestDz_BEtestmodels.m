
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

[~, tmpinds] = maxk(VolErr,20);
inds = unique(tmpinds(:));

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

o = FillFields(mNames, model, LogParam);

for zi = 1:Nz
    fprintf('Nz for ss = %d.\n', NzVec(zi));
    
    try
        [tdtmp, m, tdRTtmp, betmp, mbe, beRTtmp] = CompareTDBE_dpVec('main', o, dpVec, 'Nz', NzVec(zi));
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
end

save(FileName, 'NzVec', 'dpVec', 'tdRunTime', 'td', 'be', ...
    'tdRunTime', 'beRunTime', 'm', 'mbe', 'PlugDepth');

fprintf('Finished solution on file #%d.\nSaved in file %s.\n', fi, FileName);

end