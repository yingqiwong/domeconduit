
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
inds

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
    RunTDBE(FileName, mNames, Model(fi,:), LogParam, NzVec, dpVec)
    
end
delete(gcp('nocreate'))

