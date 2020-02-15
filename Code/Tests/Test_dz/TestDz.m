
% run dz tests on all the models in the TestBE set

clear all;
dbstop if error

AddPaths;
addpath('../CompareTDBE/');
BEFolder = '../../../../ModelTestFiles/TestBE/';
DzFolder = '../../../../ModelTestFiles/Testdz/';
addpath(BEFolder);

%% load BE file info

load TestBE_manydp_Info.mat
LogParam = zeros(1,12);

%% define some variables

Nf = length(FileNames);
NzVec = [201, 301, 401, 601, 801]

o = tdcFV('setdef');
mNames = fieldnames(o);

%% run test

NumWorkers = 10;
parpool(NumWorkers);

parfor (fi = 1:Nf, NumWorkers)
    
    fprintf('Running solution on file #%d...\n', fi);
    DzFileName = [DzFolder FileNames{fi}(1:end-4) '_dz_.mat'];
    
    o = FillFields(mNames, Model(fi,:), LogParam);
    RunConvergence(o, NzVec, DzFileName);
    
end
delete(gcp('nocreate'))

