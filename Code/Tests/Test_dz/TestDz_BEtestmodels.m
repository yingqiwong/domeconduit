
clear all;
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

[~, tmpinds] = maxk(VolErr,10);
inds = unique(tmpinds(:));

%% define some variables

Nf = length(FileNames);
NzVec = [301, 401];
Nz = length(NzVec);

% dpVec2 = dpVec(2:3);
Np = length(dpVec);
o = tdcFV('setdef');
mNames = fieldnames(o);

%% run test

for fi = inds'
    
    fprintf('Running solution on file #%d...\n', fi);
    
    td(Nz)    = struct('x', [], 'y', [], 'z', [], 'yp', [], 'xe', [], 'ye', [], 'ie', []);
    be(Nz,Np) = struct('x', [], 'y', [], 'z', [], 'yp', [], 'xe', [], 'ye', [], 'ie', []);
    
    tdRunTime = zeros(Nz,1);
    beRunTime = zeros(Nz,Np);

    o = FillFields(mNames, Model(fi,:), LogParam);
    
    for zi = 1:Nz
        fprinf('Nz for ss = %d.\n', NzVec(zi));
        
        [tdtmp, m, tdRTtmp, betmp, mbe, beRTtmp] = CompareTDBE_dpVec('main', o, dpVec, 'Nz', NzVec(zi));
        
        if isempty(tdtmp), continue; end
        
        td(zi)   = tdtmp;
        be(zi,:) = betmp;
        
        tdRunTime(zi)   = tdRTtmp;
        beRunTime(zi,:) = beRTtmp;
    end
    
    NewFileName = [BEFolder FileNames{fi}(1:end-4) '_dzbetd.mat'];
    
    save(NewFileName, 'NzVec', 'dpVec', 'tdRunTime', 'td', 'be', ...
        'tdRunTime', 'beRunTime', 'm', 'mbe');
    
    fprintf('Finished solution on file #%d.\nSaved in file %s.\n', ...
        fi, NewFileName);

    clear td be tdRunTime beRunTime m mbe
end