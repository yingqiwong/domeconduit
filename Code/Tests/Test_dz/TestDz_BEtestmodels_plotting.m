
clear all;
dbstop if error

AddPaths;
addpath('../CompareTDBE/');
BEFolder = '../../../../ModelTestFiles/TestBE/';
addpath(BEFolder);

%%

load TestBE_manydp_20190906_02_dzbetd.mat

%%
NzTD   = 0.5*(NzVec-1)+1;
mvec   = repmat(m, length(td),1);
mbevec = repmat(mbe, length(td),1);
for zi = 1:length(td)
    mvec(zi).Nz      = NzTD(zi);
    mvec(zi).dQUICKn = tdcFV('QUICK', td(zi).z');
    
    mbevec(zi).Nz      = NzTD(zi);
    mbevec(zi).dQUICKn = tdcFV('QUICK', be(zi,1).z');
end

mvec(1).PlugDepth = 118;
mbevec(1).PlugDepth = 118;

    
CompareTDBE_Plotting('PlotCompareTDBE', td(1), mvec(1), be(1,:), mbevec(1), dpVec)