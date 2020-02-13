
clear all;
dbstop if error

AddPaths;
addpath('../CompareTDBE/');
addpath(genpath('../../../../DomeEruption_MSH2004/Code/Data/Deformation/'));
BEFolder = '../../../../ModelTestFiles/TestBE/';
addpath(BEFolder);

%%

load TestBE_manydp_20190909_05_dzbetd.mat

%%
NzTD   = 0.5*(NzVec-1)+1;
mvec   = repmat(m, length(td),1);
mbevec = repmat(mbe, length(td),1);
for zi = 1:length(td)
    mvec(zi).Nz      = NzTD(zi);
    mvec(zi).dQUICKn = tdcFV('QUICK', td(zi).z');
    mvec(zi).PlugDepth = PlugDepth(zi);
    phig0 = td(zi).y(3:4:end,1)*m.slv.sy(3);
    mvec(zi).zphigc = find(phig0>m.phi_gc,1);
    
    mbevec(zi).Nz      = NzTD(zi);
    mbevec(zi).dQUICKn = tdcFV('QUICK', be(zi,1).z');
    mbevec(zi).PlugDepth = PlugDepth(zi);
    phig0 = be(zi,1).y(3:4:end,1)*mbe.slv.sy(3);
    mbevec(zi).zphigc = find(phig0>mbe.phi_gc,1);
end


%% plot td vs be solutions

zi = 4;
CompareTDBE_Plotting('PlotCompareTDBE', td(zi), mvec(zi), be(zi,:), mbevec(zi), dpVec)

%% plot td solutions for different depth discretizations

PlotComparetdSols(td, mvec, num2str(NzTD'), 5);


















