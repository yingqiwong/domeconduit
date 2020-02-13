function [] = ExtractTestDzTDBEStats (BEFolder)

addpath(BEFolder);

Files       = dir([BEFolder '*.mat']);
AllFileNames= {Files.name}';

% find only the model files
ModelFiles = AllFileNames(contains(AllFileNames, '_2019'));
Nf         = length(ModelFiles);

% output file with stats
OutFileName = [BEFolder 'TestBE_dzbetd_Stats.mat'];

% load one file to get sizes of Nz, dp
load(ModelFiles{1});
Nz = length(NzVec);
Np = length(dpVec);

Model = nan(Nf,12);

RT   = nan(Nz,Np+1,Nf); % first column is td solution
beNt = nan(Nz,Np  ,Nf);
v0   = nan(Nz,Np+1,Nf);
vEnd = nan(Nz,Np+1,Nf);
Vol  = nan(Nz,Np+1,Nf);
Def  = nan(Nz,Np+1,Nf);
CO2  = nan(Nz,Np+1,Nf);

for fi = 1:Nf
    load([ModelFiles{fi}]);
    if isempty(m), continue; end
    
    o           = GetOFromM(m);
    Model(fi,:) = struct2array(o);
    
    RT(:,:,fi) = [tdRunTime, beRunTime];
    
    for zi = 1:Nz
        if ~isempty(td(zi).x)
            if td(zi).x(end) == m.Tspan(2)
                [v0(zi,1,fi), vEnd(zi,1,fi), Vol(zi,1,fi), Def(zi,1,fi), CO2(zi,1,fi)] = ...
                    CalcData(td(zi), m, PlugDepth(zi), zphigc(zi));
            end
        end
        
        for pi = 1:Np
            if ~isempty(be(zi,pi).x)
                if be(zi,pi).x(end) == m.Tspan(2)
                    [v0(zi,1+pi,fi), vEnd(zi,1+pi,fi), Vol(zi,1+pi,fi), Def(zi,1+pi,fi), CO2(zi,1+pi,fi)] = ...
                        CalcData(be(zi,pi), mbe, PlugDepth(zi), zphigc(zi));
                    beNt(zi,pi,fi) = length(be(zi,pi).x);
                end
            end
        end
    end
    
end

save(OutFileName, 'NzVec', 'dpVec', 'Model', 'RT', 'beNt', 'v0', 'vEnd',...
    'Vol', 'Def', 'CO2', 'ModelFiles')
end

function [v0, vEnd, Vol, Def, CO2] = CalcData (td, m, PlugDepth, zphigc)

m.Nz        = length(td.z);
m.dQUICKn   = tdcFV('QUICK', td.z');
m.PlugDepth = PlugDepth;
m.zphigc    = zphigc;

v0   = td.y(2,1)*m.slv.sy(2);
vEnd = td.y(end-3,end)*m.slv.sy(2);
Vol  = td.y(end,end);

[ur, ~, ~] = CalcJRO1Def(td, m, 0);
Def = 1e3*ur(end);
[~,CO2tmp] = CalcGasEmissions(td,m,0);
CO2 = CO2tmp(1);
end

