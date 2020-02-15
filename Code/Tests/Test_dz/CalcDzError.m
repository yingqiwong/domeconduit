function [] = CalcDzError (BEFolder)

addpath(BEFolder);

Files       = dir([BEFolder '*.mat']);
AllFileNames= {Files.name}';

% find only the model files
ModelFiles = AllFileNames(contains(AllFileNames, '_2019'));
Nf         = length(ModelFiles);

% output file with stats
OutFileName = [BEFolder 'TestBE_dzError.mat'];

% load one file to get sizes of Nz, dp
load(ModelFiles{1});
Nz = length(NzVec);

Model   = nan(Nf,12);
RT      = nan(Nz,Nf); 
Solved  = nan(Nz,Nf);
v0      = nan(Nz,Nf);
VelErr  = nan(Nz,Nf);
VolErr  = nan(Nz,Nf);
DefErr  = nan(Nz,Nf);
CO2Err  = nan(Nz,Nf);
VelEnd  = nan(Nz,Nf);
VolEnd  = nan(Nz,Nf);
DefEnd  = nan(Nz,Nf);
CO2Stt  = nan(Nz,Nf);

for fi = 1:Nf
    load([ModelFiles{fi}]);
    if isempty(m), continue; end
    
    o           = GetOFromM(m);
    Model(fi,:) = struct2array(o);
    RT(:,fi)    = tdRunTime;
    
    Vel = cell(Nz,1);
    Vol = cell(Nz,1);
    Def = cell(Nz,1);
    CO2 = cell(Nz,1);
    
    for zi = 1:Nz
        if ~isempty(td(zi).x)
            if td(zi).x(end) == m.Tspan(2)
                [v0(zi,fi), Vel{zi}, Vol{zi}, Def{zi}, CO2{zi}] = ...
                    CalcData(td(zi), m, PlugDepth(zi), zphigc(zi));
                
                Solved(zi,fi) = 1;
                VelEnd(zi,fi) = Vel{zi}(end);
                VolEnd(zi,fi) = Vol{zi}(end);
                DefEnd(zi,fi) = Def{zi}(end);
                CO2Stt(zi,fi) = CO2{zi}(1);
            end
        end
    end
    
    trueInd = nan;
    for zi = Nz:-1:1
        if Solved(zi,fi)==1
            if isnan(trueInd)
                VelTrue = Vel{zi};
                VolTrue = Vol{zi};
                DefTrue = Def{zi};
                CO2True = CO2{zi};
                trueInd = 1;
            else
                VelErr(zi,fi) = 100*1/VelTrue(1)*norm(Vel{zi}-VelTrue);
                VolErr(zi,fi) = 100*1/VolTrue(end)*norm(Vol{zi}-VolTrue);
                DefErr(zi,fi) = 100*1/DefTrue(end)*norm(Def{zi}-DefTrue);
                CO2Err(zi,fi) = 100*1/CO2True(1)*norm(CO2{zi}-CO2True);
            end
        else
            continue; 
        end
    end

end

save(OutFileName, ...
    'NzVec', 'Model', 'Solved', 'RT', 'v0', 'ModelFiles',...
    'VelErr', 'VolErr', 'DefErr', 'CO2Err', ...
    'VelEnd', 'VolEnd', 'DefEnd', 'CO2Stt')
end

function [v0, vEnd, Vol, Def, CO2] = CalcData (td, m, PlugDepth, zphigc)

m.Nz        = length(td.z);
m.dQUICKn   = tdcFV('QUICK', td.z');
m.PlugDepth = PlugDepth;
m.zphigc    = zphigc;

v0   = td.y(2,1)*m.slv.sy(2);
vEnd = td.y(end-3,:)*m.slv.sy(2);
Vol  = td.y(end,:);

[ur, ~, ~] = CalcJRO1Def(td, m, 0);
Def = 1e3*ur;
[~,CO2] = CalcGasEmissions(td,m,0);
end



