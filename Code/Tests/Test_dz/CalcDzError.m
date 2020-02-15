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
RT      = nan(Nz,  Nf); 
v0      = nan(Nz-1,Nf);
VelErr  = nan(Nz-1,Nf);
VolErr  = nan(Nz-1,Nf);
DefErr  = nan(Nz-1,Nf);
CO2Err  = nan(Nz-1,Nf);

for fi = 1:Nf
    load([ModelFiles{fi}]);
    if isempty(m), continue; end
    
    o           = GetOFromM(m);
    Model(fi,:) = struct2array(o);
    RT(:,fi)    = tdRunTime;
    
    % calc data for max discretization model (assume to be true)
    if ~isempty(td(end).x)
        [v0(end,fi), VelTrue, VolTrue, DefTrue, CO2True] = CalcData(td(end), m, PlugDepth(end), zphigc(end));
    end
    
    for zi = 1:(Nz-1)
        if ~isempty(td(zi).x)
            if td(zi).x(end) == m.Tspan(2)
                [v0(zi,fi), Vel, Vol, Def, CO2] = ...
                    CalcData(td(zi), m, PlugDepth(zi), zphigc(zi));
                
                VelErr(zi,fi) = 100*1/VelTrue(1)*norm(Vel-VelTrue);
                VolErr(zi,fi) = 100*1/VolTrue(end)*norm(Vol-VolTrue);
                DefErr(zi,fi) = 100*1/DefTrue(end)*norm(Def-DefTrue);
                CO2Err(zi,fi) = 100*1/CO2True(1)*norm(CO2-CO2True);
            end
        end
        
    end

end

save(OutFileName, 'NzVec', 'Model', 'RT', 'v0', 'VelErr',...
    'VolErr', 'DefErr', 'CO2Err', 'ModelFiles')
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



