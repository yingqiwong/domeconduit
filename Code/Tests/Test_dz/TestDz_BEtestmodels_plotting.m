
clear all;
dbstop if error

AddPaths;
addpath('../CompareTDBE/');
addpath(genpath('../../../../DomeEruption_MSH2004/Code/Data/Deformation/'));
BEFolder = '../../../../ModelTestFiles/TestBE/';
addpath(BEFolder);
set(0, 'defaultaxesfontsize', 18);

%%

Files       = dir([BEFolder '*.mat']);
AllFileNames= {Files.name}';

% find only the model files
ModelFiles = AllFileNames(contains(AllFileNames, '_2019'));
Nf         = length(ModelFiles);

%% plot one file

% ind = randi(Nf,1);
ind = 3;
ModelFiles{ind}
load(ModelFiles{ind})

NzTD   = 0.5*(NzVec-1)+1;
mvec   = repmat(m, length(td),1);
mbevec = repmat(mbe, length(td),1);
SolInd = nan(length(td),1);
for zi = 1:length(td)
    if isempty(td(zi).z), continue; end
    
    mvec(zi).Nz      = NzTD(zi);
    mvec(zi).dQUICKn = tdcFV('QUICK', td(zi).z');
    mvec(zi).PlugDepth = PlugDepth(zi);
    mvec(zi).zphigc = zphigc(zi);
    
    mbevec(zi).Nz      = NzTD(zi);
    mbevec(zi).dQUICKn = tdcFV('QUICK', be(zi,1).z');
    mbevec(zi).PlugDepth = PlugDepth(zi);
    mbevec(zi).zphigc = zphigc(zi);
    
    SolInd(zi) = true;
end

zi = 2;
CompareTDBE_Plotting('PlotCompareTDBE', td(zi), mvec(zi), be(zi,:), mbevec(zi), dpVec)

PlotDzSols(td(SolInd==1), mvec(SolInd==1));
Err = PlotConvergence(td(SolInd==1), mvec(SolInd==1), NzTD(SolInd==1));

%%

figure;
colors = parula(length(dpVec)+1);
for idp = 1:length(dpVec)
    h(idp) = semilogy(NzTD, beRunTime(:,idp)./tdRunTime, 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', colors(idp,:), 'MarkerEdgeColor', colors(idp,:));
    hold on;
end
plot(xlim, [1,1], 'k-', 'linewidth', 2);
lg = legend(h, num2str(1e-6*dpVec'), 'Location', 'northeast');
title(lg, 'dp (MPa)');
hold off;
grid on;
xlabel('Nz'); ylabel('be runtime/td runtime');

%%-------------------------------------------------------------------------
%% calc stats for all the files

BEFolder = '../../../../ModelTestFiles/TestBE/';
ExtractTestDzTDBEStats(BEFolder);

%% look at all the files together

clear all
load TestBE_dzbetd_Stats.mat

Nf = length(ModelFiles);
NzTD   = 0.5*(NzVec-1)+1;
Np = length(dpVec);
Nz = length(NzVec);

VelErr = 100*abs((vEnd(:,2:end,:) - vEnd(:,1,:))./vEnd(:,1,:));
VolErr = 100*abs((Vol(:,2:end,:) - Vol(:,1,:))./Vol(:,1,:));
DefErr = 100*abs((Def(:,2:end,:) - Def(:,1,:))./Def(:,1,:));
CO2Err = 100*abs((CO2(:,2:end,:) - CO2(:,1,:))./CO2(:,1,:));

VelDzErr = squeeze(100*abs((vEnd(1:end-1,1,:)-vEnd(end,1,:))./vEnd(end,1,:)));
VolDzErr = squeeze(100*abs((Vol(1:end-1,1,:)-Vol(end,1,:))./Vol(end,1,:)));
DefDzErr = squeeze(100*abs((Def(1:end-1,1,:)-Def(end,1,:))./Def(end,1,:)));
CO2DzErr = squeeze(100*abs((CO2(1:end-1,1,:)-CO2(end,1,:))./CO2(end,1,:)));

%% plot runtimes

figure;
for idp = 1:4
    subplot(2,2,idp);
    tdRT = squeeze(RT(:,1,:));
    beRT = squeeze(RT(:,1+idp,:));
    
    loglog(NzVec, beRT./tdRT, 'o-', 'MarkerSize', 10);
    hold on;
    plot(xlim, [1,1], 'k-', 'LineWidth', 2);
    hold off;
    grid on;
    xlabel('Nz'); ylabel('be/td runtime');
    title(['dp = ' num2str(dpVec(idp)*1e-6) ' MPa']);
end

%% runtimes different view
figure;
set(gcf,'Position', [300,300,1200,600]);
colors = parula(Np);
for zi = 1:length(NzTD)
    subplot(2,3,zi);
    tdRT = 1/60*squeeze(RT(zi,1,:))';
    beRT = 1/60*squeeze(RT(zi,2:end,:));
    
    for fi = 1:Nf
        plot(tdRT(fi)*ones(Nz,1), beRT(:,fi), 'ko-', 'MarkerSize', 10); hold on;
        for idp = 1:Np
            h(idp) = plot(tdRT(fi), beRT(idp,fi), 'o', 'MarkerSize', 10, ...
                'MarkerFaceColor', colors(idp,:), 'MarkerEdgeColor', 'k');
        end
    end
    axis tight
    xlabel('td runtime (min)'); ylabel('be runtime(min)');
    title(['Nz = ' num2str(NzTD(zi)), ', Nmodels = ' num2str(sum(~isnan(tdRT)))]);
    set(gca,'YScale','log','XScale','log');
    plot(xlim,xlim,'r-');

    
    if zi == 5
        hl = legend(h, num2str(1e-6*dpVec'));
        title(hl, 'dp');
    end
end


%% plot extruded volume error as function of number of time points
figure;
colors = parula(length(dpVec));
set(gcf,'Position',[514   422   1100   386]);

zi = 3;
subplot(121);
for idp = 1:length(dpVec)
    loglog(squeeze(beNt(zi,idp,:)), squeeze(VolErr(zi,idp,:)), 'o', 'MarkerSize', 10,...
        'MarkerFaceColor', colors(idp,:), 'MarkerEdgeColor', 'k');
    hold on;
end
xlim([1,200]);ylim([0.05,10]);
plot(xlim, fliplr(xlim)*0.1, 'k-');
hold off;
hl = legend(num2str(1e-6*dpVec'), 'Location', 'southwest');
title(hl, 'dp (MPa)');
xlabel('Number of time points'); ylabel('Volume err %'); grid on;

subplot(122);
for idp = 1:length(dpVec)
    scatter(squeeze(beNt(zi,idp,:)), squeeze(VolErr(zi,idp,:)), 80, ...
        squeeze(log10(v0(zi,idp,:))), 'filled', 'MarkerEdgeColor', 'k');
    hold on;
end
xlim([1,200]); ylim([0.05,10]);
plot(xlim, fliplr(xlim)*0.1, 'k-');
set(gca,'XScale','log','YScale','log','Box', 'on');
hold off;
cb = colorbar; title(cb, 'log10(v0)'); caxis([-4,-1]);
xlabel('Number of time points'); ylabel('Volume err %'); grid on;

%% plot errors for different dp

zi = 3;

figure;
set(gcf,'Position',[514   422   1000   350]);
dpScat = 1e-6*repmat(dpVec,Nf,1);
VolErrzi = squeeze(VolErr(zi,:,:))';
DefErrzi = squeeze(DefErr(zi,:,:))';
Volzitd = repmat(squeeze(Vol(zi,1,:)),length(dpVec),1);

subplot(121);
scatter(dpScat(:), VolErrzi(:), 80, Volzitd, 'filled', 'MarkerEdgeColor','k');
cb=colorbar; grid on; caxis([0,300]); xlim([0,2.5]);
set(gca,'yscale','log'); ylim([0.05,10]);
xlabel('dp (MPa)'); ylabel('vol err (%)');
title(cb,'vol at t=4yr');
title('Extruded volume error');

subplot(122);
scatter(dpScat(:), DefErrzi(:), 80, Volzitd, 'filled', 'MarkerEdgeColor','k');
cb=colorbar; grid on; caxis([0,300]); xlim([0,2.5]);
set(gca,'yscale','log'); ylim([0.05,10]);
xlabel('dp (MPa)'); ylabel('def err (%)');
title(cb,'vol at t=4yr');
title('JRO1 radial def error');

%%
figure;
set(gcf,'Position', [400,400,1000,350]);
for fi = 1:Nf
    subplot(121);
    plot(NzTD(1:end-1), VolDzErr(:,fi), 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); 
    hold on;
    
    subplot(122);
    plot(NzTD(1:end-1), DefDzErr(:,fi), 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); 
    hold on;
    
end

subplot(121);
xlabel('Nz'); ylabel('Final volume err (%)');
title('Extruded volume error');

subplot(122);
xlabel('Nz'); ylabel('Final JRO1 def err (%)');
title('Deformation error');

%% plot model parameters for different td dz disc

o = tdcFV('setdef');
mNames = fieldnames(o);
LogParam = zeros(length(mNames),1);
LogParam(10:12) = 1;
[Nrow, Ncol] = GetSubplotRowCol(length(mNames));

tdRT =  1/60*squeeze(RT(:,1,:));
tdSolved = ~isnan(tdRT);

figure;
hAx = tight_subplot(Nrow, Ncol, [0.1,0.05]);
for mi = 1:length(mNames)
    axes(hAx(mi));
    for fi = 1:Nf
        plot(Model(fi,mi)*ones(1,5), NzTD, 'k-'); hold on;
        plot(Model(fi,mi)*ones(1,sum(tdSolved(:,fi))), NzTD(tdSolved(:,fi)), ...
            'ko', 'markerfacecolor', 'b');
        plot(Model(fi,mi)*ones(1,5-sum(tdSolved(:,fi))), NzTD(~tdSolved(:,fi)),...
            'rx', 'LineWidth', 2);
    end
    if (LogParam(mi)), set(gca,'XScale','log'); end
    title(mNames{mi});
    ylim([100,410]);
end


%% load stats for the dz error

clear all
load TestBE_dzError.mat

%%
figure;
set(gcf,'Position', [400,400,1200,350]);
NzTD   = 0.5*(NzVec-1)+1;
subplot(131); 
semilogy(NzTD(1:end-1), VolErr, 'ko-', 'MarkerFaceColor', 'b', 'MarkerSize',8); 
grid on;
title('Volume time series error'); xlim([100,302]);
xlabel('Nz'); ylabel('% error in time series');

subplot(132); 
semilogy(NzTD(1:end-1), -DefErr, 'ko-', 'MarkerFaceColor', 'b', 'MarkerSize',8); 
grid on;
xlabel('Nz'); title('JRO1 def time series error'); xlim([100,302]);

subplot(133); 
semilogy(NzTD(1:end-1), CO2Err, 'ko-', 'MarkerFaceColor', 'b', 'MarkerSize',8); 
grid on;
xlabel('Nz'); title('CO2 time series error'); xlim([100,302]);





