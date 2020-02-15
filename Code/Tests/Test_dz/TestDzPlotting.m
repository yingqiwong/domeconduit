
% run dz tests on all the models in the TestBE set
AddPaths;
DzFolder = '../../../../ModelTestFiles/Testdz/';
addpath('../CompareTDBE/');
addpath(DzFolder);

%%
clear all
DzFolder = '../../../../ModelTestFiles/Testdz/';
CalcDzError(DzFolder)

%%
clear all;
load ../../../../ModelTestFiles/Testdz/TestBE_dzError.mat
NzTD = 0.5*(NzVec-1)+1;
Nf = length(ModelFiles);

%% plot model parameters for different td dz disc

o = tdcFV('setdef');
mNames = fieldnames(o);
LogParam = zeros(length(mNames),1);
LogParam(10:12) = 1;
[Nrow, Ncol] = GetSubplotRowCol(length(mNames));

tdSolved = ~isnan(RT);

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

%% final vol, def error

VelDzErr = 100*abs((VelEnd(1:end-1,:)-VelEnd(end,:))./VelEnd(end,:));
VolDzErr = 100*abs((VolEnd(1:end-1,:)-VolEnd(end,:))./VolEnd(end,:));
DefDzErr = 100*abs((DefEnd(1:end-1,:)-DefEnd(end,:))./DefEnd(end,:));
CO2DzErr = 100*abs((CO2Stt(1:end-1,:)-CO2Stt(end,:))./CO2Stt(end,:));

MaxVol = log10(max(max(VolEnd)));
MinVol = log10(min(min(VolEnd)));
colors = 1/(MaxVol-MinVol)*(log10(nanmean(VolEnd))'-MinVol).*[1,1,1];

figure;
set(gcf,'Position', [400,400,1000,350]);
for fi = 1:Nf
    subplot(121);
    semilogy(NzTD(1:end-1), VolDzErr(:,fi), 'ko-', 'MarkerSize', 10, ...
        'MarkerFaceColor', colors(fi,:)); 
    hold on;
    
    subplot(122);
    semilogy(NzTD(1:end-1), DefDzErr(:,fi), 'ko-', 'MarkerSize', 10, ...
        'MarkerFaceColor', colors(fi,:));  
    hold on;
    
end

subplot(121); grid on;
xlabel('Nz'); ylabel('Final volume err (%)');
title('Extruded volume error'); 

subplot(122); grid on;
xlabel('Nz'); ylabel('Final JRO1 def err (%)');
title('Deformation error');


%% time series error
figure;
set(gcf,'Position', [400,400,1200,350]);
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

%% runtimes

[Nrow, Ncol] = GetSubplotRowCol(length(NzVec));

figure;
for zi = 1:length(NzVec)
    subplot(Nrow, Ncol, zi);
    histogram(RT(zi,:), 10);
    hold on;
    meanval = nanmean(RT(zi,:));
    plot(meanval*ones(1,2), ylim, 'r-', 'linewidth', 2);
    ylimits = ylim;
    text(meanval+5, ylimits(2), num2str(meanval,3), 'FontSize', 14,...
        'color', 'r', 'VerticalAlignment','top');
    title(['Nz = ' num2str(NzVec(zi))]);
end


%% plot one file

for ind = 51:Nf
    ModelFiles{ind}
load(ModelFiles{ind})

NzTD   = 0.5*(NzVec-1)+1;
mvec   = repmat(m, length(td),1);
SolInd = nan(length(td),1);
for zi = 1:length(td)
    if isempty(td(zi).z), continue; end
    
    mvec(zi).Nz      = NzTD(zi);
    mvec(zi).dQUICKn = tdcFV('QUICK', td(zi).z');
    mvec(zi).PlugDepth = PlugDepth(zi);
    mvec(zi).zphigc = zphigc(zi);
   
    SolInd(zi) = true;
end

PlotDzSols(td(SolInd==1), mvec(SolInd==1)); drawnow
% Err = PlotConvergence(td(SolInd==1), mvec(SolInd==1), NzTD(SolInd==1));
end

