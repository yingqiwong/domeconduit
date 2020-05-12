
AddPaths

%%
clear all
load pchQ_NominalModel_20190504.mat

%%

inds = [1,3,5];
pchInd = 3;
tMax = [0.05,0.5,15];
VolPower = log10(Vols)-9;

figure;
hAx = tight_subplot(1,3,0.03,0.2,0.1);
set(gcf,'Position',[360   469   813   320]);
for i = 1:length(inds)
    mi = inds(i);
    t = squeeze(tvec(pchInd,mi,:));
    p = squeeze(pchtd(pchInd,mi,:));
    q = squeeze(qtd(pchInd,mi,:));
    
    ts = linspace(0,tMax(i),21);
    ps = interp1(t,p,ts);
    qs = interp1(pchss, qss, ps, 'pchip');
    
    axes(hAx(i));
    semilogy(t,q,ts,qs,'+','MarkerSize',8);
    xlim([0,tMax(i)]);
    xlabel('Time (years)'); 
    set(gca,'Box','off');
end

axes(hAx(1)); ylabel('Magma flux (m^3/s)');
