
AddPaths;
set(0,'defaultlinelinewidth', 2, 'defaultaxesfontsize', 16);

%%
shallow = csvread('beta_ch_shallow.txt');
medium  = csvread('beta_ch_medium.txt');
deep    = csvread('beta_ch_deep.txt');

%%

ARq = logspace(-1.2,1.2,51);
colors = parula(4);

figure;
h(1) = loglog(shallow(:,1), shallow(:,2), '+-', 'color', colors(3,:));
hold on;
loglog(ARq, pchip(shallow(:,1), shallow(:,2), ARq), 'color', colors(3,:))
h(2) = loglog(medium(:,1), medium(:,2), '+-', 'color', colors(2,:));
loglog(ARq, pchip(medium(:,1), medium(:,2), ARq), 'color', colors(2,:));
h(3) = loglog(deep(:,1), deep(:,2), '+-', 'color', colors(1,:));
loglog(ARq, pchip(deep(:,1), deep(:,2), ARq), 'color', colors(1,:));
plot(ARq, (1 + 1/3*log10(ARq)), 'r-');
% beta_ch_tmp = 1/3*4*20e9*Get_beta_ch(ARq, 10e3, 4e3, 20e9);
% plot(1./ARq, beta_ch_tmp, 'y--');
hold off;

leg = legend(h, {'shallow (0.75)', 'medium', 'deep'}, 'box', 'off');
title(leg, 'depth/
xlim([0.04,20]); ylim([0.8,100]);
xlabel('Aspect ratio (height/width)');
ylabel('Compressibility [\beta_{ch}/(3/4\mu)]');
