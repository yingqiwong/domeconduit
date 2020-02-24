
AddPaths;

%%
shallow = csvread('beta_ch_shallow.txt');
medium  = csvread('beta_ch_medium.txt');
deep    = csvread('beta_ch_deep.txt');

%%

ARq = logspace(-1,1.2,51);
colors = lines(3);

figure;
loglog(shallow(:,1), shallow(:,2), '+', 'color', colors(1,:))
hold on;
loglog(ARq, pchip(shallow(:,1), shallow(:,2), ARq), 'color', colors(1,:))
loglog(medium(:,1), medium(:,2), '+', 'color', colors(2,:)); hold on;
loglog(ARq, pchip(medium(:,1), medium(:,2), ARq), 'color', colors(2,:));
loglog(deep(:,1), deep(:,2), '+', 'color', colors(3,:));
loglog(ARq, pchip(deep(:,1), deep(:,2), ARq), 'color', colors(3,:));
axis manual
plot(ARq, (1 + 1/3*log10(ARq)));
beta_ch_tmp = 1/3*4*20e9*Get_beta_ch(ARq, 10e3, 4e3, 20e9);
plot(1./ARq, beta_ch_tmp, 'k--');
hold off;
