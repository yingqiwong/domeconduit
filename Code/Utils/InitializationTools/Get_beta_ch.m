function [beta_ch] = Get_beta_ch (AR, depth, a, mu)
% [beta_ch] = Get_beta_ch (AR, depth, a, mu)
% be careful - AR in beta_ch_txtfiles is Height/width. 
% AR in tdcFV is width/height

DepthBreak = mean([0.75,3.7]);
Rs = a*(AR)^(2/3);

if depth/Rs<DepthBreak
    ARcurve = csvread('beta_ch_shallow.txt');
else
    ARcurve = csvread('beta_ch_medium.txt'); 
end

beta_ch = 3/4/mu*pchip(ARcurve(:,1), ARcurve(:,2), 1./AR);

end