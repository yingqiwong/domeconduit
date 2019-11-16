function [ur, ut, uz] = CalcJRO1Def (td, m, plot_opt)

% calculate radial displacements 10 km from the ellipsoid center ==
% approx distance of JRO1 from deformation center in Lisowski et al. (2008)
% Output in m

nu      = 0.25;         % Poisson's ratio
plunge  = 89.99;        % plunge angle of ellipsoid, deg (90 = prolate)
strike  = 0;            % strike of ellipsoid, deg (0 = aligned north)

dp = m.slv.sy(1)*(td.y(1,:) - td.y(1,1));

params = [0;0;m.ch.depth;0;m.ch.a;m.ch.a*m.ch.AR;deg2rad(plunge);deg2rad(strike)];
lambda = 2*m.ch.mu*nu/(1-2*nu);
matrl = [lambda;m.ch.mu;nu];

% initialize
Nt = length(td.x);
ur = zeros(Nt,1);     % radial
ut = zeros(Nt,1);     % tangential
uz = zeros(Nt,1);     % vertical

% calculate displacements
for ti = 1:Nt
    params(4) = dp(ti);
    [ur(ti), ut(ti), uz(ti)] = fcn_yangM(params, 8.73785e3, 0, matrl, 0);
end

if plot_opt == 1
    
    load CleanedData_Oct18.mat
    
    % load the eruption start date
    StartDate = datenum(datestr(EruptionStartDate,'dd-mmm-yy',1900));
    StartYear = datenum(datestr(datetime(2004,1,1),'dd-mmm-yy',1900));
    StartDateFrac = 2004 + 1/365*(StartDate - StartYear);
    
    % propagate error from east and north into radial displacements
    xy = llh2local([MSH_GPS.stations(4).lon; MSH_GPS.stations(4).lat; 0], [-122.1888;46.1991;0]);
    [th,~] = cart2pol(xy(1), xy(2));
    RotMat = [cos(th), sin(th); -sin(th), cos(th)];
    for ti = 1:length(MSH_GPS.stations(4).Date)
        tmp = RotMat*diag(MSH_GPS.stations(4).data(ti,4:5),0)*RotMat';
        RadErr(ti) = tmp(1);
    end
    
    RadErr(RadErr > (mean(RadErr)+std(RadErr))) = nan;
    
    figure;
    h1 = errorbar(MSH_GPS.stations(4).Date(1:end) - StartDateFrac, ...
         MSH_GPS.stations(4).rad_clean(1:end) + 16, RadErr(1:end), ...
         '.','color', 0.8*ones(3,1),'linewidth',1,...
         'markersize',5,'markerfacecolor', 'k', 'markeredgecolor', 'k', 'CapSize', 0);
    hold on; 
    h2 = plot(ConvertSecToYear(td.x), 1e3*ur, 'b-','linewidth',2);
    plot([0,0], ylim, 'k-');
    xlim([-1,4]); hold off;
    title('JRO1');
    xlabel('Time (years)'); ylabel('Radial displacements (mm)');
    legend([h1 h2], {'Data', 'Model'}); legend boxoff;
end

end