function [] = PlaytdSol (td, m, filename)

figHand = figure;
set(gcf, 'position', [520 250 1046 643], 'color', 'w');

% open movie file
vidObj = VideoWriter(filename);
vidObj.FrameRate = 10;
vidObj.Quality = 100;

open(vidObj);

tdvars = extract_y(td, m);
t = ConvertSecToYear(td.x);
z = -1e-3*td.z;

for ti = 1:length(td.x)
    
    subplot(341);
    plot(t, 1e-6*tdvars.p(end,:), 'k', t(ti), 1e-6*tdvars.p(end,ti), 'r*');
    title('Pressure (MPa)'); xlabel('Time (year)');
    
    subplot(342);
    plot(t, 1e3*tdvars.v(end,:), 'k', t(ti), 1e3*tdvars.v(end,ti), 'r*');
    title('Velocity (x10^{-3} m/s)'); xlabel('Time (year)');
    
    subplot(343);
    plot(t, 100*tdvars.phi_g(end,:), 'k', t(ti), 100*tdvars.phi_g(end,ti), 'r*');
    title('Porosity (%)'); xlabel('Time (year)');
    
    subplot(344);
    plot(t, tdvars.mw(end,:), 'k', t(ti), tdvars.mw(end,ti), 'r*');
    title('Mole fraction water'); xlabel('Time (year)');
    
    subplot(3,4,[5,9]);
    plot(1e-6*tdvars.p(:,1), z, 'k:', 1e-6*tdvars.p(:,ti), z, 'r-');
    ylabel('Depth (km)'); set(gca, 'ydir', 'reverse');
    
    subplot(3,4,[6,10]);
    plot(1e3*tdvars.v(:,1), z, 'k:', 1e3*tdvars.v(:,ti), z, 'r-');
    set(gca, 'ydir', 'reverse');
%     xlimits = xlim; xlim([0,xlimits(2)]); xlim manual
    xlim([0,3]);
    hold on;
    text(1e3*tdvars.v(end,1), 0.4, 't=0', 'fontsize', 16); 
    hold off;
    
    subplot(3,4,[7,11]);
    plot(100*tdvars.phi_g(:,1), z, 'k:', 100*tdvars.phi_g(:,ti), z, 'r-');
    set(gca, 'ydir', 'reverse');
    
    subplot(3,4,[8,12]);
    plot(tdvars.mw(:,1), z, 'k:', tdvars.mw(:,ti), z, 'r-');
    set(gca, 'ydir', 'reverse');
    
    suptitle(['t = ' num2str(ConvertSecToYear(td.x(ti)), 2), ' year']);
    
    currFrame = getframe(figHand);
    writeVideo(vidObj,currFrame);
    
end

close(vidObj);

end