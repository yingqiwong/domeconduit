function [] = PlaytdSol (td, m, Nt, filename)

figHand = figure;
set(gcf, 'position', [520 250 1046 643], 'color', 'w');

% open movie file
vidObj = VideoWriter(filename);
vidObj.FrameRate = 5;
vidObj.Quality = 100;

open(vidObj);

tdv = extract_y(td, m);
t = ConvertSecToYear(td.x);
z = -1e-3*td.z;

% lower sampling rate
if length(Nt)==1
    tNew = linspace(0,td.x(end),Nt);
else
    tNew = Nt;
end
[td2, m2] = RecalcSolutionAtNewTimes(td, m, tNew);
tdv2 = extract_y(td2,m2);
t2 = ConvertSecToYear(td2.x);
v0 = tdv.v(end,1);

zi = find(td.z>-1000,1);

for ti = 1:length(tNew)
    
    subplot(331);
    plot(t, 1e-6*tdv.p(1,:), 'k', t2(ti), 1e-6*tdv2.p(1,ti), 'r*');
    title('Chamber pressure (MPa)'); xlabel('Time (year)');
    
    subplot(332);
    plot(t, 1e3*tdv.v(end,:), 'k', t2(ti), 1e3*tdv2.v(end,ti), 'r*');
    title('Exit velocity (x10^{-3} m/s)'); xlabel('Time (year)');
    
    subplot(333);
    plot(t, 100*tdv.phi_g(end,:), 'k', t2(ti), 100*tdv2.phi_g(end,ti), 'r*');
    title('Exit porosity (%)'); xlabel('Time (year)');
    
    subplot(3,3,[4,7]);
    plot(1e-6*tdv.p(:,1), z, 'k:', 1e-6*tdv2.p(:,ti), z, 'r-');
    text(1e-6*tdv.p(zi,1), 1, ' t=0', 'fontsize', 16); 
    ylabel('Depth (km)'); set(gca, 'ydir', 'reverse');
    xlabel('Pressure (MPa)');
    
    subplot(3,3,[5,8]);
    if v0 == tdv.v(end,1)
        plot(1e3*tdv.v(:,1), z, 'k:'); hold on;
        text(1e3*tdv.v(end,1), 0.4, 't=0', 'fontsize', 16); 
    end
    if tdv2.v(end,ti)/v0<0.1
        v0 = tdv2.v(end,ti);
    end
    plot(1e3*tdv2.v(:,ti), z, 'r-');
    set(gca, 'ydir', 'reverse', 'xlim', [0,1.5*v0*1e3]);
    hold off;
    xlabel('Velocity (x10^{-3} m/s)');
    
    subplot(3,3,[6,9]);
    plot(100*tdv.phi_g(:,1), z, 'k:', 100*tdv2.phi_g(:,ti), z, 'r-');
    text(100*tdv.phi_g(zi,1), 1, ' t=0', 'fontsize', 16);
    set(gca, 'ydir', 'reverse');
    xlabel('Porosity (%)');
   
    suptitle(['t = ' num2str(t2(ti), 2), ' years']);
    
    currFrame = getframe(figHand);
    writeVideo(vidObj,currFrame);
    
end

close(vidObj);

end