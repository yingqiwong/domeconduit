
%%
AddPaths;
addpath('../../../../ModelTestFiles/');

%%
clear all; 
load JRO1TrialError.mat
%%
isol = 4;
o = GetOFromM(sols(isol).m);
opts = tdcFV('ss_init',o);
des = wrk('run',opts,100); 

%% look at two models that have the same p_top

d1 = des{20};
d2 = des{23};
% d1 = des{31};
% d2 = des{36};

variables = {'p','v','phi_g','m_w'};
pdiff = norm(d1.y(1,end)-d2.y(1,end))/o.p_top

figure(2); clf;

for i=1:4
    subplot(1,4,i);
    plot(d1.y(i,:),-d1.z,'b-',d2.y(i,:),-d2.z,'r--','linewidth',2);
    xlabel(sprintf(variables{i}));
end

legend('low v_{ch}','high v_{ch}','location','south');

figure(3); clf;
subplot(121); 
plot(d1.y(2,:),-d1.z,'b-','linewidth',2);
title('low v_{ch}');

subplot(122);
plot(d2.y(2,:),-d2.z,'r--','linewidth',2);
title('high v_{ch}');

%%
figure(4); clf;
for i=1:4
    subplot(1,4,i);
    
    if i==2
        for j=1:length(des)
            loglog(des{j}.v_ch,des{j}.y(i,end),'b*-'); hold on;
        end
    else
        
        for j=1:length(des)
            semilogx(des{j}.v_ch,des{j}.y(i,end),'b*-'); hold on;
        end
    end
    xlabel('v_{ch}');
    ylabel(sprintf(variables{i}));
    hold off;
    
end




