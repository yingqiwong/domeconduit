function tc = CalcTimeConstant(td, m, plot_opt)


tdvars = extract_y(td, m);
vend = tdvars.v(end,:);

tc_ind = find(vend < exp(-1)*vend(1), 1);

% relation between half life time and time constant for decay
if isempty(tc_ind)
    tc = Inf;
else
    tc = td.x(tc_ind);
end

% tc = interp1(vend,td.x,0.37*vend(1));

if nargin > 2 
    figure; 
    set(gcf, 'defaultlinelinewidth', 2);
    tyr = ConvertSecToYear(td.x);
    plot(tyr, vend, '-', tyr, vend(1)*exp(-td.x/tc), '--');
    legend('solution', 'exponential fit'); legend boxoff;
    xlabel('Time (years)'); ylabel('Exit velocity (m/s)');
end

end

