function [beta, beta_mag, beta_ch] = CalcCompressibility (td, m, plot_opt)

tdv = extract_y(td,m);

beta_ch = m.ch.beta_ch*ones(size(td.x));

beta = zeros(size(td.x));
for ti = 1:length(td.x)
    beta(ti) = tdcFV('CalcCompressibility', td.x(ti), tdv.p(1,ti), m);
end

beta_mag = beta - beta_ch;

if (plot_opt)
    figure;
    plot(tdv.tyr, beta, tdv.tyr, beta_mag, tdv.tyr, beta_ch);
    legend('Total Compressibility','Magma Compressibility',...
        'Chamber Compressibility','location','best'); 
    legend boxoff;
end
    
    


end