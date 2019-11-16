function [beta, beta_mag, beta_ch] = CalcCompressibilityFromMCOutput (ParamNames, ParamVals, t, pch)

o       = FillFields(ParamNames, ParamVals);
opts    = tdcFV('ss_init', o);
m       = smf_rad_dz('m_init', opts);
m.tdep  = opts.tdep;

beta = zeros(length(t),1);

for ti = 1:length(t)
    beta(ti) = tdcFV('CalcCompressibility', t(ti), pch(ti), m);
end

beta_ch  = m.ch.beta_ch;
beta_mag = beta - beta_ch;

end