

ss2 = smf_rad_dz('add_fields', ss);

dz = abs(ss.z(2) - ss.z(1));

ColWeight = m.g*dz*trapz(ss2.rho);
rhoavg = 1/m.conduit_length*dz*trapz(ss2.rho);
overpressure = ss.m.p_ch - ColWeight;

fprintf('Conduit-averaged density = %.0f kg/m3. \n', rhoavg);
fprintf('Chamber pressure = %.1f MPa. \n', 1e-6*ss.m.p_ch);
fprintf('Magma column weight = %.1f MPa. \n', 1e-6*ColWeight);
fprintf('Overpressure = %.1f MPa. \n', 1e-6*overpressure);
