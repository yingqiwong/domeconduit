function E = CalcFieldsAtOneTimeStep (td, t, m)
% can be used after integration to calculate parameters at each depth
% E = PlotResults('AddFields', td, t, m);

ti  = find(td.x>=t,1);
if isempty(ti)
    fprintf('Solution at time %.1e sec not evaluated\n',t); 
    E = [];
    return; 
end

Nvz = m.Nv*m.Nz;
yscl = td.y(1:Nvz,ti).*repmat(m.slv.sy(1:m.Nv),length(td.z),1);

is = m.blk.is;
p = yscl(is.p:m.Nv:end);

% calculate pressure gradient across faces.
% this has length Nz-1 since we only have (Nz-1) interior faces
dz = td.z(2) - td.z(1);
dpdz = 1/dz*(p(2:end) - p(1:end-1));
dpdzn = [dpdz; 2*dpdz(end) - dpdz(end-1)];
dpdzs = [2*dpdz(1)-dpdz(2); dpdz];

E = tdcFV('calc_exprs',td.x(ti),yscl,td.z',m,dpdz,1);

E.p     = p;
E.v.v   = yscl(is.v:m.Nv:end);
E.phi_g = yscl(is.phi_g:m.Nv:end);
E.mw    = yscl(is.mw:m.Nv:end);
E.pp    = dpdz;

[kmagn, kmags] = tdcFV('get_FaceVals',E.kvert/m.eta_g, m.dQUICKn);
E.vgn = -kmagn.*dpdzn;
E.vgs = -kmags.*dpdzs;
E.dpdz = dpdz;

z = td.z';
E.zMBE = z(1:end-1) + 0.5*(z(2)-z(1));
end
