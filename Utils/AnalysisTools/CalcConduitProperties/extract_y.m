function [tdvarsVec] = extract_y (tdvec, mvec)

for mi = 1:length(tdvec)
    
    td = tdvec(mi);
    m  = mvec(mi);
    
    is = m.blk.is;

    tdvars = td;
    tdvars.tyr = ConvertSecToYear(td.x);
    
    tdvars.p     = td.y(is.p:m.Nv:end,:)*m.slv.sy(is.p);
    tdvars.v     = td.y(is.v:m.Nv:end,:)*m.slv.sy(is.v);
    tdvars.phi_g = td.y(is.phi_g:m.Nv:end,:)*m.slv.sy(is.phi_g);
    tdvars.mw    = td.y(is.mw:m.Nv:end,:)*m.slv.sy(is.mw);
    
    tdvars.pp     = td.yp(is.p:m.Nv:end,:)*m.slv.sy(is.p);
    tdvars.vp     = td.yp(is.v:m.Nv:end,:)*m.slv.sy(is.v);
    tdvars.phi_gp = td.yp (is.phi_g:m.Nv:end,:)*m.slv.sy(is.phi_g);
    tdvars.mwp    = td.yp(is.mw:m.Nv:end,:)*m.slv.sy(is.mw);
    
    tdvarsVec(mi) = tdvars;
    
end
end