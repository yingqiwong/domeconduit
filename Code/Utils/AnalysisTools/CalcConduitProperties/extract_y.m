function [tdvarsVec] = extract_y (tdvec, mvec)

for mi = 1:length(tdvec)
    
    td = tdvec(mi);
    m  = mvec(mi);
    
    Nvz = m.Nv*m.Nz;
    
    is = m.blk.is;
    
    tdvars = td;
    tdvars.tyr = ConvertSecToYear(td.x);
    
    tdvars.p     = td.y(is.p:m.Nv:Nvz,:)*m.slv.sy(is.p);
    tdvars.v     = td.y(is.v:m.Nv:Nvz,:)*m.slv.sy(is.v);
    tdvars.phi_g = td.y(is.phi_g:m.Nv:Nvz,:)*m.slv.sy(is.phi_g);
    tdvars.mw    = td.y(is.mw:m.Nv:Nvz,:)*m.slv.sy(is.mw);
    
    if ~isempty(td.yp)
        tdvars.pp     = td.yp(is.p:m.Nv:Nvz,:)*m.slv.sy(is.p);
        tdvars.vp     = td.yp(is.v:m.Nv:Nvz,:)*m.slv.sy(is.v);
        tdvars.phi_gp = td.yp (is.phi_g:m.Nv:Nvz,:)*m.slv.sy(is.phi_g);
        tdvars.mwp    = td.yp(is.mw:m.Nv:Nvz,:)*m.slv.sy(is.mw);
    end
    
    tdvarsVec(mi) = tdvars;
    
end
end