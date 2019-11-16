function [Flux, Vert, Latl, rhodvdz, vdrhodz] = CalcDrhoDt (td, m, plot_opt, tquery)

% example command
% [Flux, Vert, Latl, rhodvdz, vdrhodz] = CalcDrhoDt (td, m, 1, [1:100:401]);
%
% Evaluates conservation of mass of all species together, i.e. drho/dt.
% But returns individual components of the combined mass balance equation.


% calculate dependent variables
VarList = {'rho','rho_g','phi_g','phi_l','phi_s','c1','kvert', 'ug'};
DepVar = CalcConduitProperties(td, m, VarList, 0, []);

p = td.y(m.blk.is.p:m.Nv:end,:)*m.slv.sy(m.blk.is.p);
v = td.y(m.blk.is.v:m.Nv:end,:)*m.slv.sy(m.blk.is.v);
dz = td.z(2) - td.z(1);

% for calculating vertical gas velocities
dpdz = 1/dz*(p(2:end,:) - p(1:end-1,:));
dpdzn = [dpdz; 2*dpdz(end,:) - dpdz(end-1,:)];
dpdzs = [2*dpdz(1,:)-dpdz(2,:); dpdz];

vgn = -DepVar.kvert/m.eta_g.*dpdzn;
vgs = -DepVar.kvert/m.eta_g.*dpdzs;

% initialize
Flux = zeros(size(DepVar.rho));
Vert = zeros(size(DepVar.rho));
Latl = zeros(size(DepVar.rho));

rhodvdz = zeros(size(DepVar.rho));
vdrhodz = zeros(size(DepVar.rho));

for ti = 1:length(td.x)
    
    % simple extrapolation to get velocity at base of first continuity
    % cell and velocities at center of continuity cells
    vcontface = [2*v(1,ti) - v(2,ti); v(:,ti)];
    [~, vs]   = tdcFV('get_FaceVals', v(:,ti),  m.dQUICKn);
    
    % calculate face values
    [rhon, rhos] = tdcFV('get_FaceVals', DepVar.rho(:,ti),  m.dQUICKn);
    [gasn, gass] = tdcFV('get_FaceVals', DepVar.rho_g(:,ti).*DepVar.phi_g(:,ti), m.dQUICKn);
    
    % break up flux into two parts using chain rule on 
    % d(rho*v)/dz = rho*(dv/dz) + v*(drho/dz)
    rhodvdz(:,ti) = 1/dz*DepVar.rho(:,ti).*(vcontface(2:end) - vcontface(1:end-1));
    vdrhodz(:,ti) = 1/dz*vs.*(rhon - rhos);
    
    % combined terms for all magma species
    Flux(:,ti) = 1/dz.*(rhon.*vcontface(2:end) - rhos.*vcontface(1:end-1));
    Vert(:,ti) = 1/dz.*(gasn.*vgn(:,ti) - gass.*vgs(:,ti));
    Latl(:,ti) = 2*DepVar.rho_g(:,ti).*DepVar.phi_g(:,ti).*DepVar.ug(:,ti)/m.R;
    
    
    if (m.plug_gas_loss)
        % reassign values for modified equations in the plug.  
        
        % density of non-volatile components 
        rhosl = m.rho_s*DepVar.phi_s(:,ti) + m.rho_l*DepVar.phi_l(:,ti).*DepVar.c1(:,ti);
        
        [rhosln, rhosls] =  tdcFV('get_FaceVals', rhosl,  m.dQUICKn);
        [phign, phigs] =  tdcFV('get_FaceVals', DepVar.phi_g(:,ti),  m.dQUICKn);
        
        % flux of non-volatile components in the plug
        pglflux = 1/dz*(rhosln.*vcontface(2:end) - rhosls.*vcontface(1:end-1));
        
        % flux of volatiles in the plug (degassing removes drhog/dt)
        pglgas  = DepVar.rho_g(:,ti)/dz.*(phign.*vcontface(2:end) - phigs.*vcontface(1:end-1));
        
        Flux(m.PlugDepth:end,ti) = pglflux(m.PlugDepth:end);
        Vert(m.PlugDepth:end,ti) = pglgas(m.PlugDepth:end);
        Latl(m.PlugDepth:end,ti) = 0;
    end
    
end
end
