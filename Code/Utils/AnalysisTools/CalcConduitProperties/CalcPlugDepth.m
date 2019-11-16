function PlugDepth = CalcPlugDepth (td, m)

Vel = {'vv','vf'};
VelParts = CalcConduitProperties(td,m,Vel,0,[]);
vvfrac = VelParts.vv./(VelParts.vv + VelParts.vf);

PlugDepth = zeros(1,length(td.x));

for ti = 1:length(td.x)
    PlugDepthInd = find(vvfrac(:,ti) < m.newtonian.vvfrac_thr,1);
    PlugDepth(ti) = td.z(PlugDepthInd);
end
    

end