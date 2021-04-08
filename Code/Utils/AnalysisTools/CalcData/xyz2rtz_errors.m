function [err_rtz] = xyz2rtz_errors (err_xyz, theta)
% rotates errors of coordinate system from Cartesian to cylindrical about a
% pre-defined center.  
% err_xyz is the Nd*3 of data
% Theta is the angle relative to predefined center (either scalar or Ndx1)
%
% YQW, Nov 2018

[Ndata,Ndir] = size(err_xyz);
err_rtz = zeros(size(err_xyz));

if length(theta)==1
    R = RotationMatrix(theta, Ndir);
elseif length(theta)~=Ndata
    fprintf('Size of err and theta not compatible.\n'); return; 
end

for di = 1:Ndata
    
    if length(theta)>1
        R = RotationMatrix(theta(di), Ndir);
    end
    
    sig_xyz = diag(err_xyz(di,:).^2);
    sig_rtz = R*sig_xyz*R';
    err_rtz(di,:) = sqrt(diag(sig_rtz));
end

end

function [Rot] = RotationMatrix (theta, Ndir)

if Ndir == 2
    Rot = [cos(theta), sin(theta); -sin(theta), cos(theta)];
else
    Rot = [cos(theta), sin(theta), 0;...
          -sin(theta), cos(theta), 0;...
           0, 0, 1];
end
end
