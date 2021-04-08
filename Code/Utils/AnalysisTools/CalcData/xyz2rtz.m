function [rtz] = xyz2rtz (xyz, theta)
% rotates coordinate system from Cartesian to cylindrical about a
% pre-defined center.  Theta is the angle relative to the predefined center
%
% YQW, Nov 2018

[Nsta,Ndir] = size(xyz);
Nth         = length(theta);

if Nth == 1
    rtz = (RotationMatrix(theta,Ndir)*xyz(:,1:Ndir)')';
else
    rtz = zeros(size(xyz));
    for it = 1:Nth
        rtz(it,:) = (RotationMatrix(theta(it),Ndir)*xyz(it,1:Ndir)')';
    end
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
