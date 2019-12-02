function varargout = ellipse_plot(X, Y, EMAX, EMIN, AZ)
%
%   ellipse_plot(X, Y, EMAX, EMIN, AZ)
%
%   (X,Y)   - Center location of ellipse
%   EMAX    - semi-major axis
%   EMIN    - semi-minor axis
%   AZ      - Angle of rotation
%

np = 100;
angle = [0:np]*2*pi/np;
pts = [X;Y]*ones(size(angle)) + [cos(AZ) -sin(AZ); sin(AZ) cos(AZ)]*[cos(angle)*EMAX; sin(angle)*EMIN];
a = plot( pts(1,:), pts(2,:) );

varargout{1} = a;
