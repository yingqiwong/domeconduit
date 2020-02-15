function [tdOut, mOut] = RecalcSolutionAtNewTimes (td, m, tNew)
% tNew should be in seconds

tdOut = td;
mOut  = m;

tdOut.x = tNew;
tdOut.y = interp1(td.x, td.y', tdOut.x, 'pchip', 'extrap')';

end