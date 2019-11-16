function [tss] = CalcSteadyStateTime (td, m, tol)

tdvars = extract_y(td, m);

% vex = tdvars.v(end,:);
% vexdiff = diff(vex)./vex(1:end-1);
% tssind = find(abs(vexdiff)<tol, 1);

% find the steady state with the condition that dvdt is less than some
% fraction (tol) of dvdt at t=0
dvdtex = tdvars.vp(end,:);
tssind = find(abs(dvdtex./dvdtex(2))<tol, 1);

if isempty(tssind)
    tss = inf;
else
    tss = td.x(tssind);
end

end
