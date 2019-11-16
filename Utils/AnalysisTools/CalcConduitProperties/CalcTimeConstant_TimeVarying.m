function [tcExitVel, tcChPress] = CalcTimeConstant_TimeVarying (td, m)

% find time constant of exit velocity
tdvars = extract_y(td, m);
tyr = ConvertSecToYear(td.x);

tcExitVel = ConvertSecToYear(-tdvars.v(end,:)./tdvars.vp(end,:));




% find time constant of chamber pressure decay
tcChPress = zeros(size(tcExitVel));
dp = tdvars.p(1,:) - tdvars.p(1,1);
dpdt = tdvars.pp(1,:);

for ti = 2:length(td.x)
    % implicit equation for time constant for pressure decay.
    tmp = fzero(@(tc) dp(ti)/dpdt(ti) + tc - tc*exp(td.x(ti)/tc), ConvertYearToSec(1));
    tcChPress(ti) = ConvertSecToYear(tmp);
end

tcChPress(1) = nan;
tcChPress(tcChPress<0) = nan;

end


