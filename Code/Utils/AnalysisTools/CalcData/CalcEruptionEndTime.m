function [tEnd, Ind] = CalcEruptionEndTime (td, m, VelThres)
% crude way of estimating the end of an eruption based on the exit 
% velocity going below some threshold (VelThres)

tEnd = nan; 

tdvars = extract_y(td,m);

Ind = find(tdvars.v(end,:) < VelThres, 1);

if ~isempty(Ind)
    tEnd = td.x(Ind);
end


end