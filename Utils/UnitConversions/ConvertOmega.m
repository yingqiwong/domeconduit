function [OmegaOut, UnitOut] = ConvertOmega (OmegaIn, UnitIn)
% function to convert Omega between units of m3/s/Pa to m3/day/MPa

if strcmp(UnitIn, 'm3/s/Pa')
    
    OmegaOut = OmegaIn*1e6*24*3600;
    UnitOut = 'm3/day/MPa';
    
elseif strcmp(UnitIn, 'm3/day/MPa')
    
    OmegaOut = OmegaIn/1e6/24/3600;
    UnitOut = 'm3/s/Pa';
    
end

end