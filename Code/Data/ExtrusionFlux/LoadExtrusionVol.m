function [FluxNoDate, CumVol, CumVolErr] = LoadExtrusionVol ()

T = readtable('ExtrusionFluxData.xlsx');

date = table2array(T(:,1));
FluxDate = datenum(datestr(date,'dd-mmm-yy',1900));
StartDate = datenum(datestr(EruptionStartDate,'dd-mmm-yy',1900));
FluxNoDate = 1/365.25*(FluxDate - StartDate);

CumVol = table2array(T(:,6));
CumVolErr = table2array(T(:,7));
ExtrRate = table2array(T(:,4));



end