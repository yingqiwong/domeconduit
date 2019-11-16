function [VolRate] = LoadExtrusionRate ()
% loads dense rock equivalent of extrusion rate in m3/s

T = readtable('ExtrusionFluxData.xlsx');

date = table2array(T(:,1));
FluxDate = datenum(datestr(date,'dd-mmm-yy',1900));
StartDate = datenum(datestr(EruptionStartDate,'dd-mmm-yy',1900));
FluxNoDate = 1/365.25*(FluxDate - StartDate);

Vol = [0;table2array(T(:,5))];
Err = [0;table2array(T(:,7))];
Date2 = [0;FluxNoDate];

ExtrRate = zeros(length(FluxNoDate),1);
ExtrRateErr = zeros(size(ExtrRate));

for ti = 2:length(FluxNoDate)-1
    ExtrRate(ti-1) = 1e6/(365.25*24*3600)*0.9*(Vol(ti+1) - Vol(ti-1))./(Date2(ti+1) - Date2(ti-1));
    ExtrRateErr(ti-1) = 1e6/(365.25*24*3600)*0.9*(Err(ti+1) + Err(ti-1))./(Date2(ti+1) - Date2(ti-1));
end

VolRate.Date = FluxNoDate(1:end-2);
VolRate.ExtrRate = ExtrRate(1:end-2);
VolRate.ExtrRateErr = ExtrRateErr(1:end-2);

end