function [GasNoDate, CO2, CO2Err, SO2, SO2Err] = LoadGasEmissions ()

T = readtable('GasEmissionsTimeSeries.xlsx');
ErrFrac = 0.2;

date = table2array(T(:,1));
GasDate = datenum(datestr(date,'dd-mmm-yy',1900));
StartDate = datenum(datestr(EruptionStartDate,'dd-mmm-yy',1900));
GasNoDate = 1/365*(GasDate - StartDate);

CO2    = table2array(T(:,2));
CO2Err = ErrFrac*table2array(T(:,2));
SO2    = table2array(T(:,4));
SO2Err = ErrFrac*table2array(T(:,4));

end