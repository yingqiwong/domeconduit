function [Poro] = LoadPorosityData ()
% loads porosity data from Smith et al (2011) EPSL

T = readtable('PorositySmith2011.xlsx');

Dates = table2array(T(:,2:3));
PoroDates = [datenum(datestr(Dates(:,1),'dd-mmm-yy',1900)), datenum(datestr(Dates(:,2),'dd-mmm-yy',1900))];
StartDate = datenum(datestr(EruptionStartDate,'dd-mmm-yy',1900));
Poro.Date = 1/365.25*(PoroDates - StartDate);
Poro.Data = table2array(T(:,7));
Poro.Err  = table2array(T(:,8));

end