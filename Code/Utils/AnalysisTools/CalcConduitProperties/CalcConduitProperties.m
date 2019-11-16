function VarOut = CalcConduitProperties (td, m, VarNameList, plot_opt, tquery)
% VarOut = CalcConduitProperties (td, m, VarNameList, plot_opt, tquery);
% Initialize structure with fields specified in VarNameList


VarOut = [];
Nvars = length(VarNameList);
for vi = 1:Nvars
    VarOut = setfield(VarOut, VarNameList{vi}, []);
end


% loop to calculate the properties desired at all times
for ti = 1:length(td.x)
    
    E = CalcFieldsAtOneTimeStep(td, td.x(ti), m);
    
    for vi = 1:Nvars
        if isfield(E, VarNameList{vi})
            VarOut.(VarNameList{vi})(:,ti) = E.(VarNameList{vi});
        elseif isfield(E.nv, VarNameList{vi})
            VarOut.(VarNameList{vi})(:,ti) = E.nv.(VarNameList{vi});
        elseif isfield(E.h2o, VarNameList{vi})
            VarOut.(VarNameList{vi})(:,ti) = E.h2o.(VarNameList{vi});
        elseif isfield(E.co2, VarNameList{vi})
            VarOut.(VarNameList{vi})(:,ti) = E.co2.(VarNameList{vi});
        elseif isfield(E.v, VarNameList{vi})
            VarOut.(VarNameList{vi})(:,ti) = E.v.(VarNameList{vi});
        end
    end
    
end

if plot_opt == 1
    PlotProps(VarOut, td, m, tquery);
end


end

function [] = PlotProps(Vars, td, m, tquery)

if isempty(tquery)
    tquery = round(linspace(1,length(td.x), 9));
else
    if td.x(tquery(end)) > td.x(end)
        tquery = round(linspace(1,length(td.x), length(tquery)));
    end
end

Fields = fieldnames(Vars);
Nvars = length(Fields);
Ntime = length(tquery);
zCNT = td.z;
zMBE = td.z(1:end-1) + 0.5*(td.z(2) - td.z(1));

figure;
set(gcf,'DefaultAxesColorOrder', parula(length(tquery)+1));  
PlotInd = 0;

for vi = 1:Nvars
    
    PlotInd = PlotInd + 1;
    
    subplot(1,Nvars,PlotInd);
    
    if size(Vars.(Fields{vi}),1) == m.Nz-1 
        z = zMBE;
    elseif size(Vars.(Fields{vi}),1) == m.Nz
        z = zCNT;
    end
    
    plot(Vars.(Fields{vi})(:,tquery),z);
    title(sprintf('%s', Fields{vi}));
    legend(num2str(ConvertSecToYear(td.x(tquery))'), 'location', 'best'); 
    legend boxoff;
    
end

end