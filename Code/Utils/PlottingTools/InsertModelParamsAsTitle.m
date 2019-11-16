function [] = InsertModelParamsAsTitle (o)

[FldNames, FldVals] = ConvertModelParamsForPlotting(o);

TitleText = '';

for ivar = 1:length(FldVals)
    TitleText = [TitleText, FldNames{ivar}, ' = ', num2str(FldVals(ivar),3), '.  '];
    
    if ivar == round(length(FldVals)/2)
        TitleText = [TitleText, '\newline'];
    end
end

suptitle(TitleText);



end