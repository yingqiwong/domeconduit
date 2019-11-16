function [dydtMat] = CalcRHS (td, m, plot_opt, tquery)

dydtMat = zeros(size(td.y));

for ti = 1:length(td.x)
    dydtMat(:,ti) = tdcFV('des', td.x(ti), td.y(:,ti), td.z', m);
end


if plot_opt == 1

    if nargin == 3, tquery = round(linspace(1,length(td.x), 5)); end
    plotRHS (td, m, dydtMat, tquery);
end

end

function [] = plotRHS (td, m, dydtMat, tquery)

if isempty(tquery)
    tquery = round(linspace(1,length(td.x), 9));
else
    if td.x(tquery(end)) > td.x(end)
        tquery = round(linspace(1,length(td.x), length(tquery)));
    end
end


figure;
set(gcf, 'defaultAxescolororder', parula(length(tquery)+1));
subplot(141); plot(dydtMat(2:m.Nv:end,tquery), td.z'); title('MBE');
ylabel('Depth (m)');
subplot(142); plot(dydtMat(1:m.Nv:end,tquery), td.z'); title('Mass, SL');
subplot(143); plot(dydtMat(3:m.Nv:end,tquery), td.z'); title('Mass, H2O');
subplot(144); plot(dydtMat(4:m.Nv:end,tquery), td.z'); title('Mass, CO2');
legend(num2str(ConvertSecToYear(td.x(tquery)')),'location','southeast'); 
legend boxoff;

suptitle('RHS of equations');

end

