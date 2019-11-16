function [] = AddPlugGasLossDepth (z, depthInd, linestyle)

hold on;
plot(xlim, z(depthInd)*ones(1,2), linestyle, 'color', 'k') 
hold off;
end