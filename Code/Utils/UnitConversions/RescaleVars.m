function [yOut, mOut] = RescaleVars (yIn, mIn, NewScale)
% rescale values

mOut = mIn;
yDim = yIn.*repmat(mIn.slv.sy,mIn.Nz,1);

mOut.slv.sy = NewScale;
yOut = yDim./repmat(mOut.slv.sy,mIn.Nz,1);

end