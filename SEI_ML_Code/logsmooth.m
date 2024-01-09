function [xOut,yOut]=logsmooth(xIn,yIn,N)
    yOut=[];
    xIn=xIn(find(not(xIn==0)));
    xMax=ceil(log10(max(xIn)));
    xMin=floor(log10(min(xIn)));
    xOut=logspace(xMin,xMax,(xMax-xMin)*N);
    for i=(1:length(xOut)-1)
        yOut=[yOut; mean(yIn(find(and(xIn>xOut(i),xIn<xOut(i+1)))))];        
    end
    xOut=xOut(1:end-1);
end