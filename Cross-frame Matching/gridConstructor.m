maxX=max(pos1(:,1));
maxY=max(pos1(:,2));
minY=min(pos1(:,2));
minX=min(pos1(:,1));

XSpan=maxX-minX
YSpan=maxY-minY

nodeY=linspace(0,1,50);
nodeX=linspace(0,1,50);

nodeX=minX+nodeX*XSpan;
nodeY=minY+nodeY*YSpan;

[xq,yq]=meshgrid(nodeX,nodeY);