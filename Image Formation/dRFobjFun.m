% Objective function for minimizing dRF - difference between the
% experimental and calibration reference frames.

function g=dRFobjFun(dRF,sz01,sz02,c1p,c2p,n)


e3=[0 0 1]';
siz=size(sz01,2);
sz01=[sz01;repmat(dRF,1,siz)];
sz02=[sz02;repmat(dRF,1,siz)];

c1p=c1p+dRF;
c2p=c2p+dRF;

u_i1=(sz01-repmat(c1p,1,siz))./repmat(normMat(sz01-repmat(c1p,1,siz)),3,1);
u_i2=(sz02-repmat(c2p,1,siz))./repmat(normMat(sz02-repmat(c2p,1,siz)),3,1);

l1=repmat(dRF,1,siz)./dot(-u_i1,repmat(e3,1,siz));
l2=repmat(dRF,1,siz)./dot(-u_i2,repmat(e3,1,siz));

sz01p=sz01+repmat(l1,3,1).*u_i1;
sz02p=sz02+repmat(l2,3,1).*u_i2;

[~,l3]=TriangulateFromDewarpedImageVec(sz01p,sz02p,c1p,c2p,e3,n);

g=sum(l3.^2);