function [ObjFunVal] = RefParObjFun(u,sz01p,sz02p,c1p,c2p,n)

c1p=[u(1);u(2);u(3)];
c2p=[u(4);u(5);u(6)];
dRF=u(7);

e3=[0 0 1]';
siz=size(sz01p,2);
sz01p=[sz01p(1:2,:);repmat(dRF,1,siz)];
sz02p=[sz02p(1:2,:);repmat(dRF,1,siz)];

c1p=c1p+dRF;
c2p=c2p+dRF;

u_i1=(sz01p-repmat(c1p,1,siz))./repmat(normMat(sz01p-repmat(c1p,1,siz)),3,1);
u_i2=(sz02p-repmat(c2p,1,siz))./repmat(normMat(sz02p-repmat(c2p,1,siz)),3,1);

l1=repmat(dRF,1,siz)./dot(-u_i1,repmat(e3,1,siz));
l2=repmat(dRF,1,siz)./dot(-u_i2,repmat(e3,1,siz));

sz01pp=sz01p+repmat(l1,3,1).*u_i1;
sz02pp=sz02p+repmat(l2,3,1).*u_i2;

[l3,l3]=TriangulateFromDewarpedImageVec(sz01pp,sz02pp,c1p,c2p,e3,n);

ObjFunVal=sum(l3.^2);

end

