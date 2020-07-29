% Triangulation
%
% Vectorized version of the original code.
%
% Provided that we know the location of the centers of perspective from
% calibration, we may reconstruct the rays directed from the centers
% to the z=0 points of the dewarped image,
% Further we may find the refracted rays and find the optimal point of
% their intersection.

% Nomenclature
% c - center of perspective
% u_i - unit vector of the incident ray
% u_t - unit vector of the refracted ray

% Input
% p01, p02 - positions of the particle on the dewarped (with polynomial) image
% n - refractive index of the medium


function [OP,varargout]=TriangulateFromDewarpedImageVec(p01,p02,c1p,c2p,e3,n,dRF)

sp=size(p01,2);

% Adjust for difference in Reference Frames
if size(p01,1)==3
    p01=p01(1:2,:);
    p02=p02(1:2,:);
end

p01=[p01;repmat(dRF,1,sp)];
p02=[p02;repmat(dRF,1,sp)];

c1p=c1p+dRF;
c2p=c2p+dRF;

u_i1=(p01-repmat(c1p,1,sp))./repmat(normMat(p01-repmat(c1p,1,sp)),3,1);
u_i2=(p02-repmat(c2p,1,sp))./repmat(normMat(p02-repmat(c2p,1,sp)),3,1);

l1=repmat(dRF,1,sp)./dot(-u_i1,repmat(e3,1,sp));
l2=repmat(dRF,1,sp)./dot(-u_i2,repmat(e3,1,sp));

p01=p01+repmat(l1,3,1).*u_i1;
p02=p02+repmat(l2,3,1).*u_i2;


% if size(p01,1)==2
%      p01(3,:)=zeros(1,sp);
%      p02(3,:)=zeros(1,sp);
% end

u_i1=(p01-repmat(c1p,1,sp))./repmat(normMat(p01-repmat(c1p,1,sp)),3,1);     % Find the unit incident ray
u_i2=(p02-repmat(c2p,1,sp))./repmat(normMat(p02-repmat(c2p,1,sp)),3,1);


u_t1=1/n*u_i1+repmat((1/n*dot(-u_i1,repmat(e3,1,sp))-sqrt(1-(1/n)^2*(1-(dot(-u_i1,repmat(e3,1,sp))).^2))),3,1).*repmat(e3,1,sp); % Find the unit refracted ray from Snell's law
u_t2=1/n*u_i2+repmat((1/n*dot(-u_i2,repmat(e3,1,sp))-sqrt(1-(1/n)^2*(1-(dot(-u_i2,repmat(e3,1,sp))).^2))),3,1).*repmat(e3,1,sp);

u_r3=cross(u_t1,u_t2)./repmat(normMat(cross(u_t1,u_t2)),3,1);   % Vector of the perpendicular between the two crossing refracted rays
l3=dot((p02-p01),u_r3);                                         % Length of the perpendicular between the two crossing refracted rays

l1=(repmat(-l3,3,1).*cross(u_r3,u_t2)+cross(p02,u_t2)-cross(p01,u_t2))./(u_r3.*repmat(normMat(cross(u_t1,u_t2)),3,1));
l1=l1(2,:);
l2=(repmat(l1(1,:),3,1).*u_t1-p02+p01+repmat(l3,3,1).*u_r3)./u_t2;
l2=l2(2,:);

dc=normMat(p01-repmat(c1p,1,sp))./(normMat(p01-repmat(c1p,1,sp))+normMat(p02-repmat(c2p,1,sp)));    % Decentering factor

OP=p01+repmat(l1,3,1).*u_t1+u_r3.*repmat(l3,3,1).*repmat(dc,3,1);       % Optimal position of the particle. The distance of the position 
                                % from each of the rays is weighted by the correspondent distance to the projection centers.
                                % This way we try to compensate for the angular error
varargout{1}=l3;                   

end