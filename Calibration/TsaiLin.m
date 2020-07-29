% This function is calculating the initial estimates using simplified linear systems for Tsai calibration 
% parameters which are later to be fed into the non-linear minimization scheme.
% The code is following "Tsai's camera calibration method revisited" by Berthold Horn (2000), based on Tsai (1987)

% Peyruz Gasimov
% Piri Research Group
% June, 2016

% Input
% xi - raw image coordinates (pixels)
% s - known dot positions in physical space
% ps - pixel size, mm
% Note: the correspondence of xi and s columns is n eeded.


function [fini,skini,tini,qini]=TsaiLin(xi,s, pxsxini,pxsyini,x0ini)

%% Preprocessing of the input. 
% To run the calibraion we need to normalize the 
% image coordinates and convert them into mm. We assume the pixel size
% indicated in the camera documentation is true.

gsdn=(xi(1:2,:)-repmat(x0ini,1, size(xi,2))).*repmat([pxsxini;pxsyini],1, size(xi,2)); 

%% In this first linear step we estimate the Rotation Matrix (DCM) and tx and ty
% of the translation vector. We then convert the DCM to an equivalent
% quaternion.

for i=1:size(gsdn,2)
    Alin1(i,1)=s(1,i)*gsdn(2,i);
    Alin1(i,2)=s(2,i)*gsdn(2,i);
    Alin1(i,3)=s(3,i)*gsdn(2,i);
    Alin1(i,4)=gsdn(2,i);
    Alin1(i,5)=-s(1,i)*gsdn(1,i);
    Alin1(i,6)=-s(2,i)*gsdn(1,i);
    Alin1(i,7)=-s(3,i)*gsdn(1,i);

    Blin1(i,1)=gsdn(1,i);
end

xlin1=pinv(Alin1)*Blin1;    % Linear regression

c=1/sqrt(xlin1(5)^2+xlin1(6)^2+xlin1(7)^2);  % Scale factor
skini=c*sqrt(xlin1(1)^2+xlin1(2)^2+xlin1(3)^2); % skew factor

xlin1=xlin1*c;

% Recover translation vector. 

% ty was assumed to be =1 when solving the above
% system of eqs. The recovered scale factor reveals how much off we are
% with this assumption. However due to the fact that c>0 in all cases, we
% might have sign ambiguity. In this particular case ty=c
% We temporarily assign zero to tz

tini=[xlin1(4)/skini; c; 0];

a=[xlin1(1)/skini; xlin1(2)/skini; xlin1(3)/skini]; 
b=[xlin1(5)/skini; xlin1(6)/skini; xlin1(7)/skini];

% Forcing orthonormality on the estimated basis
syms k;
eqn=0==dot(a,b)+k*(dot(a,a)+dot(b,b))+k^2*dot(a,b);
solx = solve(eqn,k);
solx=double(solx);
[~, kix]=min(abs(solx));

a1=a+solx(kix)*b;
b1=b+solx(kix)*a;

a=a1/norm(a1);
b=b1/norm(b1);

Rtini=[a';b';cross(a',b')];

% Calculating quaternion from the rotation matrix (requires Quaternion
% class available at the MATLAB File Exchange site)

qini=quaternion.rotationmatrix(Rtini);

%% Second linear step. Estimate the principal distance and the rest of the tranlation vector.

for i=1:size(gsdn,2)
    Alin2(i,1)=(skini*Rtini(1,1)+Rtini(2,1))*s(1,i)+(skini*Rtini(1,2)+Rtini(2,2))*s(2,i)+(skini*Rtini(1,3)+Rtini(2,3))*s(3,i)+skini*tini(1)+tini(2);
    Alin2(i,2)=-(gsdn(1,i)+gsdn(2,i));

    Blin2(i,1)=(Rtini(3,1)*s(1,i)+Rtini(3,2)*s(2,i)+Rtini(3,3)*s(3,i))*(gsdn(1,i)+gsdn(2,i));
end

xlin2=pinv(Alin2)*Blin2;

fini=xlin2(1);
tini(3)=xlin2(2);
