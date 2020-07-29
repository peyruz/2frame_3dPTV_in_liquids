% Test 3dPTV Processor on  2d Burger's vortex dataset
clear
clc
% Generate a sample of randomly positioned particles
N1=1500;
unmatched=-50;
N2=N1+unmatched;

Width=10;   %Width of the vortex

pos1=Width*rand(N1,2)-0.5*Width*ones(N1,2);
% pos1(:,3)=zeros(N1,1);

r=sqrt(pos1(:,1).^2+pos1(:,2).^2);

XplusYplus=pos1(:,1)>0 & pos1(:,2)>0;
XminusYplus=pos1(:,1)<0 & pos1(:,2)>0;
XminusYminus=pos1(:,1)<0 & pos1(:,2)<0;
XplusYminus=pos1(:,1)>0 & pos1(:,2)<0;

theta=zeros(N1,1);
theta(XplusYplus)=atan(pos1(XplusYplus,2)./pos1(XplusYplus,1));
theta(XminusYplus)=pi-atan(abs(pos1(XminusYplus,2)./pos1(XminusYplus,1)));
theta(XminusYminus)=pi+atan(abs(pos1(XminusYminus,2)./pos1(XminusYminus,1)));
theta(XplusYminus)=atan(pos1(XplusYminus,2)./pos1(XplusYminus,1));

% theta=theta';

Gamma=200;  % mm^2 s^{-1}
Sigma=0.01; % s^{-1}
nu=1;       % mm^2 s^{-1}

dt=0.05;
Nstepsff=5; % Time of time steps forward

Vel(:,1)=(-sin(theta)).*repmat(Gamma, N1,1)./(2*pi*r).*(ones(N1,1)-exp(-Sigma/4*r.^2./repmat(nu,N1,1)));
Vel(:,2)=(cos(theta)).*repmat(Gamma, N1,1)./(2*pi*r).*(ones(N1,1)-exp(-Sigma/4*r.^2./repmat(nu,N1,1)));
pos2=pos1+Vel*dt;


for i=1:Nstepsff;
    
    r=sqrt(pos2(:,1).^2+pos2(:,2).^2);
    
    XplusYplus=pos2(:,1)>0 & pos2(:,2)>0;
    XminusYplus=pos2(:,1)<0 & pos2(:,2)>0;
    XminusYminus=pos2(:,1)<0 & pos2(:,2)<0;
    XplusYminus=pos2(:,1)>0 & pos2(:,2)<0;
    
    theta=zeros(N1,1);
    theta(XplusYplus)=atan(pos2(XplusYplus,2)./pos2(XplusYplus,1));
    theta(XminusYplus)=pi-atan(abs(pos2(XminusYplus,2)./pos2(XminusYplus,1)));
    theta(XminusYminus)=pi+atan(abs(pos2(XminusYminus,2)./pos2(XminusYminus,1)));
    theta(XplusYminus)=atan(pos2(XplusYminus,2)./pos2(XplusYminus,1));
    
    Vel(:,1)=(-sin(theta)).*repmat(Gamma, N1,1)./(2*pi*r).*(ones(N1,1)-exp(-Sigma/4*r.^2./repmat(nu,N1,1)));
    Vel(:,2)=(cos(theta)).*repmat(Gamma, N1,1)./(2*pi*r).*(ones(N1,1)-exp(-Sigma/4*r.^2./repmat(nu,N1,1)));
    pos2=pos2+Vel*dt;
    
end

if unmatched>0                                                          % If the second snapshot has to have more particles than the first
    pos2=[pos2; Width*rand(unmatched,2)-0.5*Width*ones(unmatched,2)];   % We add randomly generated particles to the matched ones
elseif unmatched<0                                                      % If second snapshot had less particles than the first
    pos2=pos2(1:N2,:);                                                  % We crop the second snapshot collection
end

pos2(:,3)=zeros(N2,1);
pos1(:,3)=zeros(N1,1);
