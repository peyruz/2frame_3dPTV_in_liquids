function [xi1, xi2]=Image_Formation_ver4_1(c1, c2, A1,A2, R1, R2, es13,es23, s,oi1,oi2,dcShift,RInt,n)
% The following function calculates the trace of the ray from a particle
% sampling point onto the sensor
%
%
% The modelled environment may be visualized by running the
% PlotCameras.m script
%
% Peyruz Gasimov Rafael
% Piri Research Group
% June 2016
%
% Changelog:
% ver 2, Sep, 2016
% - an additional feature is added wherein particles
%   submersed into a substance (with air interface at z=0) can now be considered
%
% ver 3, Sep 22, 2016
% - case with refractive interface (RInt=1) is parallelized.
% - case of RInt=0 is vectorized
%
% ver 4, Aug, 2017
% - the inverse refraction problem is now solved without using a symbolic
%   toolbox (87x speedup)
% - corrected bug - wrong formula for decentering error
% - added affinity and shear distortion terms
% - added unquatified distortion based on the real data

%% Nomenclature
% u_tj      -    unit vector of the refracted ray for camera j
% The rest of the variables are defined on the way if necessary.

% Observe the consistensy of units. Since the system settigs are defined in
% mm, the matrix of particle positions should also be given in mm.

%% Input
% c1, c2 - positions of the centers of projections for each camera;
% A1, A2 - positions of the principal point in world reference frame
% (RF)
% R1, R2 - rotation matrices representing the rotation of camera coordinate
% system wrt the world RF
% es13, es23 - norms to the image plane
% s - coordinates of the calibration pattern features (dots) in world RF
% oi1, oi2 - coordinates of the image origin for the two camera sensors in
% World RF
% RInt - a boolean variable denoting the presense of refracting interface at z=0

%% Output
% xi1, xi2 - raw image coordinates in pixels
%%

k1=-0.000002;     % Radial Distortion terms
k2=-0.00000004;

tau1=0.000003;   % Tangential distortion terms
tau2=0.000025;

A11=0.06;           % Affinity distortion 
A12=0.04;
B11=5.56e-4;        % Shear distortion
B12=-2.09e-4;

A21=0.06;           % Camera 2 
A22=0.04;
B21=5.56e-4;                                     
B22=-2.09e-4;

% dcShift=[2 0]';  % Difference between the principal point and the center of distortion

Z=size(s,2);

if size(dcShift,1)>2 || size(dcShift,2)>2
    error('Invalid input for Center of Distortion shift. length(dcShift)=2.');
end

if size(dcShift,2)==2
    dcShift=dcShift';
end

%% Reconstruction of incident rays

if RInt==1               % If the interface is present, we will need to numerically solve for the unit incident ray vector
    
    s1=s(1,:);          % We divide the matrix into vectors for smooth operation of the parallel loop
    s2=s(2,:);
    s3=s(3,:);
    
    % To avoid singularity (see the equation of inverse refraction problem)
    s2(s2==c1(2) | s2==c2(2))=s2(s2==c1(2) | s2==c2(2))+0.000001; 
    
    % Need to initialize the arrays somehow for the parfor to work. These values will be cropped later.
    u_i1=[-9999.99;-9999.99;-9999.99];  
    u_i2=[-9999.99;-9999.99;-9999.99];
    
    
    parfor i=1:Z
        ini1=[s1(i); s2(i); s3(i)]-c1;
        ini1=ini1/norm(ini1);   % initial guess for u_i1
        ini1=ini1(2);
        
        ini2=[s1(i); s2(i); s3(i)]-c2;
        ini2=ini2/norm(ini2);   % initial guess for u_i2
        ini2=ini2(2);
        
        solx1 = fzero(@(i12)InvRefract(i12,c1,s1(i),s2(i),s3(i),n), ini1);
        u_i12=real(double(solx1));
        u_i11=u_i12*(s1(i)-c1(1))/(s2(i)-c1(2));
        u_i13=-sqrt(1-u_i11^2-u_i12^2);
        
        
        solx2 = fzero(@(i22)InvRefract(i22,c2,s1(i),s2(i),s3(i),n), ini2);
        u_i22=real(double(solx2));
        u_i21=u_i22*(s1(i)-c2(1))/(s2(i)-c2(2));
        u_i23=-sqrt(1-u_i21^2-u_i22^2);
        
        u_i1=[u_i1,[u_i11; u_i12; u_i13]];
        u_i2=[u_i2,[u_i21; u_i22; u_i23]];
        
    end
    u_i1=u_i1(1:3,2:(Z+1));     % Cropping out the first column
    u_i2=u_i2(1:3,2:(Z+1));
else
    u_i1=(s-repmat(c1,1,Z))/norm(s-repmat(c1,1,Z));
    u_i2=(s-repmat(c2,1,Z))/norm(s-repmat(c2,1,Z));
end

%% Ray tracing to the sensor plane

l_i1=repmat(dot((A1-c1),es13),1,Z)./dot(u_i1,repmat(es13,1,Z));     % Plane-line intersection solved
l_i2=repmat(dot((A2-c2),es23),1,Z)./dot(u_i2,repmat(es23,1,Z));

% Undistorted image coordinates
g1=repmat(c1,1,Z)+repmat(l_i1,3,1).*u_i1;                           
g2=repmat(c2,1,Z)+repmat(l_i2,3,1).*u_i2;

% Undistorted unscaled sensor coordinates
gs1=R1'*(g1-repmat(A1,1,Z));
gs2=R2'*(g2-repmat(A2,1,Z));
gs1=gs1(1:2,:);
gs2=gs2(1:2,:);

% % No distortion case (otherwise - uncomment)
% gsd1=gs1;
% gsd2=gs2;

% Recalculate the unscaled and undistorted coordinates wrt to the center of
% distortion:
gs1cd=gs1+dcShift;
gs2cd=gs2+dcShift;

% Radii defined based on undistorted coordinates shifted accroding to the center of distortion 
r1=sqrt(gs1cd(1,:).^2+gs1cd(2,:).^2);   
r2=sqrt(gs2cd(1,:).^2+gs2cd(2,:).^2);

% Add distortion
gsd1(1,:)=gs1(1,:)+gs1cd(1,:).*(k1*r1.^2+k2*r1.^4)+tau1*(r1.^2+2*gs1cd(1,:).^2)+2*tau2*gs1cd(1,:).*gs1cd(2,:)+A11*gs1cd(1,:)+B11*gs1cd(2,:);
gsd1(2,:)=gs1(2,:)+gs1cd(2,:).*(k1*r1.^2+k2*r1.^4)+tau2*(r1.^2+2*gs1cd(2,:).^2)+2*tau1*gs1cd(1,:).*gs1cd(2,:)+A12*gs1cd(2,:)+B12*gs1cd(1,:);

gsd2(1,:)=gs2(1,:)+gs2cd(1,:).*(k1*r2.^2+k2*r2.^4)+tau1*(r2.^2+2*gs2cd(1,:).^2)+2*tau2*gs2cd(1,:).*gs2cd(2,:)+A21*gs2cd(1,:)+B21*gs2cd(2,:);
gsd2(2,:)=gs2(2,:)+gs2cd(2,:).*(k1*r2.^2+k2*r2.^4)+tau2*(r2.^2+2*gs2cd(2,:).^2)+2*tau1*gs2cd(1,:).*gs2cd(2,:)+A22*gs2cd(2,:)+B22*gs2cd(1,:);

gsd1(3,:)=zeros(1,size(gsd1,2));
gsd2(3,:)=zeros(1,size(gsd2,2));
  

%% Generate the pixel version
xi1=R1'*((R1*gsd1+repmat(A1-oi1,1,Z)))/5.5e-3;
xi2=R2'*((R2*gsd2+repmat(A2-oi2,1,Z)))/5.5e-3;

% % Add unaccounted distortion
% xi1(1,:)=xi1(1,:)+distortInXCam1(xi1(1,:),xi1(2,:));
% xi1(2,:)=xi1(2,:)+distortInYCam1(xi1(1,:),xi1(2,:));
% xi2(1,:)=xi2(1,:)+distortInXCam2(xi2(1,:),xi2(2,:));
% xi2(2,:)=xi2(2,:)+distortInYCam2(xi2(1,:),xi2(2,:));

%% Add Error of finding the Center of Mass
xi1=xi1(1:2,:);     % Crop the size of the arrays from (3,Z) to (2,Z)
xi2=xi2(1:2,:);

xi1=xi1+normrnd(0,0.01,[2,Z]);
xi2=xi2+normrnd(0,0.01,[2,Z]);
end



