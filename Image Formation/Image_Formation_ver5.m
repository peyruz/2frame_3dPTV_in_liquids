% Image formation by Steger camera model

% Source:
% A Comprehensive and Versatile Camera Model for Cameras with Tilt Lenses
% by Carsten Steger
% Int J Comput Vis (2017) 123:121–159

% Input
% c1, c2  -     mm, position of the projection center (entrance pupil)
% Av1, Av2
% s -           collection of real-world positions to be imaged
% od -          optical distance, distance from the projection center
%               (entrance pupil) to the sensor
% Rv -          Rotation (alibi) matrices of the virtual untilted plane
% Rs -          Rotation (alibi) matrices of the sensor plane
% es.., ev.. -  mm, basis vectors of the sensor and virtual plane respectively
%               (e.g. es13 - basis unit vector of sensor plane for camera
%               1, direction 3 (z) )
% fv -          mm, principal distance to the virtual image plane (in world RF)
%               
% Rint -        refractive interface flag. (Rint=1 if refractive interface
%               is present at z=0)
% n -           refractive index of the medium (assumed that the model is
%               residing in air)
% 
% Output
% xi -          image coordinates of the imaged world positions (eg xi1 - in camera
%               one)
% 
% Other nomenclature
% fv    -       mm, pricipal distance to the virtual image plane
% od    -       mm, optical distance. Distance from the entrance pupil to
%               the image sensor along the optical axis (if od>1000000, the camera
%               can treated as telecentric see source for details)
% distCenWRF -  mm, position of the distortion center on the sensor plane in
%               world RF. distortion center lies on the line of the intersection
%               of the translated virtual image plane and the sensor plane. This line of
%               intersection if normal to the optical axis.

% Name of variable is concatenated according to g+(u,d)+(v,t,s)+(W,L)
% u,d - undistorted vs. distorted
% v,t,s - virtual vs. tilted and sensor image planes
% W,L - world vs. local (planar) reference frame

% by Peyruz Gasimov
% Oct 2017

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
% - added unquantified distortion based on the real data
%
% ver 5, Oct, 2017
% - transition from pinhole model to thick-lens model
% - explicit inclusion of lens tilt into the model
% - distortion parameters are now variable between cameras
% - unquantified distortion removed

function [xi1,xi2]=Image_Formation_ver5(c1, c2, fv1, fv2, od1,od2, Rs1, Rs2,Rv1,Rv2, es13,es23, ev13, ev23, s,oi1,oi2,RInt,n,LensType)

% String argument check
if ~( strcmp(LensType,'Telc') || strcmp(LensType,'Entc') )
    error('Invalid LensType. Accepted options are: ''Telc'' ot ''Entc''. ')
end

% Distortion parameters (No distortion case)
% Camera 1
k11=0;      % Radial Distortion terms
k12=-0.0;
tau11=0.0;     % Tangential distortion terms
tau12=0.0;
A11=0.0;           % Affinity distortion 
A12=0.0;
B11=0;        % Shear distortion
B12=-0;

% Camera 2 
k21=-0.0;     
k22=-0.0;
tau21=0.0;   
tau22=0.0;
A21=0.0;           
A22=0.0;
B21=0;                                     
B22=-0;

% % Distortion parameters
% % Camera 1
% k11=-0.000002;      % Radial Distortion terms
% k12=-0.00000004;
% tau11=0.000003;     % Tangential distortion terms
% tau12=0.000025;
% A11=0.06;           % Affinity distortion 
% A12=0.04;
% B11=5.56e-4;        % Shear distortion
% B12=-2.09e-4;
% 
% % Camera 2 
% k21=-0.000002;     
% k22=-0.00000004;
% tau21=0.000003;   
% tau22=0.000025;
% A21=0.06;           
% A22=0.04;
% B21=5.56e-4;                                     
% B22=-2.09e-4;

% Pixel sizes, mm
pxs1=[0.0055; 0.0055];
pxs2=[0.0055; 0.0055];

% Calculate the world positions of the virtual principal point and
% the distortion center
Av1=c1+ev13*fv1;
Av2=c2+ev23*fv2;
distCenW1=c1+ev13*od1;
distCenW2=c2+ev23*od2;

Z=size(s,2);

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
    u_i1=(s-repmat(c1,1,Z))./normMat(s-repmat(c1,1,Z));
    u_i2=(s-repmat(c2,1,Z))./normMat(s-repmat(c2,1,Z));
end

%% Ray tracing to the virtual image plane
l_i1=repmat(dot((Av1-c1),ev13),1,Z)./dot(u_i1,repmat(ev13,1,Z));     % Plane-line intersection solved
l_i2=repmat(dot((Av2-c2),ev23),1,Z)./dot(u_i2,repmat(ev23,1,Z));

% Undistorted image coordinates (World RF)
guvW1=repmat(c1,1,Z)+repmat(l_i1,3,1).*u_i1;                           
guvW2=repmat(c2,1,Z)+repmat(l_i2,3,1).*u_i2;

% Undistorted unscaled sensor coordinates (passive transformation into
% virtual image coords)
guvL1=Rv1'*(guvW1-repmat(Av1,1,Z));
guvL2=Rv2'*(guvW2-repmat(Av2,1,Z));

% Radii defined based on undistorted coordinates shifted according to the center of distortion 
r1=sqrt(guvL1(1,:).^2+guvL1(2,:).^2);   
r2=sqrt(guvL2(1,:).^2+guvL2(2,:).^2);

% Add distortion
gdvL1(1,:)=guvL1(1,:)+guvL1(1,:).*(k11*r1.^2+k12*r1.^4)+tau11*(r1.^2+2*guvL1(1,:).^2)+2*tau12*guvL1(1,:).*guvL1(2,:)+A11*guvL1(1,:)+B11*guvL1(2,:);
gdvL1(2,:)=guvL1(2,:)+guvL1(2,:).*(k11*r1.^2+k12*r1.^4)+tau12*(r1.^2+2*guvL1(2,:).^2)+2*tau11*guvL1(1,:).*guvL1(2,:)+A12*guvL1(2,:)+B12*guvL1(1,:);

gdvL2(1,:)=guvL2(1,:)+guvL2(1,:).*(k21*r2.^2+k22*r2.^4)+tau21*(r2.^2+2*guvL2(1,:).^2)+2*tau22*guvL2(1,:).*guvL2(2,:)+A21*guvL2(1,:)+B21*guvL2(2,:);
gdvL2(2,:)=guvL2(2,:)+guvL2(2,:).*(k21*r2.^2+k22*r2.^4)+tau22*(r2.^2+2*guvL2(2,:).^2)+2*tau21*guvL2(1,:).*guvL2(2,:)+A22*guvL2(2,:)+B22*guvL2(1,:);

gdvL1(3,:)=zeros(1,size(gdvL1,2));
gdvL2(3,:)=zeros(1,size(gdvL2,2));
  
% Calculate Image Space projecting rays
gdvW1=Av1+Rv1*gdvL1;    % World RF
gdvW2=Av2+Rv2*gdvL2;

gdtW1=gdvW1-ev13*(od1-fv1);
gdtW2=gdvW2-ev23*(od2-fv2);

if strcmp(LensType,'Telc')
    u_pr1=repmat(ev13,1,Z);
    u_pr2=repmat(ev23,1,Z);
else
u_pr1=gdtW1-c1;
u_pr2=gdtW2-c2;
u_pr1=u_pr1./normMat(u_pr1);
u_pr2=u_pr2./normMat(u_pr2);
end
lts1=distCenW1-gdtW1;
lts1=dot(lts1, repmat(-es13,1,Z));
lts1=lts1./dot(u_pr1,repmat(-es13,1,Z));

lts2=distCenW2-gdtW2;
lts2=dot(lts2, repmat(-es23,1,Z));
lts2=lts2./dot(u_pr2,repmat(-es23,1,Z));

% Project the virtual image onto the tilted sensor plane
gdsW1=gdtW1+lts1.*u_pr1;
gdsW2=gdtW2+lts2.*u_pr2;

%% Generate the pixel version
xi1=(Rs1'*(gdsW1-oi1))./[pxs1; 1];
xi2=(Rs2'*(gdsW2-oi2))./[pxs2; 1];

%% Add Error of finding the Center of Mass
xi1=xi1(1:2,:);     % Crop the size of the arrays from (3,Z) to (2,Z)
xi2=xi2(1:2,:);

% xi1=xi1+normrnd(0,0.01,[2,Z]);
% xi2=xi2+normrnd(0,0.01,[2,Z]);
end
