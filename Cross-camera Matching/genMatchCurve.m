% Generate a Match Curve
% 
% About Match Curves.
% Match curve is a curve on the refractive interface where given an
% intersection of line of sight of Camera 1 we can expect to see the
% intersection of line of sight of Camera 2, provided the lines of sight
% converge in a point below the interface.
% In other words, Match Curve is a prediction of where the z=0 plane image
% of Camera 2 should fall, if we know the z=0 plane image of Camera 1
% particle. Thus, Match Curve is an epipolar line of a kind which due to
% the presense of a refractive interface is not a line but a curve.
% Match curve is used for cross-camera matching of particles. Ideally, the
% z=0 Cam2 image should be exactly on the Match Curve (with its exact
% position along the curve determined by the depths of convergence of lines
% of sight of Camera 1 and 2). In reality, the z=0 Cam2 image will lie at
% some distance from the curve.
% 
% The function calculates a number of points on the Match Curve and fits a
% 2d polynomial to it.
% 
% Input
% z0 -                  coordinates of z=0 image of Camera 1 <2,1>. Supplied from dewarped
%                       image.
% c1,c2 -               position vectors of the two cameras' projection centers
% n1,n2 -               refractive indices of the media. (for air-Polyurethane case,
%                       n1=1; n2= 1.5
% minDepth, maxDepth -  minimum and maximum depths possible to be occupied
%                       by particle. Note that 0>minDepth>maxDepth.
% 
% Output
% PolyCoef -            Polynomial coefficients of the fit quadratic curve
% x1, x2 -              Beginning and end of the curve segment given as x
%                       coordinates
% MatchCurvePoints -    (optional) additionally output the Match Curve
%                       Points used for the fit
% 
% Example usage:
% [PolyCoef,x1,x2]=genMatchCurve(z0,c1,c2,n1,n2,minDepth,maxDepth)
% [PolyCoef,x1,x2,MatchCurvePoints]=genMatchCurve(z0,c1,c2,n1,n2,minDepth,maxDepth)

% by Peyruz Gasimov,Sep 2017

function [PolyFitCoef,x1,x2,varargout]=genMatchCurve(z0,c1,c2,n1,n2,minDepth,maxDepth)
z0=[z0;0];

% All calculations are for refractive interface at z=0 plane.
e3=[0 0 1]';

u_i1=(z0-c1)/norm(z0-c1);  % unit incidence ray vector
u_t1=n1/n2*u_i1+e3*( (n1/n2)*dot(-e3,u_i1)-sqrt(1-(n1/n2)^2*(1-dot(-e3,u_i1)^2)));

% This factor below controls the density of points generated to fit the match curve
DensityFactor=2;    

% The number of fit points will vary based on the inclination of the
% refracted ray
pointNum=round((maxDepth-minDepth)/u_t1(3)*DensityFactor); 

sampleDepths=linspace(minDepth, maxDepth,pointNum);

% Preallocate
samplePositions=zeros(3,pointNum);
MatchCurvePoints=zeros(3,pointNum);

for ii=1:pointNum
samplePositions(:,ii)=z0+( sampleDepths(ii)/u_t1(3) )*u_t1;
end

for ii=1:pointNum
    
% Solve inverse refraction problem
ini2=samplePositions(:,ii)-c2;
ini2=ini2/norm(ini2);   % initial guess for u_i2
ini2=ini2(2);

solx2 = fzero(@(i22)InvRefract(i22,c2,samplePositions(1,ii),samplePositions(2,ii),samplePositions(3,ii),n2/n1), ini2);
u_i22=real(double(solx2));
u_i21=u_i22*(samplePositions(1,ii)-c2(1))/(samplePositions(2,ii)-c2(2));
u_i23=-sqrt(1-u_i21^2-u_i22^2);

u_i2=[u_i21; u_i22; u_i23];

MatchCurvePoints(:,ii)=c2+u_i2*abs(c2(3)/u_i2(3));

end

% Match Curve boundaries
x1=MatchCurvePoints(1,1);
x2=MatchCurvePoints(1,end);

% Fit a quadratic polynomial
[PolyFitCoef]=polyfit(MatchCurvePoints(1,:), MatchCurvePoints(2,:),2);

% Variable Output
if nargout==4
    varargout{1}=MatchCurvePoints(1:2,:);
end

    

