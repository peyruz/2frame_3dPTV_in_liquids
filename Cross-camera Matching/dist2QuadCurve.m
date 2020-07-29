% Calculate the distances of given points to a given quadratic curve
% segment on a plane

% Input
% p - array of points (N,2), where N is the number of points
% a, b, c - coefficients of the curve given as y=a*x^2+b*x+c
% x1, x2 - lower and upper boundaries in x of the curve segment

% Output
% dist - vector of distances of the points to the curve segment

% Note, since a segment is considered, if the closest point on the curve is
% outside of the given boundaries, then one of the boundary points will be
% the closest point.
% if unbounded curve is to be used, x1 and x2 can be set to -inf and inf
% respectively
% Because of the way the curve is parameterized, the curve has to represent
% a bijective or surjective x->y mapping.

% Example usage:
% [dist] = dist2QuadCurve(p,a,b,c,x1,x2)
% [dist, xint, yint] = dist2QuadCurve(p,a,b,c,x1,x2)

% by Peyruz Gasimov, Sep 2017

function [dist, varargout] = dist2QuadCurve(p,a,b,c,x1,x2)

% Find the closest point on the curve (curve boudaries x1 and x2 ignored)
xint =((p(:,1) - b.*c + b.*p(:,2))./(4.*a.^2) - b.^3./(8.*a.^3) + (((b.^2 + 2.*a.*c - 2.*a.*p(:,2) + 1)./(6.*a.^2)...
    - b.^2./(4.*a.^2)).^3 + ((p(:,1) - b.*c + b.*p(:,2))./(4.*a.^2) - b.^3./(8.*a.^3) + (b.*(b.^2 + 2.*a.*c - 2.*a.*p(:,2) + 1))./(8.*a.^3)).^2).^(1./2)...
    + (b.*(b.^2 + 2.*a.*c - 2.*a.*p(:,2) + 1))./(8.*a.^3)).^(1./3) - b./(2.*a) - ((b.^2 + 2.*a.*c - 2.*a.*p(:,2) + 1)./(6.*a.^2)...
    - b.^2./(4.*a.^2))./((p(:,1) - b.*c + b.*p(:,2))./(4.*a.^2) - b.^3./(8.*a.^3) + (((b.^2 + 2.*a.*c - 2.*a.*p(:,2) + 1)./(6.*a.^2) - ...
    b.^2./(4.*a.^2)).^3 + ((p(:,1) - b.*c + b.*p(:,2))./(4.*a.^2) - b.^3./(8.*a.^3) + (b.*(b.^2 + 2.*a.*c - 2.*a.*p(:,2) + 1))./(8.*a.^3)).^2).^(1./2)...
    + (b.*(b.^2 + 2.*a.*c - 2.*a.*p(:,2) + 1))./(8.*a.^3)).^(1./3);
    
% Apply curve boundaries
xint(xint<x1)=x1;
xint(xint>x2)=x2;

% Calculate y of intersection points
yint=a*xint.^2+b*xint+c;

% Preallocate
dist=zeros(size(p,1),1);

% Calculate distances to the curve segment
for ii=1:size(p,1)
dist(ii)=norm([xint(ii)-p(ii,1),yint(ii)-p(ii,2)]);
end

% Variable output
if nargout==3
    varargout{1}=xint;
    varargout{2}=yint;
end

end