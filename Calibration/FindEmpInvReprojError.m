% Inverse Reprojection Error of an Empirical Calibration

function [invEmpReprojErr]=FindEmpInvReprojError(CalPhys,PolyCalSystem,PolyCoefx,PolyCoefy)

invReproj(:,1)=PolyCalSystem*PolyCoefx;
invReproj(:,2)=PolyCalSystem*PolyCoefy;

invEmpReprojErr=sqrt((CalPhys(:,1)-PolyCalSystem*PolyCoefx).^2+(CalPhys(:,2)-PolyCalSystem*PolyCoefy).^2);

fprintf('Mean Inverse Reprojection Error: %2.4f mm \n',mean(invEmpReprojErr));

% Visualization
tri=delaunay(CalPhys(:,1),CalPhys(:,2));
invReprojErrFig(invReproj(:,1),invReproj(:,2),invReproj(:,1)-CalPhys(:,1),invReproj(:,2)-CalPhys(:,2),...
                   tri, invEmpReprojErr, invEmpReprojErr)