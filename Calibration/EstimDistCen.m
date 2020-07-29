% Estimate distortion center of camera

% It is a simple principle used frequently previously We will project into
% world RF an array of straight horizontal and vertical image lines in
% image space using the inverse empirical calibration polynomial we had
% obtained. We will then fit a straight line to the resulting world curves
% and calculate the fitting error. Assuming all distortion sources but
% radial negligible, we can estimate a horizontal and vertical position of
% a line which yields a minimum total fitting error.

% The method is not expected to give an accurate answer, but only a good
% estimate which then can be used in nonlinear regression

% by Peyruz Gasimov, Oct 2017

function [distCenPx]=EstimDistCen(invPolyCoef1x, invPolyCoef1y,xini,yini)

% % Generate points along straight lines in world RF
xh=linspace(1000,6000,150)';
yv=linspace(200,4200,100)';

% Solver options
optionsLsq = optimoptions('lsqnonlin','Display','off',...
                            'TypicalX', 1000,...
                            'MaxFunctionEvaluations',500,'MaxIterations', 100,'StepTolerance',1e-10,'FunctionTolerance',1e-8,'UseParallel',true);
                        
%% Minimize the fit error in horizontal direction
if xini>2500
    lbx=xini-2500;
else
    lbx=0;
end

if xini<(6592-2500)
    ubx=xini+2500;
else
    ubx=6592;
end
                       

fx = @(xv)distCenObjFun_x(xv,yv);

msx = MultiStart('UseParallel',true);
problemx = createOptimProblem('lsqnonlin','objective',fx,'x0',xini,'lb',lbx,'ub',ubx,'options',optionsLsq);

% Number of staring locations
guessNumx=5;

distCenPx(1) = run(msx,problemx,guessNumx);

%% Minimize the fit error in vertical direction

if yini>2000
    lby=yini-2000;
else
    lby=0;
end

if yini<(4400-2000)
    uby=yini+2000;
else
    uby=4400;
end
       
                        
fy = @(yh)distCenObjFun_y(yh,xh);

msy = MultiStart('UseParallel',true);
problemy = createOptimProblem('lsqnonlin','objective',fy,'x0',yini,'lb',lby,'ub',uby,'options',optionsLsq);

% Number of staring locations
guessNumy=5;

distCenPx(2) = run(msy,problemy,guessNumy);

%% Nested functions
% Horizontal direction Objective Function
function [objfunVal]=distCenObjFun_x(xv,yv)

% Image the line
[invPolyCalSystem]=invPolyCal2d7ver2(repmat(xv,size(yv,1),1),yv);
lineW(:,1)=invPolyCalSystem*invPolyCoef1x;
lineW(:,2)=invPolyCalSystem*invPolyCoef1y;

% Fit a straight line to the image
[~,s] = polyfit(lineW(:,1),lineW(:,2),1);

objfunVal=s.normr;
end


% Vertical direction Objective Function
function [objfunVal]=distCenObjFun_y(yh,xh)
    
% Image the line
[PolyCalSystem]=invPolyCal2d7ver2(xh,repmat(yh,size(xh,1),1));
lineW(:,1)=PolyCalSystem*invPolyCoef1x;
lineW(:,2)=PolyCalSystem*invPolyCoef1y;

% Fit a straight line to the image
[~,s] = polyfit(lineW(:,1),lineW(:,2),1);

objfunVal=s.normr;
end
 

end

