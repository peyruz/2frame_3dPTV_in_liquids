% Circular Control bias removal
%
% The algorithm follows the idea described in:
% A Comprehensive and Versatile Camera Model for Cameras with Tilt Lenses
% by Carsten Steger
% Int J Comput Vis (2017) 123:121–159
%
% Input:
% CalCam                -     	<N,2> biased dot image centroids
% CalPhys               -     	<N,2> true position of dots in world RF
% Contours              -      	{1,K}{1,N}(M,2) a cell array (K - number of images,
%                               M - number of pixels containing the contour).
% iterNum               -       <1> number of iterations of bias correction (use 1, Headnote [1])
%
% Output:
% CalCamCorr            -      <N,2> corrected, unbiased dot image centers
% contValidityMaskAll   -      <size(CalCamMeanValid)>
% invPolyCoef -     unbiased inverse calibration polynomial coefficients
% PolyCoef -        unbiased direct calibration polynomial coefficients
%
% Headnote [1]  :   Using only one iteration is recommended because in
%                   absolute majority of cases it is enough. Moreover,
%                   using more iterations may actually lead to detrimental
%                   results, which may be explained by the errors of the
%                   empirical calibration.
%                   
% by Peyruz Gasimov. Oct, 2017


function [CalCamCorr,contValidityMaskAll, invPolyCoefxCorr, invPolyCoefyCorr, PolyCoefxCorr, PolyCoefyCorr]=UnbiasedEmpCal(CalCamMeanValid, CalPhys, Contours, iterNum)

imNum=length(Contours);

%% Initial (biased) empirical calibration
[PolyCalSystemCent]=invPolyCal2d7ver2(CalCamMeanValid(:,1),CalCamMeanValid(:,2));
invPolyCoefx=(PolyCalSystemCent\eye(size(PolyCalSystemCent,1)))*CalPhys(:,1);
invPolyCoefy=(PolyCalSystemCent\eye(size(PolyCalSystemCent,1)))*CalPhys(:,2);

PolyCalSystemW0=PolyCal2d7ver2(CalPhys(:,1),CalPhys(:,2));
% PolyCoefx=(PolyCalSystemW0\eye(size(PolyCalSystemW0,1)))*CalCamMeanValid(:,1);
% PolyCoefy=(PolyCalSystemW0\eye(size(PolyCalSystemW0,1)))*CalCamMeanValid(:,2);

CalCamCorr=zeros(size(CalCamMeanValid));
PolyCalSystemCont=cell(imNum,1);
% If the detected contour consists of less than 'minDotNum'  points, we
% ignore it and assign nan to the circle center
minDotNum=10;

%% Bias correction
contValidityMask=cell(1,imNum);

for kk=1:iterNum
    a=cell(1,1,imNum);
    b=cell(1,1,imNum);
    
    % Find the biased centers' projections (ixCol takes care of sorting)
    for jj=1:imNum
        
        for ii=1:size(CalCamMeanValid,1)
            
            if max( size( Contours{jj}{ii} ) ) <=minDotNum || all(all(isnan(Contours{jj}{ii})))
                
                a{jj}(ii)=nan;
                b{jj}(ii)=nan;
                
                if  kk==1
                    contValidityMask{jj}(ii)=false;
                    fprintf('Poor quality edges: Image# %i, particle# %i. Ignored. \n', jj, ii);
                end
                continue
                
            else
                if iterNum>1 % Speed-up, sacrificing memory, for multiple iterations
                    if kk==1
                        PolyCalSystemCont{jj}{ii}=invPolyCal2d7ver2(Contours{jj}{ii}(:,1),Contours{jj}{ii}(:,2));
                        contValidityMask{jj}(ii)=true;
                    end
                    % Fit a circle to the contour mapped into the world RF.
                    circW = CircleFitByPratt([PolyCalSystemCont{jj}{ii}*invPolyCoefx PolyCalSystemCont{jj}{ii}*invPolyCoefy]);
                    
                elseif iterNum==1 % Save memory (and slightly speed-up) for single-iteration runs
                    PolyCalSystemCont=invPolyCal2d7ver2(Contours{jj}{ii}(:,1),Contours{jj}{ii}(:,2));
                    contValidityMask{jj}(ii)=true;
                    % Fit a circle to the contour mapped into the world RF.
                    circW = CircleFitByPratt([PolyCalSystemCont*invPolyCoefx PolyCalSystemCont*invPolyCoefy]);
                end
                
                % Record the circle center
                a{jj}(ii)=circW(1);
                b{jj}(ii)=circW(2);
                
            end
        end
    end
    
    contValidityMaskAll=all(cell2mat(contValidityMask'),1);
    % Find the mean of the centers accross the images (reduce the effect of
    % image noise)
    a=mean(cell2mat(a),3);
    b=mean(cell2mat(b),3);
    
    % Point 'be' in physical space corresponds to the centroid of the particle in the image space.
    be=[-a(contValidityMaskAll)'+2*CalPhys(contValidityMaskAll,1) -b(contValidityMaskAll)'+2*CalPhys(contValidityMaskAll,2)];
    PolyCalSystemW4=PolyCal2d7ver2(be(:,1),be(:,2));
    PolyCoefxCorr=(PolyCalSystemW4\eye(size(PolyCalSystemW4,1)))*CalCamMeanValid(contValidityMaskAll,1);
    PolyCoefyCorr=(PolyCalSystemW4\eye(size(PolyCalSystemW4,1)))*CalCamMeanValid(contValidityMaskAll,2);
    
    % From this new map we can find the true images of the dot centers
    CalCamCorr(contValidityMaskAll,1)=PolyCalSystemW0(contValidityMaskAll,:)*PolyCoefxCorr;
    CalCamCorr(contValidityMaskAll,2)=PolyCalSystemW0(contValidityMaskAll,:)*PolyCoefyCorr;
    
    % Calculate a new unbiased inverse mapping
    PolyCalSystemCentCorr=invPolyCal2d7ver2(CalCamCorr(contValidityMaskAll,1),CalCamCorr(contValidityMaskAll,2));
    invPolyCoefxCorr=(PolyCalSystemCentCorr\eye(size(PolyCalSystemCentCorr,1)))*CalPhys(contValidityMaskAll,1);
    invPolyCoefyCorr=(PolyCalSystemCentCorr\eye(size(PolyCalSystemCentCorr,1)))*CalPhys(contValidityMaskAll,2);
    
    % Assign for the next iteration
    if kk~=iterNum
        invPolyCoefx=invPolyCoefxCorr;
        invPolyCoefy=invPolyCoefyCorr;
    end
end


%% Display the reprojection error statistics
FindEmpInvReprojError(CalPhys(contValidityMaskAll,:),PolyCalSystemCentCorr,invPolyCoefxCorr,invPolyCoefyCorr);


%% Plot the contour inverse projection into the world RF (biased and
% corrected) and compare against the true contours (by viscircles)
figure
oo=1;   % Image number
mm=30;  % Particle number
while isnan(Contours{oo}{mm})
    mm=mm+1;
    if mm==length(Contours{oo})
        error('Invalid contours. Execution stopped.');
    end
end

PolyCalSystemCont_try=invPolyCal2d7ver2(Contours{oo}{mm}(:,1),Contours{oo}{mm}(:,2));
% Draw the circle of the calibration dot (0.5mm radius)
viscircles(CalPhys(mm,1:2),0.5);
hold on
% plot the biased contours
scatter(PolyCalSystemCont_try*invPolyCoefx,PolyCalSystemCont_try*invPolyCoefy);
% plot the unbiased contours
scatter(PolyCalSystemCont_try*invPolyCoefxCorr,PolyCalSystemCont_try*invPolyCoefyCorr);
axis equal

end



