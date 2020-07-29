% Find particle center via a 2d Gaussian fit
% 
% Following Mann J, Ott S, Andersen JS (1999) Experimental study of relative,
%           turbulent diffusion. Risø National Laboratory Report Risø-R-
%           1036(EN) pp 19-20
% 
% Summary:
% The method is useful when more than one particle is contained in a
% connected image object. This is determined by counting intensity maxima.
% The number of maxima is assumed to be the number of objects contained in
% the image. The image is the approximated a summation of n Gaussians,
% where n is the number of maxima. The initial guesses need to be adjusted
% by the user based on the particular usage case. 
%
% Notable Assumptions:
%   1.  The PSF can be adequately modeled by a Gaussian function.
%   2.  In this implementation I assume that sigmaX and sigmaY are same for all
%       particles. This assumption is supported by the fact that the particle
%       tracers used are quite monodisperse. 
%   3.  Similarly, due to the fact that the
%       particles belonging to the same connected object are obviously close to
%       each other in physical space (provided the thickness of the model is
%       limited as is the case in fracture tracer flow studies) and thus
%       experience similar amount of light and thus should have similar I0
%       or the total emmited light intensity.
%   4.  The noise can be modelled by a constant (b)
%
% Input:
% im - image containing the connected object
% 
% Output:
% xc, yc - center of the fit Gaussian
% varargout{1} - (optional) number of peaks in the supplied image
% 
% Example usage:
%   [xc,yc]=GaussFit2d(im);
% Additionally output fit error of the Gaussians
%   [xc,yc, fitError]=GaussFit2d(im);
% Additionally output the number of peaks found
%   [xc,yc, fitError, peakNum]=GaussFit2d(im);
% 
% by Peyruz Gasimov, Oct 28 2017

function [xc,yc, varargout]=GaussFit2d(im)


% Find intensity peaks
conn=8;
peakMask = imregionalmax(im,conn);
peakNum=sum(peakMask(:));

% Generate Weight matrix
latW1=0.8;
latW2=0.2;
Weights=latW2*ones(size(im));
[id,jd]=find(peakMask);
idt=id;
idt(idt~=1)=idt(idt~=1)-1;
idb=id;
idb(idb~=size(im,1))=idb(idb~=size(im,1))+1;
jdl=jd;
jdl(jdl~=1)=jdl(jdl~=1)-1;
jdr=jd;
jdr(jdr~=size(im,2))=jdr(jdr~=size(im,2))+1;
Weights(sub2ind_2d(size(im,1),idt,jd))=latW1;
Weights(sub2ind_2d(size(im,1),idt,jdl))=latW1;
Weights(sub2ind_2d(size(im,1),idt,jdr))=latW1;
Weights(sub2ind_2d(size(im,1),id,jdl))=latW1;
Weights(sub2ind_2d(size(im,1),id,jdr))=latW1;
Weights(sub2ind_2d(size(im,1),idb,jdl))=latW1;
Weights(sub2ind_2d(size(im,1),idb,jd))=latW1;
Weights(sub2ind_2d(size(im,1),idb,jdr))=latW1;
Weights(peakMask)=1;

% Initial guesses
I0ini=1000;                     % the total luminenscence a standalone particle
sigmaXini=1;                    % stddev of the Gaussians
sigmaYini=1;                    
[ycini,xcini]=find(peakMask);   % center coordinates
bini=0;                         % constant intensity bias

xini=zeros(peakNum,1);
xini(1)=I0ini;
xini(2)=sigmaXini;
xini(3)=sigmaYini;
xini(4)=bini;
xini(5:5+peakNum-1)=xcini;
xini(5+peakNum:5+2*peakNum-1)=ycini;

% It is important as it turned out to allow the center locations range at least -1
% and +1 from the middle of the pixel as opposed to -0.5:+0.5 as would seem
% intuitively. In other words, the brightest pixel does not necessarily
% contain the center a particle and instead may be a result of summation of
% two adjacent particles' intensities.
lb=[0 0 0 0 xcini'-1 ycini'-1]';
ub=[6000 100 100 1000 xcini'+1 ycini'+1]';

% Calculate matrix indices for the image for gaussian calculation (see
% below)
[mixX, mixY]=meshgrid(1:size(im,2), 1:size(im,1));

% Fit to the image
optionsLsq=optimoptions('lsqnonlin','Display','none');

fun=@(x) gaussFitObjfun(x);

[x, fitError] = lsqnonlin(fun,xini,lb,ub,optionsLsq);

xc=x(5:5+peakNum-1);
yc=x(5+peakNum:5+2*peakNum-1);

% varargout estimation
if nargout>=3
    varargout{1}=fitError;
end

if nargout==4
    varargout{2}=peakNum;
end

% The nested objective function
    function [ObjFunVal]=gaussFitObjfun(x)
        
        I0=x(1);
        sigmaX=x(2);
        sigmaY=x(3);
        b=x(4);
        xc=x(5:5+peakNum-1);
        yc=x(5+peakNum:5+2*peakNum-1);
        
        imEst=zeros(size(im));
        
        for ii=1:peakNum       
            imEst=imEst+...
                I0/(2*pi*sigmaX*sigmaY)*exp(-0.5*(((mixX-xc(ii))/sigmaX).^2+((mixY-yc(ii))/sigmaY).^2)); 
        end
        
        imEst=imEst+b;
        
        ObjFunVal=mat2vec(double(imEst)-double(im).*Weights);
        
    end
end
