% Generate array of physical coors for ReadCalIm function to correspond the
% coordinates of physical calibration pattern with its image coordinates
% This function is a part of software package to calibrate a 3dPTV setup
% Note 1: the number of calibration points needs to be below 6000. Else, is
% the number of dots = x>6000, then set DotNumBound to any number exceeding
% x.
% Note 2: if only a single plane is desired, input 'false' instead of the
% remaining CalDotMasks, e.g. :
% [CalibPatternPhys]=GenPhysDotArray(dg,CalElev,CalDotMask1,false,false,RowNum1,RowNum2,RowNum3)
% will generate the array containing only the mid-plane dot coordinates

% Input:
% dg - grid constant. Spacing between the centers of the dots of the
% calibration pattern
% CalElev - elevations at which the images had been taken. By convention,
% the first image is taken at the postion in between images 2 and 3.
% CalDotMask - logical mask produced by ReadCalIm. Exposes those dots which
% had been recongnized by the image processor as such
% RowNum - number of rows in the recognized patterns

% Output:
% List of physical coordinates od the calibration dots.

% by Peyruz Gasimov. Aug 2017

function [CalibPatternPhys]=GenPhysDotArray(dg,CalElev,CalCamMask1,CalCamMask2,CalCamMask3,RowNum1,RowNum2,RowNum3)

DotNumBound=6000;

[xg1,yg1]=meshgrid([0:dg:300],linspace(0,-(RowNum1-1)*dg,RowNum1));
[xg2,yg2]=meshgrid([0:dg:300],linspace(0,-(RowNum2-1)*dg,RowNum2));
[xg3,yg3]=meshgrid([0:dg:300],linspace(0,-(RowNum3-1)*dg,RowNum3));

zg1=ones(DotNumBound,1)*CalElev(1);
zg2=ones(DotNumBound,1)*CalElev(2);
zg3=ones(DotNumBound,1)*CalElev(3);

CalibPatternPhys =  [[mat2vec(xg1(CalCamMask1)) mat2vec(yg1(CalCamMask1)) mat2vec(zg1(CalCamMask1))];...
                    [mat2vec(xg2(CalCamMask2)) mat2vec(yg2(CalCamMask2)) mat2vec(zg2(CalCamMask2))];...
                    [mat2vec(xg3(CalCamMask3)) mat2vec(yg3(CalCamMask3)) mat2vec(zg3(CalCamMask3))]];

