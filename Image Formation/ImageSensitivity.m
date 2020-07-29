clear
clc
% Calculate image sensitivity
system_settings4;

xmin=10;   % Set up the dimensions of the calibration pattern
xmax=200;
ymin=-110;
ymax=-15;
dotNx=70;   % Number of dots in x/y
dotNy=50;

depth=-6;   % Depth
dz=1;

[s0c]=GenerateRecPattern(xmin,xmax,dotNx,ymin,ymax,dotNy,depth);
[s1c]=GenerateRecPattern(xmin,xmax,dotNx,ymin,ymax,dotNy,depth-dz);

[s0xg, s0yg]=meshgrid(linspace(10,200,70), linspace(-15, -110,50));

% Capture images
RInt=0; % in air
[xi01air, ~]=Image_Formation_ver4(c1, c2, A1,A2, R1, R2, es13,es23, s0c,oi1,oi2,RInt,n);
[xi11air, ~]=Image_Formation_ver4(c1, c2, A1,A2, R1, R2, es13,es23, s1c,oi1,oi2,RInt,n);

RInt=1; % in in liquid
n=1.5;
[xi01liq, ~]=Image_Formation_ver4(c1, c2, A1,A2, R1, R2, es13,es23, s0c,oi1,oi2,RInt,n);
[xi11liq, ~]=Image_Formation_ver4(c1, c2, A1,A2, R1, R2, es13,es23, s1c,oi1,oi2,RInt,n);

% diff in xw, yw and zw are known

% Calculate difference in xi, yi and zi for air shots
% X-x sensitivity
xi01airx=reshape(xi01air(1,1:dotNx*dotNy),dotNy,dotNx);
diff_xi01airXx=abs(xi01airx(:,1:end-1)-xi01airx(:,2:end));
diff_xi01airXx=[diff_xi01airXx, NaN(50,1)];

% X-y sensitivity
diff_xi01airXy=abs(xi01airx(1:end-1,:)-xi01airx(2:end,:));
diff_xi01airXy=[diff_xi01airXy; NaN(1,70)];

% Y-y
xi01airy=reshape(xi01air(2,1:dotNx*dotNy),dotNy,dotNx);
diff_xi01airYy=abs(xi01airy(1:end-1,:)-xi01airy(2:end,:));
diff_xi01airYy=[diff_xi01airYy; NaN(1,70)];

% Y-x
diff_xi01airYx=abs(xi01airy(:,1:end-1)-xi01airy(:,2:end));
diff_xi01airYx=[diff_xi01airYx, NaN(50,1)];

% X-z and Y-z
xi11airx=reshape(xi11air(1,1:dotNx*dotNy),dotNy,dotNx);
xi11airy=reshape(xi11air(2,1:dotNx*dotNy),dotNy,dotNx);
diffx_z_air=abs(xi11airx-xi01airx);
diffy_z_air=abs(xi11airy-xi01airy);

% Calculate difference in xi, yi and zi for liquid shots
% X-x sensitivity
xi01liqx=reshape(xi01liq(1,1:dotNx*dotNy),dotNy,dotNx);
diff_xi01liqXx=abs(xi01liqx(:,1:end-1)-xi01liqx(:,2:end));
diff_xi01liqXx=[diff_xi01liqXx, NaN(50,1)];

% X-y sensitivity
diff_xi01liqXy=abs(xi01liqx(1:end-1,:)-xi01liqx(2:end,:));
diff_xi01liqXy=[diff_xi01liqXy; NaN(1,70)];

% Y-y
xi01liqy=reshape(xi01liq(2,1:dotNx*dotNy),dotNy,dotNx);
diff_xi01liqYy=abs(xi01liqy(1:end-1,:)-xi01liqy(2:end,:));
diff_xi01liqYy=[diff_xi01liqYy; NaN(1,70)];

% Y-x
diff_xi01liqYx=abs(xi01liqy(:,1:end-1)-xi01liqy(:,2:end));
diff_xi01liqYx=[diff_xi01liqYx, NaN(50,1)];

% X-z and Y-z
xi11liqx=reshape(xi11liq(1,1:dotNx*dotNy),dotNy,dotNx);
xi11liqy=reshape(xi11liq(2,1:dotNx*dotNy),dotNy,dotNx);
diffx_z_liq=abs(xi11liqx-xi01liqx);
diffy_z_liq=abs(xi11liqy-xi01liqy);

