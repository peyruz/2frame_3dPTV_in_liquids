% Script for plotting the system of two cameras and their FOV's. Start the
% script after the systems_settings have been initiated.

scatter3([A1(1),c1(1),c2(1),A2(1)],[A1(2),c1(2),c2(2),A2(2)],[A1(3),c1(3),c2(3),A2(3)])
hold on

sc=10   % Scale
lw=3    % LineWidth
quiver3(oi1(1),oi1(2),oi1(3),es11(1),es11(2),es11(3),sc,'LineWidth',lw,'color','blue')
quiver3(oi1(1),oi1(2),oi1(3),es12(1),es12(2),es12(3),sc,'LineWidth',lw,'color','red')
quiver3(oi1(1),oi1(2),oi1(3),es13(1),es13(2),es13(3),sc,'LineWidth',lw,'color','black')

quiver3(oi2(1),oi2(2),oi2(3),es23(1),es23(2),es23(3),sc,'LineWidth',lw,'color','black')
quiver3(oi2(1),oi2(2),oi2(3),es22(1),es22(2),es22(3),sc,'LineWidth',lw,'color','red')
quiver3(oi2(1),oi2(2),oi2(3),es21(1),es21(2),es21(3),sc,'LineWidth',lw,'color','blue')

oi12=oi1+Hs*es12
oi13=oi12+Ws*es11
oi14=oi1+es11*Ws

% scatter3([oi12(1),oi13(1),oi14(1)],[oi12(2),oi13(2),oi14(2)],[oi12(3),oi13(3),oi14(3)])


oi22=oi2+es22*Hs
oi23=oi22+es21*Ws
oi24=oi2+es21*Ws

% scatter3([oi22(1),oi23(1),oi24(1)],[oi22(2),oi23(2),oi24(2)],[oi22(3),oi23(3),oi24(3)])

axis equal
% axis([0 450 -200 0 0 300])

patch([oi1(1),oi12(1),oi13(1),oi14(1)],[oi1(2),oi12(2),oi13(2),oi14(2)],[oi1(3),oi12(3),oi13(3),oi14(3)],'g')
patch([oi2(1),oi22(1),oi23(1),oi24(1)],[oi2(2),oi22(2),oi23(2),oi24(2)],[oi2(3),oi22(3),oi23(3),oi24(3)],'g')
%%
% Plot FOV

%% Solve for the points of intersection of the sensor corner rays and
% xy-plane
l11=-dot(e3,oi1)/dot((c1-oi1),e3)
l12=-dot(e3,oi12)/dot((c1-oi12),e3)
l13=-dot(e3,oi13)/dot((c1-oi13),e3)
l14=-dot(e3,oi14)/dot((c1-oi14),e3)

l24=-dot(e3,oi24)/dot((c2-oi24),e3)
l23=-dot(e3,oi23)/dot((c2-oi23),e3)
l22=-dot(e3,oi22)/dot((c2-oi22),e3)
l21=-dot(e3,oi2)/dot((c2-oi2),e3)

f11=(c1-oi1)*l11+oi1
f12=(c1-oi12)*l12+oi12
f13=(c1-oi13)*l13+oi13
f14=(c1-oi14)*l14+oi14

f21=(c2-oi2)*l21+oi2
f22=(c2-oi22)*l22+oi22
f23=(c2-oi23)*l23+oi23
f24=(c2-oi24)*l24+oi24

%%
plot3([oi2(1),f21(1)],[oi2(2),f21(2)],[oi2(3),f21(3)])
plot3([oi22(1),f22(1)],[oi22(2),f22(2)],[oi22(3),f22(3)])
plot3([oi23(1),f23(1)],[oi23(2),f23(2)],[oi23(3),f23(3)])
plot3([oi24(1),f24(1)],[oi24(2),f24(2)],[oi24(3),f24(3)])

plot3([oi1(1),f11(1)],[oi1(2),f11(2)],[oi1(3),f11(3)])
plot3([oi12(1),f12(1)],[oi12(2),f12(2)],[oi12(3),f12(3)])
plot3([oi13(1),f13(1)],[oi13(2),f13(2)],[oi13(3),f13(3)])
plot3([oi14(1),f14(1)],[oi14(2),f14(2)],[oi14(3),f14(3)])

patch([f21(1),f22(1),f23(1),f24(1)],[f21(2),f22(2),f23(2),f24(2)],[f21(3),f22(3),f23(3),f24(3)],'b', 'FaceAlpha',0.5)
patch([f11(1),f12(1),f13(1),f14(1)],[f11(2),f12(2),f13(2),f14(2)],[f11(3),f12(3),f13(3),f14(3)],'b', 'FaceAlpha',0.5)