% Rotate a swarm of particles into a horizontal position

function [rotated]=swarmRotate(swarm)

% Fit a line through the swarm (XZ view)
c = polyfit(swarm(1,:),swarm(3,:),1);

% Calculate the angle of the rotation
angle=atan(c(1));

% Rotation Matrix
Ry=[cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];

% Rotate
rotated=Ry*swarm;

end
