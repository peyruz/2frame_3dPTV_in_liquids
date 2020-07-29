     
function R=vecAngle2RotMatrix(v,ang)

R=  [1+(1-cos(ang))*(v(1)*v(1)-1)             -v(3)*sin(ang)+(1-cos(ang))*v(1)*v(2)	v(2)*sin(ang)+(1-cos(ang))*v(1)*v(3);
	 v(3)*sin(ang)+(1-cos(ang))*v(1)*v(2)	 1+(1-cos(ang))*(v(2)*v(2)-1)            -v(1)*sin(ang)+(1-cos(ang))*v(2)*v(3);
 	-v(2)*sin(ang)+(1-cos(ang))*v(1)*v(3)	 v(1)*sin(ang)+(1-cos(ang))*v(2)*v(3)	1+(1-cos(ang))*(v(3)*v(3)-1)];

end