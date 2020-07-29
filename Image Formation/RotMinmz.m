function ObjFun=RotMinmz(q,rd)

q=quaternion(q);
ec=[0 0 1]';

ObjFun=(acos(dot(RotateVector(q,ec),ec))-rd)^2+(norm(q)-1)^2;

end
