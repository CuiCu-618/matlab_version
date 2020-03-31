function km = rigid_gointed(prop,gamma,etype,iel,coord)
% !
% ! This subroutine forms the stiffness matrix of a
% ! general beam/column element (1-, 2- or 3-d).
% !
two = 2;
d4 = 4;
d6 = 6;
d12 = 12;
d180 = 180;
ndim = size(coord,2);
switch ndim
    case 1
        ei = prop(1,etype(iel));
        ell = coord(2,1)-coord(1,1);
        km(1,1) = d12*ei/(ell*ell*ell); 
        km(3,3) = km(1,1);
        km(1,2) = d6*ei/(ell*ell); 
        km(2,1) = km(1,2); 
        km(1,4) = km(1,2);
        km(4,1) = km(1,4); 
        km(1,3) = -km(1,1); 
        km(3,1) = km(1,3); 
        km(3,4) = -km(1,2);
        km(4,3) = km(3,4); 
        km(2,3) = km(3,4); 
        km(3,2) = km(2,3);
        km(2,2) = d4*ei/ell;
        km(4,4) = km(2,2); 
        km(2,4) = two*ei/ell; 
        km(4,2) = km(2,4);
    case 2
        ea = prop(1,etype(iel));
        ei = prop(2,etype(iel));
        x1 = coord(1,1);
        y1 = coord(1,2);
        x2 = coord(2,1);
        y2 = coord(2,2);
        ell = sqrt((y2-y1)^2+(x2-x1)^2);
        c = (x2-x1)/ell;
        s = (y2-y1)/ell;
        e1 = ea/ell;
        e2 = d12*ei/(ell*ell*ell);
        e3 = ei/ell;
        e4 = d6*ei/(ell*ell);
        km(1,1) = c*c*e1+s*s*e2;
        km(4,4) = km(1,1);
        km(1,2) = s*c*(e1-e2);
        km(2,1) = km(1,2);
        km(4,5) = km(1,2);
        km(5,4) = km(4,5);
        km(1,3) = -s*e4;
        km(3,1) = km(1,3);
        km(1,6) = km(1,3);
        km(6,1) = km(1,6);
        km(3,4) = s*e4;
        km(4,3) = km(3,4);
        km(4,6) = km(3,4);
        km(6,4) = km(4,6);
        km(1,4) = -km(1,1); 
        km(4,1) = km(1,4);
        km(1,5) = s*c*(-e1+e2);
        km(5,1) = km(1,5);
        km(2,4) = km(1,5);
        km(4,2) = km(2,4);
        km(2,2) = s*s*e1+c*c*e2;
        km(5,5) = km(2,2);
        km(2,5) = -km(2,2);
        km(5,2) = km(2,5);
        km(2,3) = c*e4;
        km(3,2) = km(2,3);
        km(2,6) = km(2,3);
        km(6,2) = km(2,6);
        km(3,3) = d4*e3;
        km(6,6) = km(3,3);
        km(3,5) = -c*e4;
        km(5,3) = km(3,5);
        km(5,6) = km(3,5);
        km(6,5) = km(5,6);
        km(3,6) = two*e3;
        km(6,3) = km(3,6);
    case 3
end

end