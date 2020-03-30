function km = pin_jointed(ea,coord)
% !
% ! This subroutine forms the stiffness matrix of a
% ! general rod element (1-, 2- or 3-d).
% !
ndim = size(coord,2);
switch ndim
    case 1
        ell = coord(2,1)-coord(1,1);
        km(1,1) = 1;
        km(1,2) = -1;
        km(2,1) = -1;
        km(2,2) = 1;
    case 2
        x1 = coord(1,1);
        y1 = coord(1,2);
        x2 = coord(2,1);
        y2 = coord(2,2);
        ell = sqrt((y2-y1)^2+(x2-x1)^2);
        cs = (x2-x1)/ell;
        sn = (y2-y1)/ell;
        a = cs*cs;
        b = sn*sn;
        c = cs*sn;
        km(1,1) = a;
        km(3,3) = a;
        km(1,3) = -a;
        km(3,1) = -a;
        km(2,2) = b;
        km(4,4) = b;
        km(2,4) = -b;
        km(4,2) = -b;
        km(1,2) = c;
        km(2,1) = c;
        km(3,4) = c;
        km(4,3) = c;
        km(1,4) = -c;
        km(4,1) = -c;
        km(2,3) = -c;
        km(3,2) = -c;
    case 3
        x1 = coord(1,1);
        y1 = coord(1,2);
        z1 = coord(1,3);
        x2 = coord(2,1);
        y2 = coord(2,2);
        z2 = coord(2,3);
        xl = x2-x1;
        yl = y2-y1;
        zl = z2-z1;
        ell = sqrt(xl*xl+yl*yl+zl*zl);
        xl = xl/ell;
        yl = yl/ell;
        zl = zl/ell;
        a = xl*xl;
        b = yl*yl;
        c = zl*zl;
        d = xl*yl;
        e = yl*zl;
        f = zl*xl;
        km(1,1) = a;
        km(4,4) = a;
        km(2,2) = b;
        km(5,5) = b;
        km(3,3) = c;
        km(6,6) = c;
        km(1,2) = d;
        km(2,1) = d;
        km(4,5) = d;
        km(5,4) = d;
        km(2,3) = e;
        km(3,2) = e;
        km(5,6) = e;
        km(6,5) = e;
        km(1,3) = f;
        km(3,1) = f;
        km(4,6) = f;
        km(6,4) = f;
        km(1,4) = -a;
        km(4,1) = -a;
        km(2,5) = -b;
        km(5,2) = -b;
        km(3,6) = -c;
        km(6,3) = -c;
        km(1,5) = -d;
        km(5,1) = -d;
        km(2,4) = -d;
        km(4,2) = -d;
        km(2,6) = -e;
        km(6,2) = -e;
        km(3,5) = -e;
        km(5,3) = -e;
        km(1,6) = -f;
        km(6,1) = -f;
        km(3,4) = -f;
        km(4,3) = -f;
end
km = km*ea/ell;     
end