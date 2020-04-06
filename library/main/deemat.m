function dee = deemat(dee,e,v)
% !
% ! This subroutine returns the elastic dee matrix for ih=3 (plane strain),
% ! ih=4 (axisymmetry or plane strain elastoplasticity) or ih=6
% ! (three dimensions).
% !
pt5 = 0.5;
one = 1;
two = 2;
ih = size(dee,1);
v1 = one - v;
c = e/((one+v)*(one-two*v));
switch ih
    case 3
        dee(1,1)=v1*c;
        dee(2,2)=v1*c;
        dee(1,2)=v*c;
        dee(2,1)=v*c;
        dee(3,3)=pt5*c*(one-two*v);
    case 4
        dee(1,1)=v1*c;
        dee(2,2)=v1*c;
        dee(4,4)=v1*c;
        dee(3,3)=pt5*c*(one-two*v);
        dee(1,2)=v*c;
        dee(2,1)=v*c;
        dee(1,4)=v*c;
        dee(4,1)=v*c;
        dee(2,4)=v*c;
        dee(4,2)=v*c;
    case 6
        v2 = v/(one-v);
        vv=(one-two*v)/(one-v)*pt5;
        for i = 1:3
            dee(i,i)=one;
        end
        for i = 4:6
            dee(i,i)=vv;
        end
        dee(1,2)=v2;
        dee(2,1)=v2;
        dee(1,3)=v2;
        dee(3,1)=v2;
        dee(2,3)=v2;
        dee(3,2)=v2;
        dee=dee*e/(two*(one+v)*vv);
    otherwise
        fprintf("wrong size for dee matrix")
end
end