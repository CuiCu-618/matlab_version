function [coord,num] = geom_freesurf(iel,nxe,fixed_seep,fixed_down,down,width,angs,surf,coord,num)
% !
% ! This subroutine forms the coordinates and steering vector
% ! for 4-node quads numbering in the x-direction
% ! (Laplace's equation, variable mesh, 1-freedom per node).
% !
pt5 = 0.5;
small = 0.0001;
d180 = 180;
angr = angs*pi/d180;
iq = floor((iel-1)/nxe)+1;
ip = iel-(iq-1)*nxe;
num(1) = iq*(nxe+1)+ip;
num(2) = (iq-1)*(nxe+1)+ip;
num(3) = num(2)+1;
num(4) = num(1)+1;
if iq < fixed_seep+1
    b1=(surf(ip)-down)/(fixed_seep+1);
    b2=(surf(ip+1)-down)/(fixed_seep+1);
    coord(1,2)=down+(fixed_seep+1-iq)*b1;
    coord(2,2)=down+(fixed_seep+2-iq)*b1;
    coord(3,2)=down+(fixed_seep+2-iq)*b2; 
    coord(4,2)=down+(fixed_seep+1-iq)*b2;
else
    b1=(fixed_down+fixed_seep-iq)/(fixed_down-1); 
    b2=(fixed_down+fixed_seep-iq+1)/(fixed_down-1); 
    coord(1,2)=down*b1;
    coord(2,2)=down*b2;
    coord(3,2)=coord(2,2);
    coord(4,2)=coord(1,2);
end
if abs(angr(ip)-pi*pt5) < small
    fac1 = 0;
else
    tan1 = tan(angr(ip));
    fac1 = 1/tan1;
end
if abs(angr(ip+1)-pi*pt5) < small
    fac2 = 0;
else
    tan2 = tan(angr(ip+1));
    fac2 = 1/tan2;
end
coord(1,1)=width(ip)+coord(1,2)*fac1; 
coord(2,1)=width(ip)+coord(2,2)*fac1;
coord(3,1)=width(ip+1)+coord(3,2)*fac2; 
coord(4,1)=width(ip+1)+coord(4,2)*fac2;

end