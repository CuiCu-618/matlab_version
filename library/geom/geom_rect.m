function [coord,num] = geom_rect(element,iel,x_coords,num,y_coords,dir)
% !
% ! This subroutine forms the coordinates and connectivity for a
% ! rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
% ! or quadrilateral elements (4, 8 or 9-node) counting in the
% ! x- or y-dir. 
% !
pt5 = 0.5;
two = 2;
d3 = 3;
nxe = size(x_coords,1)-1;
nod = size(num,1);
if element == "triangle"
    nye = (size(y_coords,1)-1)*2;
    if dir == "x" || dir == "r"
        jel = 2*nxe*((iel-1)/(2*nxe));
        ip = (iel-jel+1)/2;
        iq = 2*((iel-1)/(2*nxe)+1)-1+((iel/2)*2)/iel;
    else
        jel = (iel-1)/nye;
        ip = jel+1;
        iq = iel-nye*jel;
    end
end

end