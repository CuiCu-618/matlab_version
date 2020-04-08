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
%         jel = 2*nxe*((iel-1)/(2*nxe));
%         ip = (iel-jel+1)/2;
%         iq = 2*((iel-1)/(2*nxe)+1)-1+((iel/2)*2)/iel;
        ip = ceil(iel/(2*nxe));
        iq = iel-(ip-1)*2*nxe;
    else
%         jel = (iel-1)/nye;
%         ip = jel+1;
%         iq = iel-nye*jel;
        ip = ceil(iel/(nye));
        iq = iel-(ip-1)*nye;
    end
    switch nod
        case 3
            if mod(iq,2) ~= 0
                if dir == "x" || dir == "r"
                    num(1) = (nxe+1)*(iq-1)/2+ip;
                    num(2) = num(1)+1;              
                    num(3) = (nxe+1)*(iq+1)/2+ip;
                else
                    num(1) = (ip-1)*(nye+2)/2+(iq+1)/2;
                    num(2) = num(1)+(nye+2)/2;
                    num(3) = num(1)+1;
                end
                coord(1,1) = x_coords(ip);
                coord(1,2) = y_coords((iq+1)/2);
                coord(2,1) = x_coords(ip+1);   
                coord(2,2) = y_coords((iq+1)/2);
                coord(3,1) = x_coords(ip);   
                coord(3,2) = y_coords((iq+3)/2);
            else
                if dir == "x" || dir == "r"
                    num(1) = (nxe+1)*iq/2+ip+1;     
                    num(2) = num(1)-1;               
                    num(3) = (nxe+1)*(iq-2)/2+ip+1;
                else
                    num(1) = ip*(nye+2)/2+(iq+2)/2;
                    num(2) = (ip-1)*(nye+2)/2+(iq+1)/2+1;
                    num(3) = num(1)-1;
                end
                coord(1,1) = x_coords(ip+1);
                coord(1,2) = y_coords((iq+2)/2);
                coord(2,1) = x_coords(ip);   
                coord(2,2) = y_coords((iq+2)/2);
                coord(3,1) = x_coords(ip+1);
                coord(3,2) = y_coords(iq/2);
            end
        case 6
            if mod(iq,2) ~= 0
                if dir == "x" || dir == "r"
                    num(1) = (iq-1)*(2*nxe+1)+2*ip-1;
                    num(2) = num(1)+1; 
                    num(3) = num(1)+2; 
                    num(4) = (iq-1)*(2*nxe+1)+2*nxe+2*ip+1;
                    num(5) = (iq+1)*(2*nxe+1)+2*ip-1;
                    num(6) = num(4)-1;
                else
                    num(1) = 2*(nye+1)*(ip-1)+iq;
                    num(2) = 2*(nye+1)*(ip-1)+nye+1+iq;
                    num(3) = 2*(nye+1)*ip+iq;
                    num(4) = num(2)+1;
                    num(5) = num(1)+2; 
                    num(6) = num(1)+1;
                end
                coord(1,1) = x_coords(ip);
                coord(1,2) = y_coords((iq+1)/2);
                coord(3,1) = x_coords(ip+1);   
                coord(3,2) = y_coords((iq+1)/2);
                coord(5,1) = x_coords(ip);   
                coord(5,2) = y_coords((iq+3)/2);
            else
                if dir == "x" || dir == "r"
                    num(1) = iq*(2*nxe+1)+2*ip+1;
                    num(2) = num(1)-1; 
                    num(3) = num(1)-2; 
                    num(4) = (iq-2)*(2*nxe+1)+2*nxe+2*ip+1;
                    num(5) = (iq-2)*(2*nxe+1)+2*ip+1;
                    num(6) = num(4)+1; 
                else
                    num(1) = 2*(nye+1)*ip+iq+1; 
                    num(2) = 2*(nye+1)*(ip-1)+nye+iq+2;
                    num(3) = 2*(nye+1)*(ip-1)+iq+1;
                    num(4) = num(2)-1; 
                    num(5) = num(1)-2;
                    num(6) = num(1)-1;
                end
                coord(1,1) = x_coords(ip+1);
                coord(1,2) = y_coords((iq+2)/2);
                coord(3,1) = x_coords(ip);   
                coord(3,2) = y_coords((iq+2)/2);
                coord(5,1) = x_coords(ip+1); 
                coord(5,2) = y_coords(iq/2);
            end
            coord(2,:) = pt5*(coord(1,:)+coord(3,:));
            coord(4,:) = pt5*(coord(3,:)+coord(5,:));
            coord(6,:) = pt5*(coord(5,:)+coord(1,:));
        case 15
            if mod(iq,2) ~= 0
                if dir == "x" || dir == "r"
                    fac1 = 4*(4*nxe+1)*(iq-1)/2;
                    num(1) = fac1+4*ip-3;
                    num(2) = num(1)+1;
                    num(3) = num(1)+2;
                    num(4) = num(1)+3;
                    num(5) = num(1)+4;
                    num(6) = fac1+ 4*nxe+1+4*ip;
                    num(7) = fac1+ 8*nxe+1+4*ip;
                    num(8) = fac1+12*nxe+1+4*ip;
                    num(9) = fac1+16*nxe+1+4*ip;
                    num(10) = num(8)-1;
                    num(11) = num(7)-2;
                    num(12) = num(6)-3;
                    num(13) = num(12)+1;
                    num(14) = num(12)+2;
                    num(15) = num(11)+1;
                else
                    fac1 = 4*(2*nye+1)*(ip-1)+2*iq-1; 
                    num(1) = fac1;
                    num(2) = fac1+2*nye+1;
                    num(3) = fac1+4*nye+2; 
                    num(4) = fac1+6*nye+3; 
                    num(5) = fac1+8*nye+4;
                    num(6) = fac1+6*nye+4; 
                    num(7) = fac1+4*nye+4; 
                    num(8) = fac1+2*nye+4;
                    num(9) = fac1+4; 
                    num(10) = fac1+3; 
                    num(11) = fac1+2;
                    num(12) = fac1+1;
                    num(13) = fac1+2*nye+2;
                    num(14) = fac1+4*nye+3;
                    num(15) = fac1+2*nye+3;  
                end
                coord(1,1) = x_coords(ip);
                coord(1,2) = y_coords((iq+1)/2);
                coord(5,1) = x_coords(ip+1);   
                coord(5,2) = y_coords((iq+1)/2);
                coord(9,1) = x_coords(ip);   
                coord(9,2) = y_coords((iq+3)/2);
            else
                if dir == "x" || dir == "r"
                    fac1 = 4*(4*nxe+1)*(iq-2)/2;
                    num(1) = fac1+16*nxe+5+4*ip;
                    num(2) = num(1)-1;
                    num(3) = num(1)-2;
                    num(4) = num(1)-3;
                    num(5) = num(1)-4;
                    num(6) = fac1+12*nxe+1+4*ip;
                    num(7) = fac1+8*nxe+1+4*ip;
                    num(8) = fac1+4*nxe+1+4*ip;
                    num(9) = fac1+4*ip+1;
                    num(10) = num(8)+1;
                    num(11) = num(7)+2;
                    num(12) = num(6)+3;
                    num(13) = num(12)-1;
                    num(14) = num(12)-2;
                    num(15) = num(11)-1;
                else
                    fac1 = 4*(2*nye+1)*(ip-1)+2*iq+8*nye+5; 
                    num(1) = fac1; 
                    num(2) = fac1-2*nye-1;
                    num(3) = fac1-4*nye-2; 
                    num(4) = fac1-6*nye-3; 
                    num(5) = fac1-8*nye-4;
                    num(6) = fac1-6*nye-4;  
                    num(7) = fac1-4*nye-4;
                    num(8) = fac1-2*nye-4;
                    num(9) = fac1-4;
                    num(10) = fac1-3; 
                    num(11) = fac1-2;
                    num(12) = fac1-1;
                    num(13) = fac1-2*nye-2; 
                    num(14) = fac1-4*nye-3;
                    num(15) = fac1-2*nye-3; 
                end
                coord(1,1) = x_coords(ip+1);
                coord(1,2) = y_coords((iq+2)/2);
                coord(5,1) = x_coords(ip);   
                coord(5,2) = y_coords((iq+2)/2);
                coord(9,1) = x_coords(ip+1);
                coord(9,2) = y_coords(iq/2);
            end
            coord(3,:) = pt5*(coord(1,:)+coord(5,:));
            coord(7,:) = pt5*(coord(5,:)+coord(9,:));
            coord(11,:) = pt5*(coord(9,:)+coord(1,:));
            coord(2,:) = pt5*(coord(1,:)+coord(3,:));
            coord(4,:) = pt5*(coord(3,:)+coord(5,:));
            coord(6,:) = pt5*(coord(5,:)+coord(7,:));
            coord(8,:) = pt5*(coord(7,:)+coord(9,:));
            coord(10,:) = pt5*(coord(9,:)+coord(11,:));
            coord(12,:) = pt5*(coord(11,:)+coord(1,:));
            coord(15,:) = pt5*(coord(7,:)+coord(11,:));
            coord(14,:) = pt5*(coord(3,:)+coord(7,:));
            coord(13,:) = pt5*(coord(2,:)+coord(15,:));
        otherwise
            fprintf("Wrong number of nodes for triangular element, nod = 3 , 6 or 15\n")
    end
else
    nye = size(y_coords,1)-1;
    if dir == "x" || dir == "r"
        ip = ceil(iel/nxe);
        iq = iel-(ip-1)*nxe;
    else
        ip = ceil(iel/nye);
        iq = iel-(ip-1)*nye;
    end
    switch nod
        case 4
            if dir == "x" || dir == "r"
                num(1) = iq*(nxe+1)+ip;		        		
                num(2) = (iq-1)*(nxe+1)+ip;				
                num(3) = num(2)+1;					
                num(4) = num(1)+1;
            else
                num(1) = (ip-1)*(nye+1)+iq+1;
                num(2) = num(1)-1;
                num(3) = ip*(nye+1)+iq;
                num(4) = num(3)+1;
            end
            coord(1:2,1) = x_coords(ip);
            coord(3:4,1) = x_coords(ip+1);
            coord(1,2) = y_coords(iq+1);
            coord(2:3,2) = y_coords(iq);
            coord(4,2) = coord(1,2);
        case 8
            if dir == "x" || dir == "r"
                num(1) = iq*(3*nxe+2)+2*ip-1;                
                num(2) = iq*(3*nxe+2)+ip-nxe-1;		  
                num(3) = (iq-1)*(3*nxe+2)+2*ip-1;		   
                num(4) = num(3)+1;
                num(5) = num(4)+1;
                num(6) = num(2)+1;
                num(7) = num(1)+2;
                num(8) = num(1)+1;
            else
                num(1) = (ip-1)*(3*nye+2)+2*iq+1;
                num(2) = num(1)-1;
                num(3) = num(1)-2;
                num(4) = (ip-1)*(3*nye+2)+2*nye+iq+1;
                num(5) = ip*(3*nye+2)+2*iq-1;
                num(6) = num(5)+1;
                num(7) = num(5)+2;
                num(8) = num(4)+1;
            end
            coord(1:3,1) = x_coords(ip);
            coord(5:7,1) = x_coords(ip+1);
            coord(4,1) = pt5*(coord(3,1)+coord(5,1));
            coord(8,1) = pt5*(coord(7,1)+coord(1,1));
            coord(1,2) = y_coords(iq+1);
            coord(7:8,2) = y_coords(iq+1);
            coord(3:5,2) = y_coords(iq);
            coord(2,2) = pt5*(coord(1,2)+coord(3,2));
            coord(6,2) = pt5*(coord(5,2)+coord(7,2));
        otherwise
            fprintf("Wrong number of nodes for quadrilateral element, nod = 4 or 8 \n")

    end
end

end