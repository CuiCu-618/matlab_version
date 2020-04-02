function beamdis(loads,nf,ratmax,interp,nels,ell,argv)
% !
% ! This subroutine produces a PostScript output file "*.dis" displaying
% ! the deformed 1-D finite element mesh.
% !
scale = 72;
pt5 = 0.5;
opt5 = 1.5;
fpt5 = 5.5;
d8 = 8;
ept5 = 8.5;
d11 = 11;
one = 1;
two = 2;
thr = 3;

nn = nels+1;
xcoord = zeros(nn,1);
xmin = 0;
xmax = 0;
xnow = 0;
ymin = 0;
ymax = 0;

xcoord(1) = xnow;

for i = 2:nn
    xnow = xnow + ell(i-1);
    xcoord(i) = xnow;
end
xmax = xcoord(nn);

for i = 1:nn
    if loads(nf(1,i)) < ymin
        ymin = loads(nf(1,i));
    end
    if loads(nf(1,i)) > ymax
        ymax = loads(nf(1,i));
    end
end

width = xmax-xmin;
height = ymax-ymin;
dmax = ratmax*width;
if height > width
    dmax = ratmax*height;
end

vmax = 0;
for i = 1:nn
    if abs(loads(nf(1,i))) > vmax
        vmax = abs(loads(nf(1,i)));
    end
end
dismag = dmax/vmax;

ymin = 0;
ymax = 0;
for i = 1:nn
    if dismag*loads(nf(1,i)) < ymin
        ymin = dismag*loads(nf(1,i));
    end
    if dismag*loads(nf(1,i)) > ymax
        ymax = dismag*loads(nf(1,i));
    end
    
    if loads(nf(1,i)) < ymin
        ymin = loads(nf(1,i));
    end
    if loads(nf(1,i)) > ymax
        ymax = loads(nf(1,i));
    end
        
end

width = xmax-xmin;
height = ymax-ymin;
% !
% !        allow 1.5" margin minimum on each side of figure
% !
% !        portrait mode 
% !
if height >= d11/ept5*width
% !
% !               height governs the scale
% !
    sxy = scale*d8/height;
    xo = scale*pt5*(ept5-d8*width/height);
    yo = scale*opt5;
else
% !
% !               width governs the scale
% !
    sxy = scale*fpt5/width;
    xo = scale*opt5;
    yo = scale*pt5*(d11-fpt5*width/height);
end
% !
% !
% !         start PostScript output
% !
file1 = strcat(argv,'.dis');
file_1 = fopen(file1,'w');
fprintf(file_1,'%%PS-Adobe-1.0\n');
fprintf(file_1,'%%DocumentFonts: none\n');
fprintf(file_1,'%%Pages: 1\n');
fprintf(file_1,'%%EndComments\n');
fprintf(file_1,'/m {moveto} def\n');
fprintf(file_1,'/l {lineto} def\n');
fprintf(file_1,'/s {stroke} def\n');
fprintf(file_1,'/c {closepath} def\n');
fprintf(file_1,'%%EndProlog\n');
fprintf(file_1,'%%Page: 0 1\n');
fprintf(file_1,'gsave\n');
% !
% !                 draw the deformed mesh
% !
fprintf(file_1,'%9.2f %9.2f  translate\n',xo,yo);
fprintf(file_1,'%9.2f  setlinewidth\n',0.5);

wnew = loads(nf(1,1));

for i = 1:nels
    for j = 1:interp
        wold = wnew;
        localx = j*ell(i)/interp;
        globalx = localx+xcoord(i);
        ll = ell(i);
        na = (one/(ll^3))*(ll^3-thr*ll*localx^2+two*localx^3);
        nb = (one/(ll^2))*(localx*ll^2-two*ll*localx^2+localx^3);
        nc = (one/(ll^3))*(thr*ll*localx^2-two*localx^3);
        nd = (one/(ll^2))*(localx^3-ll*localx^2);
        wnew = na*loads(nf(1,i))+nb*loads(nf(2,i))+...
            nc*loads(nf(1,1+i))+nd*loads(nf(2,i+1));
        
        x = sxy*((globalx-ell(i)/interp)-xmin);
        y = sxy*(dismag*wold-ymin);
        fprintf(file_1,'%9.2f %9.2f  m\n',x,y);
        x = sxy*(globalx-xmin);
        y = sxy*(dismag*wnew-ymin);
        fprintf(file_1,'%9.2f %9.2f  l\n',x,y);
        fprintf(file_1,'c s\n');
    end
end

fprintf(file_1,'grestore\n');
fprintf(file_1,'showpage\n');







end