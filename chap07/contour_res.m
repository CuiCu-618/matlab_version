[X1,Y1] = meshgrid(x_coords,y_coords);
Z1 = zeros(nxe,nye);
k = 1;
for i = 1:nye+1
    for j = 1:nxe+1
        Z1(i,j) = loads(k);
        k = k+1;
    end
end

[X2,Y2] = meshgrid(x_coords+6,y_coords);
k = 1;
for i = 1:nye+1
    for j = nxe+1:-1:1
        Z2(i,j) = loads(k);
        k = k+1;
    end
end

X = [X1,X2];
Y = [Y1,Y2];
Z = [Z1,Z2];

contour(X,Y,Z,10)
[DX,DY] = gradient(Z,.2,.2);
hold on
quiver(X,Y,DX,-DY)
