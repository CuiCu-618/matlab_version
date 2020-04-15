function a = cross_product(b,c,a)
% !
% ! This subroutine forms the cross product of two vectors, a = b x c
% !
ib = size(b,1);
ic = size(c,1);
for i = 1:ib
    for j = 1:ic
        a(i,j) = b(i)*c(j);
    end
end


end