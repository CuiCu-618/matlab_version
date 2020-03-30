function axial = glob_to_axial(glob,coord)
% !
% ! This subroutine transforms the global end reactions
% ! into an axial force for rod elements (2- or 3-d).
% !
ndim = size(coord,2);
add = 0;
for i = 1:ndim
    add = add + (coord(2,i)-coord(1,i))^2;
end
ell = sqrt(add);
axial = 0;
for i = 1:ndim
    axial = axial + (coord(2,i)-coord(1,i))/ell*glob(ndim+i);
end
end