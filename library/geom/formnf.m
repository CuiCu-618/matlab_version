function nf = formnf(nf)
% !
% ! This subroutine forms the nf matrix.
% !
m = 0;
for j = 1:size(nf,2)
    for i = 1:size(nf,1)
        if nf(i,j) ~= 0
            m = m + 1;
            nf(i,j) = m;
        end
    end
end
end