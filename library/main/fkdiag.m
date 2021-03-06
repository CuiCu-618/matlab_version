function kdiag = fkdiag(kdiag,g)
% !
% ! This subroutine computes the skyline profile.
% !
idof = size(g,1);
for i = 1:idof
    iwp1 = 1;
    if g(i) ~= 0
        for j = 1:idof
            if g(j) ~= 0
                im = g(i)-g(j)+1;
                if im > iwp1
                    iwp1 = im;
                end
            end
        end
        k = g(i);
        if iwp1 > kdiag(k)
            kdiag(k) = iwp1;
        end
    end
end
end