function pb = formtb(pb,km,g)
% !
% ! This subroutine assembles an unsymmetrical band matrix pb from
% ! element constituent matrices km.
% !
idof = size(km,1);
iw = (size(pb,2)-1)/2;
for i = 1:idof
    if g(i) ~= 0
        for j = 1:idof
            if g(j) ~= 0
                icd = g(j)-g(i)+iw+1;
                pb(g(i),icd) = pb(g(i),icd)+km(i,j);
            end
        end
    end     
end
end