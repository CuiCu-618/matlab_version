function ku = formku(ku,km,g)
% !
% ! This subroutine assembles element matrices into symmetrical
% ! global matrix (stored as an upper rectangle).
% !
nforf = size(km,1);
for i=1:nforf
    if g(i)~=0
        for j=1:nforf
            if g(j)~=0
                icd=g(j)-g(i)+1;
                if icd >= 1
                    ku(g(i),icd) = ku(g(i),icd)+km(i,j);
                end
            end
        end
    end
end
end