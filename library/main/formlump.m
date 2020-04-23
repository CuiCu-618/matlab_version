function diag = formlump(diag,emm,g)
% !
% ! This subroutine forms the lumped global mass matrix as a vector diag.
% !
ndof = size(emm,1);
for i = 1:ndof
    diag(g(i)) = diag(g(i)) + emm(i,i);
end
end