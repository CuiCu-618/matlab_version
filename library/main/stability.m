function [iters,kv,gv,tol,eval,evec] = stability(kv,gv,tol,limit,kdiag)
% !
% ! This subroutine computes the smallest eigenvalue in a beam
% ! stability analysis.
% !
neq = max(max(kdiag));
x0 = zeros(neq+1,1);
x1 = zeros(neq+1,1);
kv = sparin(kv,kdiag);

x0(1) = 1;

MAX_IT = 1000;
for iters = 1:MAX_IT
    x1 = linmul_sky(gv,x0,kdiag);
    x1 = spabac(kv,x1,kdiag);
    big = max(x1(1:end-1));
    if abs(min(x1(1:end-1))) > big
        big = min(x1(1:end-1));
    end
    x1 = x1/big;
    converged = (max(abs(x1(1:end-1)-x0(1:end-1)))/max(abs(x1(1:end-1))) < tol);
    x0 = x1;
    if converged || iters == limit
        break;
    end
end
x1(1:end-1) = x1(1:end-1)/sqrt(sum(x1(1:end-1).^2));
evec = x1;
eval = 1/big;
end