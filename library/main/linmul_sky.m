function loads = linmul_sky(kv,disps,kdiag)
% !
% ! This subroutine forms the product of symmetric matrix stored as
% ! a skyline and a vector.
% !
n = size(disps,1);
loads = zeros(n,1);
for i = 1:n
    x = 0;
    lup = kdiag(i);
    if i == 1
        low = lup;
    end
    if i ~= 1
        low = kdiag(i-1) + 1;
    end
    for j = low:lup
        x = x+kv(j)*disps(i+j-lup);
    end
    loads(i) = x;
    if i == 1
        continue;
    end
    lup = lup - 1;
    for j = low:lup
        k = i+j-lup-1;
        loads(k) = loads(k)+kv(j)*disps(i);
    end
end
end