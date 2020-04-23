function [a,d,e] = bandred(a,d,e)
% !
% ! This subroutine transforms a real symmetric band matrix a,
% ! of order n and band width iw, to tridiagonal form by an appropriate
% ! sequence of Jacobi rotations. During the transformation the
% ! property of the band matrix is maintained. The method yields
% ! a tridiagonal matrix, the diagonal elements of which are in
% ! d(n) and off-diagonal elements in e(n).
% !
one = 1;
two = 2;
n = size(a,1);
iw = size(a,2)-1;
n2 = n-2;
if n2 >= 1
    for k = 1:n2
        maxr = iw;
        if n-k < iw
            maxr = n-k;
        end
        for irr = 2:maxr
            ir = 2+maxr-irr;
            kr = k+ir;
            for j = kr:iw:n
                if j ~= kr
                    if abs(g) < 1e-15
                        break
                    end
                    jm = j-iw;
                    b=-a(jm-1,iw+1)/g;
                    iugl=j-iw;
                else
                    if abs(a(k,ir+1)) < 1e-15
                        break
                    end
                    b=-a(k,ir)/a(k,ir+1);
                    iugl=k;
                end
                s = one/sqrt(one+b*b);
                c = b*s;
                c2 = c*c;
                s2 = s*s;
                cs = c*s;
                u = c2*a(j-1,1)-two*cs*a(j-1,2)+s2*a(j,1);
                u1 = s2*a(j-1,1)+two*cs*a(j-1,2)+c2*a(j,1);
                a(j-1,2) = cs*(a(j-1,1)-a(j,1))+(c2-s2)*a(j-1,2);
                a(j-1,1) = u;
                a(j,1) = u1;
                j2 = j-2;
                for l = iugl:j2
                    jl = j-l;
                    u = c*a(l,jl)-s*a(l,jl+1);
                    a(l,jl+1) = s*a(l,jl)+c*a(l,jl+1);
                    a(l,jl) = u;
                end
                jm = j-iw;
                if j ~= kr
                    a(jm-1,iw+1) = c*a(jm-1,iw+1)-s*g;
                end
                maxl = iw-1;
                if n-j < iw-1
                    maxl = n-j;
                end
                if maxl > 0
                    for l = 1:maxl
                        u = c*a(j-1,l+2)-s*a(j,l+1);
                        a(j,l+1) = s*a(j-1,l+2)+c*a(j,l+1);
                        a(j-1,l+2) = u;
                    end
                end
                if j+iw <= n
                    g = -s*a(j,iw+1);
                    a(j,iw+1) = c*a(j,iw+1);
                end
            end
        end
    end
end
e(1) = 0;
d(1:n) = a(1:n,1);
if 2 <= n
    for i = 2:n
        e(i) = a(i-1,2);
    end
end
end
