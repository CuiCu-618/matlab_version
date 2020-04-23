function [d,e,ifail] = bisect(d,e,acheps,ifail)
% !
% ! This subroutine finds the eigenvalues of a tridiagonal matrix,
% ! given with its diagonal elements in the array d(n) and
% ! its subdiagonal elements in the last n - 1 stores of the
% ! array e(n), using ql transformations. The eigenvalues are
% ! overwritten on the diagonal elements in the array d in
% ! ascending order. The subroutine will fail if any one
% ! eigenvalue takes more than 30 iterations.
% !
pt5 = 0.5;
two = 2;
one = 1;
n = size(d,1);
if n ~= 1
    for i = 2:n
        e(i-1) = e(i);
    end
end
e(n) = 0;
b = 0;
f = 0;
for l = 1:n
    j = 0;
    h = acheps*(abs(d(l))+abs(e(l)));
    if b < h
    % ! look for small sub diagonal element
        for m = l:n
            if abs(e(m)) <= b
                break
            end
        end
    end
    if m ~= l
        for j = 1:30
            % ! form shift
            g = d(l);
            h = d(l+1) - g;
            if abs(h)<abs(e(l))
                p=h*pt5/e(l);
                r=sqrt(p*p+one);
                h=p+r;
                if p < 0
                    h=p-r;
                end
                d(l)=e(l)/h;
            else
                p=two*e(l)/h;
                r=sqrt(p*p+one);
                d(l)=e(l)*p/(one+r);
            end
            h=g-d(l);
            i1=l+1;
            if i1<=n
             for i=i1:n
               d(i)=d(i)-h;
             end
            end
            f=f+h;
            % ! ql transformation
           p=d(m);
           c=one;
           s=0;
           m1=m-1;
            for ii=l:m1
                 i=m1-ii+l;
                 g=c*e(i);
                 h=c*p;
         if abs(p)>=abs(e(i))
           c= e(i)/p;
           r=sqrt(c*c+one);
           e(i+1)=s*p*r;
           s=c/r;
           c=one/r;
         else
           c=p/e(i);
           r=sqrt(c*c+one);
           e(i+1)=s*e(i)*r;
           s=one/r;
           c=c/r;
         end
         p=c*d(i)-s*g;
         d(i+1)=h+s*(c*g+s*d(i));
            end
       e(l)=s*p;
       d(l)=c*p;
       if abs(e(l))<=b 
           break
       end
        end
        if j == 30
            ifail = 1;
        end
    end
   p=d(l)+f;
   %! order eigenvalue
   aux=0;
   if l~=1
     for ii=2:l
       i=l-ii+2;
       if p>=d(i-1)
         aux=1;
         break
       end
       d(i)=d(i-1);
     end
   end
   if aux==0
     i=1;
   end
   d(i) = p;
   ifail=0;
end
end