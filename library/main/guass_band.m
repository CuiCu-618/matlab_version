function [pb,work] = guass_band(pb,work)
% !
% ! This subroutine performs gaussian reduction of an unsymmetric
% ! banded matrix pb. Array work used as working space.
% !
n = size(pb,1);
iwp1 = (size(pb,2)-1)/2+1;
zero = 0;
iq = 2*iwp1-1;
iqp = iwp1;
iwp11 = iwp1-1;
small = 1e-10;
for i=1:iwp11
    for j=1:iq
        if j>=iwp1+i
            pb(i,j)=zero;
            pb(n-i+1,j)=zero;
        else
            pb(i,j)=pb(i,j+iwp1-i);
        end
    end
end
for k=1:n
   l=k+iwp1-1;
   if l>n 
       l=n;
   end
   ip=0;
   s=small;
   for i=k:l
     if abs(pb(i,1))<=s 
         continue
     end
     s=abs(pb(i,1));
     ip=i;
   end
   if ip==0
     fprintf("singular")
     break
   end
   if k==n
       break
   end
   work(iwp1,k)=ip;
   iqp=iqp-1;
   j=iwp1+ip-k;
   if iqp<j
       iqp=j;
   end
   if j~=iwp1
     for j=1:iqp
       s=pb(k,j);
       pb(k,j)=pb(ip,j);
       pb(ip,j)=s;
     end
   end
   k1=k+1;
   for i=k1:l
     s=pb(i,1)/pb(k,1);
     for j=2:iq
       if j>iqp
         pb(i,j-1)=pb(i,j);
       else
         pb(i,j-1)=pb(i,j)-s*pb(k,j);
       end
     end
     pb(i,iq)=0;
     work(i-k,k)=s;
   end
end
end