function loads = solve_band(pb,work,loads)
% !
% ! This subroutine performs Gaussian forward and back-substitution
% ! on the reduced unsymmetric band matrix pb.
% !
pt5 = 0.5;
iwp1 = (size(pb,2)-1)/2+1;
n = size(pb,1);
iq = 2*iwp1-1;
n1 = n-1;
 for iv=1:n1
   i=floor(work(iwp1,iv)+pt5);
   if i~=iv
     s=loads(iv);
     if i == 0
         loads(iv) = 0;
     else
        loads(iv)=loads(i);
        loads(i)=s;
     end
     
   end
   l=iv+iwp1-1;
   if l>n
       l=n;
   end
   iv1=iv+1;
   for i=iv1:l
     loads(i)=loads(i)-work(i-iv,iv)*loads(iv);
   end
 end
 loads(n)=loads(n)/pb(n,1);
 iv=n-1;
 while iv~=0
   s=loads(iv);
   l=iq;
   if iv+l-1>n
       l=n-iv+1;
   end
   for i=2:l
     s=s-pb(iv,i)*loads(iv+i-1);
     loads(iv)=s/pb(iv,1);
   end
 iv=iv-1;
 end
end