function mm = rod_mm(mm,length)
% !
% ! This subroutine forms the consistent mass matrix of a 1-d "rod" element.
% !
one = 1;
d3 = 3;
d6 = 6;
mm(1,1)=one/d3;
mm(1,2)=one/d6;
mm(2,1)=one/d6;
mm(2,2)=one/d3;
mm=mm*length;
end