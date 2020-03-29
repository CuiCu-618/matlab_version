function km = rod_km(km,ea,length)
% !
% ! This subroutine forms the stiffness matrix of a 1-d "rod" element.
% !
km(1,1) = 1;
km(2,2) = 1;
km(1,2) = -1;
km(2,1) = -1;
km = km * ea / length;
end