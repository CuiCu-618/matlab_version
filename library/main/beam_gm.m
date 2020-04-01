function gm = beam_gm(ell)
% !
% ! This subroutine forms the beam geometric matrix for stability analysis.
% !
pt1 = 0.1;
opt2 = 1.2;
two = 2;
d15 = 15;
d30 = 30;
gm(1,1) = opt2/ell;
gm(1,2) = pt1;
gm(2,1) = pt1;
gm(1,3) = -opt2/ell;
gm(3,1) = -opt2/ell;
gm(1,4) = pt1;
gm(4,1) = pt1;
gm(2,2) = two*ell/d15;
gm(2,3) = -pt1;
gm(3,2) = -pt1;
gm(2,4) = -ell/d30;
gm(4,2) = -ell/d30;
gm(3,3) = opt2/ell;
gm(3,4) = -pt1;
gm(4,3) = -pt1;
gm(4,4) = two*ell/d15;
end