function [converged,oldlds] = checon(loads,oldlds,tol)
% !
% ! This subroutine sets converged to .FALSE. if relative change in loads
% ! and oldlds is greater than tol and updates oldlds.
% !
converged = true;
converged = max(abs(loads-oldlds))/max(abs(loads)) <= tol;
oldlds = loads;
end