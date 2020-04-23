function [omega,gamma,s] = form_s(gg,ell,kappa,omega,gamma,s)
% !
% ! This subroutine forms the s vector in bicgstab(l)
% !
one = 1;
zero = 0;
hh = zeros(ell-1,ell-1);
gamma0 = zeros(ell+1,1);
p = zeros(ell-1,1);
q = zeros(ell-1,1);
gamma1 = zeros(ell+1,1);
hh = -gg(2:ell,2:ell);
p = hh\gg(2:ell,1);
q = hh\gg(2:ell,ell+1);
gamma0(1) = one;
gamma0(ell+1) = zero;
gamma0(2:ell) = p;
gamma1(1) = zero;
gamma1(ell+1) = one;
gamma1(2:ell) = q;
ngamma0 = gamma0'*(gg*gamma0);
ngamma1 = gamma1'*(gg*gamma1);
omega = gamma0'*(gg*gamma1);
cosine = abs(omega)/sqrt(abs(ngamma0*ngamma1));
omega = omega/ngamma1;
if cosine < kappa
    omega = kappa/cosine*omega;
end
gamma = gamma0 - omega*gamma1;
s(1:ell) = gamma(2:ell+1);
s(ell+1) = zero;
end