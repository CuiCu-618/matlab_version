% PROGRAM p101
% !-------------------------------------------------------------------------
% ! Program 10.1 Eigenvalue analysis of elastic beams using 2-node
% !              beam elements. Lumped mass.
% !-------------------------------------------------------------------------
%% !---------------------input and initialisation--------------------------
nod = 2;
nodof = 2;
ndof = 4;
nprops = 2;
d12 = 12;
pt5 = 0.5;
penalty = 1e20;
etol = 1e-30;

nels = 5;
np_types = 1;
nn = nels+1;

nf = zeros(nodof,nn);
km = zeros(ndof,ndof);
num = zeros(nod,1);
g = zeros(ndof,1);
g_g = zeros(ndof,nels);
prop = zeros(nprops,np_types);
mm = zeros(ndof,ndof);
ell = zeros(nels,1);
etype = zeros(nels,1);

prop(:,:) = [1/12,1];
etype(:,:) = 1;

ell(:) = 0.8*ones(nels,1);
nf(:,:) = 1;
nr = 1;
k = 1;
nf(:,k) = [0;0];
nf = formnf(nf);
neq = max(max(nf));
diag = zeros(neq+1,1);
udiag = zeros(neq+1,1);         % working space vector or eigenvector (in Lanczos)
kdiag = zeros(neq,1);
rrmass = zeros(neq+1,1);        % reciprocal square rooted global lumped mass matrix
%% !---------------------loop the elements to find global array sizes------
nband = 0;
for iel = 1:nels
    num(:,:) = [iel,iel+1];
    g(:) = num_to_g(num,nf);
    g_g(:,iel) = g;
    if nband < bandwidth_(g)
        nband = bandwidth_(g);
    end
    kdiag = fkdiag(kdiag,g);
end
for i = 2:neq
    kdiag(i) = kdiag(i) + kdiag(i-1);
end
fprintf(" There are %d equations\n The half-bandwidth (including diagonal) is %d\n The skyline storage is %d \n",...
            neq,nband+1,kdiag(neq))
%% !---------------------global stiffness and mass matrix assembly---------
ku = zeros(neq,nband+1);        % global stiffness matrices
kv = zeros(kdiag(neq));
kh = zeros(kdiag(neq));         % element stiffness matrices
for iel = 1:nels
    g(:) = g_g(:,iel);
    g(g==0) = neq+1;
    mm(:,:) = 0;
    mm(1,1) = pt5*prop(2,etype(iel))*ell(iel);
    mm(3,3) = mm(1,1);
    mm(2,2) = mm(1,1)*ell(iel)^2/d12;
    mm(4,4) = mm(2,2);
    diag = formlump(diag,mm,g);
    km = beam_km(prop(1,etype(iel)),ell(iel));
    g(g==neq+1) = 0;
    ku = formku(ku,km,g);
end
%% !---------------------reduce to standard eigenvalue problem-------------
rrmass(1:end-1) = 1./sqrt(diag(1:end-1));
for i = 1:neq
    if i <= neq-nband
        k = nband+1;
    else
        k = neq-i+1;
    end
    for j = 1:k
        ku(i,j) = ku(i,j)*rrmass(i)*rrmass(i+j-1);
    end
end
%% !---------------------convert to skyline form---------------------------
kh(1) = ku(1,1);
k = 1;
for i = 2:neq
    idiag = kdiag(i)-kdiag(i-1);
    for j = 1:idiag
        k = k+1;
        kh(k) = ku(i+j-idiag,1-j+idiag);
    end
end
%% !---------------------extract the eigenvalues---------------------------
[ku,diag,udiag] = bandred(ku,diag,udiag);
ifail = 1;
[diag,udiag,ifail] = bisect(diag,udiag,etol,ifail);
fprintf(" The eigenvalues are:\n")
fprintf("%13.4e\n",diag(1:end-1))
%% !---------------------extract the eigenvectors--------------------------
nmodes = 3;






