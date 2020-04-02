% PROGRAM p46              
% !-------------------------------------------------------------------------
% ! Program 4.6 Stability (buckling) analysis of elastic beams using 2-node
% !             beam elements (elastic foundation optional).
% !-------------------------------------------------------------------------
%% ---------------------------initialisation-------------------------------
limit = 100;
ndof = 4;
nels = 4;
nod = 2;
nodof = 2;
nprops = 2;
np_types = 1;
nr = 2;
tol = 1e-5;

nn = nels+1;
nf = zeros(nodof,nn);
ell = zeros(nels,1);
num = zeros(nod,1);
g = zeros(ndof,1);
g_g = zeros(ndof,nels);
etype = zeros(nels,1);
prop = zeros(nprops,np_types);
km = zeros(ndof,ndof);
gm = zeros(ndof,ndof);          % element geometric matrix
mm = zeros(ndof,ndof);

prop(:) = [1,800];
etype(:) = 1;
ell(:) = [0.25,0.25,0.25,0.25];
nf(:,:) = 1;
k = [1,5];
nf(:,k) = [0,0;1,1];
nf = formnf(nf);
neq = max(max(nf));
kdiag = zeros(neq,neq);
evec = zeros(neq+1,0);          % eigenvector (mode shape)
%% !-------------loop the elements to find global arrays sizes-------------
for iel = 1:nels
    num(:,1) = [iel,iel+1];
    g(:) = num_to_g(num,nf);
    g_g(:,iel) = g;
    kdiag = fkdiag(kdiag,g);
end
for i = 2:neq
    kdiag(i) = kdiag(i) + kdiag(i-1);
end
kv = zeros(kdiag(neq),1);
gv = zeros(kdiag(neq),1);       % global geometric matrix
fprintf(" There are %d equations and the skyline storage is %d \n",...
                neq,kdiag(neq))
%% !---------------------global stiffness matrix assembly------------------
for iel = 1:nels
    km = beam_km(prop(1,etype(iel)),ell(iel));
    if nprops > 1
        mm = beam_mm(prop(2,etype(iel)),ell(iel));
    end
    gm = beam_gm(ell(iel));
    g = g_g(:,iel);
    kv = fsparv(kv,km+mm,g,kdiag);
    gv = fsparv(gv,gm,g,kdiag);
end
%% !----------------------- loads and/or displacements---------------------
[iters,kv,gv,tol,eval,evec] = stability(kv,gv,tol,limit,kdiag);
fprintf(" The buckling load = %13.4e \n", eval)
evec(end) = 0;
fprintf(" The buckling mode\n")
fprintf("  Node  Translation    Rotation \n")
for k = 1:nn
    if nf(1,k) == 0
        nf(1,k) = neq+1;
        evec(nf(1,k)) = 0;
    end
    if nf(2,k) == 0
        nf(2,k) = neq+1;
        evec(nf(2,k)) = 0;
    end
    fprintf("   %d   %13.4e  %13.4e\n",k,evec(nf(1,k)),evec(nf(2,k)))
end
fprintf(" Converged in  %d  iterations\n", iters)













