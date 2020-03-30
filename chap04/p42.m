% PROGRAM p42
% !-------------------------------------------------------------------------
% ! Program 4.2 Analysis of elastic pin-jointed frames using 2-node rod
% !             elements in 2- or 3-dimensions
% !-------------------------------------------------------------------------
%% ---------------------------initialisation-------------------------------
fixed_freedom = 0;
loaded_nodes = 1;
ndim = 2;                       % number of dimensions
nels = 10;
nod = 2;
nnprops = 1;
np_types = 1;
nr = 2;
penalty = 1e20;

nodof = ndim;
ndof = nod*nodof;
nn = 6;

nf = zeros(nodof,nn);
km = zeros(ndof,ndof);
coord = zeros(nod,ndim);        % element nodal coordinates
g_coord = zeros(ndim,nn);       % nodal coordinates for all elements
eld = zeros(ndof,1);
action = zeros(ndof,1);
g_num = zeros(nod,nels);        % global element node numbers matrix
num = zeros(nod,1);
g = zeros(ndof,1);
g_g = zeros(ndof,nels);
etype = zeros(nels,1);
prop = zeros(nprops,np_types);

prop(:) = 2e5;
etype(:) = 1;
g_coord(:,:) = [0,4,4,8,8,12;3,0,3,3,0,0];
g_num(:,:) = [1,1,3,3,3,2,2,5,4,5;2,3,4,5,2,4,5,4,6,6];
nf(:,:) = 1;
k = [1,2];
nf(:,k) = [0,1;0,0];
nf = formnf(nf);
neq = max(max(nf));
kdiag = zeros(neq,1);
loads = zeros(neq+1,1);
%% !--------------------loop the elements to find global array sizes-------
for iel = 1:nels
    num = g_num(:,iel);
    g(:) = num_to_g(num,nf);
    g_g(:,iel) = g;
    kdiag = fkdiag(kdiag,g);
end













