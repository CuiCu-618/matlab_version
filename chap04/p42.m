% PROGRAM p42
% !-------------------------------------------------------------------------
% ! Program 4.2 Analysis of elastic pin-jointed frames using 2-node rod
% !             elements in 2- or 3-dimensions
% !-------------------------------------------------------------------------
%% ---------------------------initialisation-------------------------------
fixed_freedoms = 1;
loaded_nodes = 1;
ndim = 3;                       % number of dimensions
nels = 4;
nod = 2;
nprops = 1;
np_types = 1;
nr = 4;
penalty = 1e20;

nodof = ndim;
ndof = nod*nodof;
nn = 5;

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

prop(:) = 5e5;
etype(:) = 1;
% g_coord(:,:) = [0,4,4,8,8,12;3,0,3,3,0,0];
g_coord(:,:) = [0,1.25,3.5,4,2;0,3,2,1,1.5;0,0,0,0,3];
% g_num(:,:) = [1,1,3,3,3,2,2,5,4,5;2,3,4,5,2,4,5,4,6,6];
g_num(:,:) = [1,2,3,4;5,5,5,5];
nf(:,:) = 1;
k = [1,2,3,4];
% nf(:,k) = [0,1;0,0];
nf(:,k) = [0,0,0,0;0,0,0,0;0,0,0,0];
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
for i = 2:neq
    kdiag(i) = kdiag(i)+kdiag(i-1);
end
kv = zeros(kdiag(neq),1);
fprintf(" There are %d equations and the skyline storage is %d \n",...
                neq,kdiag(neq))
%% !---------------------global stiffness matrix assembly------------------
for iel = 1:nels
    num = g_num(:,iel);
    coord = g_coord(:,num)';
    km = pin_jointed(prop(1,etype(iel)),coord);
    g = g_g(:,iel);
    kv = fsparv(kv,km,g,kdiag);
end
%% !----------------------- loads and/or displacements---------------------
k = 5;
% loads(nf(:,k)) = [0,-10]';
loads(nf(:,k)) = [20,-20,30]';
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    node(:) = 5;
    no = zeros(fixed_freedoms,1);
    sense = zeros(fixed_freedoms,1);
    sense(:) = 2;
    value = zeros(fixed_freedoms,1);
    value(:) = -0.005;
    for i = 1:fixed_freedoms
        no(i) = nf(sense(i),node(i));
    end
    kv(kdiag(no)) = kv(kdiag(no)) + penalty;
    loads(no) = kv(kdiag(no)) * value;
end
%% !----------------------equation solution ------------------------------- 
kv = sparin(kv,kdiag);
loads = spabac(kv,loads,kdiag);
fprintf("  Node       Displacement(s) \n")
for k = 1:nn
    if nf(1,k) == 0
        nf(1,k) = 4;
        loads(nf(1,k)) = 0;
    end
    if nf(2,k) == 0
        nf(2,k) = 4;
        loads(nf(2,k)) = 0;
    end
    if nf(3,k) == 0
        nf(3,k) = 4;
        loads(nf(2,k)) = 0;
    end
    fprintf("   %d   %13.4e  %13.4e\n",k,loads(nf(1,k)),loads(nf(2,k)))
end
%% !----------------------retrieve element end actions---------------------
fprintf(" Element Actions \n")
for iel = 1:nels
    num = g_num(:,iel);
    coord = g_coord(:,num)';
    g = g_g(:,iel);
    g(g==0) = 4;
    eld = loads(g);
    km = pin_jointed(prop(1,etype(iel)),coord);
    action = km*eld;
    axial = glob_to_axial(action,coord);
    fprintf(" %d  %13.4e  %13.4e  %13.4e  %13.4e \n",...
            iel,action(1),action(2),action(3),action(4))
    fprintf("     Axial force = %11.4e\n",axial)
end







