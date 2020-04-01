% PROGRAM p44
% !-------------------------------------------------------------------------
% ! Program 4.4 Analysis of elastic rigid-jointed frames using 2-node
% !             beam/rod elements in 2- or 3-dimensions.
% !
%% ---------------------------initialisation-------------------------------
fixed_freedoms = 0;
loaded_nodes = 4;
ndim = 3;
nels = 3;
nod = 2;
nn = 4;
nprops = 4;
np_types = 1;
nr = 2;

if ndim == 2
    nodof = 3;
end
if ndim == 3
    nodof = 6;
end
ndof = nod * nodof;

nf = zeros(nodof,nn);
km = zeros(ndof,ndof);
coord = zeros(nod,ndim);
g_coord = zeros(ndim,nn);
eld = zeros(ndof,1);
action = zeros(ndof,1);
g_num = zeros(nod,nels);
num = zeros(nod,1);
g = zeros(ndof,1);
gamma = zeros(nels,1);          % rotation of element about local axis
g_g = zeros(ndof,nels);
prop = zeros(nprops,np_types);
etype = zeros(nels,1);

% prop(:,:) = [5e9,1e9;6e4,2e4];
prop(:,:) = [4e6,1e6,0.3e6,0.3e6]';
% etype(:) = [1,1,1,2,2,2];
etype(:) = 1;
% g_coord(:,:) = [0,6,6,12,12,14;0,0,-4,0,-5,0];
% g_num(:,:) = [1,2,4,3,3,5;2,4,6,2,4,4];
g_coord(:,:) = [0,5,5,5;5,5,5,0;5,5,0,0];
g_num(:,:) = [1,3,4;2,2,3];
nf(:,:) = 1;
% k = [1,3,5];
k = [1,4];
% nf(:,k) = [0,0,0;0,0,0;1,0,0];
nf(:,k) = zeros(6,2);
nf = formnf(nf);
neq = max(max(nf));
kdiag = zeros(neq,1);
loads = zeros(neq+1,1);
gamma(:,:) = [0,0,90]';
%% !-------------loop the elements to find global arrays sizes-------------
for iel = 1:nels
    num(:,1) = g_num(:,iel);
    g(:) = num_to_g(num,nf);
    g_g(:,iel) = g;
    kdiag = fkdiag(kdiag,g);
end
for i = 2:neq
    kdiag(i) = kdiag(i) + kdiag(i-1);
end
fprintf(" There are %d equations and the skyline storage is %d \n",...
                neq,kdiag(neq))
%% !---------------------global stiffness matrix assembly------------------
kv = zeros(kdiag(neq),1); 
for iel = 1:nels
    num(:,1) = g_num(:,iel);
    coord(:,:) = g_coord(:,num)';
    km = rigid_gointed(prop,gamma,etype,iel,coord);
    g = g_g(:,iel);
    kv = fsparv(kv,km,g,kdiag);
end
%% !----------------------- loads and/or displacements---------------------
% k = [1,2,4,6];
k = 2;
nf(nf == 0) = neq+1;
% loads(nf(:,k)) = [0,0,0,0;-60,-180,-140,-20;-60,-80,133.33,6.67];
loads(nf(:,k)) = [0,-100,0,0,0,0]';
loads(neq+1) = 0;
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    no = zeros(fixed_freedoms,1);
    sense = zeros(fixed_freedoms,1);
    value = zeros(fixed_freedoms,1);
    for i = 1:fixed_freedoms
        no(i) = nf(sense(i),node(i));
    end
    kv(kdiag(no)) = kv(kdiag(no)) + penalty;
    loads(no) = kv(kdiag(no)) .* value;
end
%% !----------------------equation solution ------------------------------- 
kv = sparin(kv,kdiag);
loads = spabac(kv,loads,kdiag);
fprintf("  Node   Displacements and Rotation(s)\n")
for k = 1:nn
    if nf(:,k) == 0
        nf(:,k) = neq+1;
        loads(nf(:,k)) = 0;
    end
    fprintf("   %d   %13.4e  %13.4e  %13.4e %13.4e  %13.4e  %13.4e\n",...
        k,loads(nf(1,k)),loads(nf(2,k)),loads(nf(3,k)),...
          loads(nf(4,k)),loads(nf(5,k)),loads(nf(6,k)))
end
%% !----------------------retrieve element end actions---------------------
fprintf(" Element Actions \n")
for iel = 1:nels
    num(:,1) = g_num(:,iel);
    coord(:,:) = g_coord(:,num)';
    g = g_g(:,iel);
    g(g==0) = neq+1;
    eld = loads(g);
    km = rigid_gointed(prop,gamma,etype,iel,coord);
    action = km*eld;
    fprintf(" %d  %13.4e  %13.4e  %13.4e  %13.4e %13.4e  %13.4e \n",...
           iel,action(1),action(2),action(3),action(4),action(5),action(6))
end











