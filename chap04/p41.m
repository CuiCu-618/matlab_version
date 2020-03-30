% PROGRAM p41
% !------------------------------------------------------------------------
% ! Program 4.1 One dimensional analysis of axially loaded elastic rods
% !             using 2-node rod elements.
% !------------------------------------------------------------------------
%% ---------------------------initialisation-------------------------------
ndof = 2;                       % number of freedoms per element
nod = 2;                        % number of nodes per element
nodof = 1;                      % number of freedoms per node
nprops = 1;                     % number of material properties
np_types = 1;                   % number of different property types
nels = 4;                       % number of elements
nn = nels + 1;                  % number of nodes
nr = 1;                         % number of restrained nodes
loaded_nodes = 5;               % number of loaded nodes
fixed_freedoms = 0;             % number of ﬁxed displacements

g = zeros(ndof,1);              % element steering vector
num = zeros(nod,1);             % element node numbers vector
nf = zeros(nodof,nn);           % nodal freedom matrix
etype = zeros(nels,1);          % element property type vector
ell = zeros(nels,1);            % element lengths vector
eld = zeros(ndof,1);            % element displacement vector
km = zeros(ndof,ndof);          % element stiffness matrix
action = zeros(ndof,1);         % element nodal action vector
g_g = zeros(ndof,nels);         % global element steering matrix
prop = zeros(nprops,np_types);  % number of material properties

prop(:) = 100000.0;
etype(:) = 1;
nf(:,:) = 1;
k = 5;                          % simple counter
nf(:,k) = 0;
ell(:) = 0.25;
penalty = 1e20;
nf = formnf(nf);
neq = max(max(nf));             % number of degrees of freedom in the mesh
kdiag = zeros(neq,1);           % diagonal term location vector
loads = zeros(neq+1,1);         % global load (displacement) vector
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
fprintf(" There are %d equations and the skyline storage is %d \n",...
                neq,kdiag(neq))
%% !---------------------global stiffness matrix assembly------------------
kv = zeros(kdiag(neq),1);         % global stiffness matrix
for iel = 1:nels
    % km : element stiffness matrix
    km = rod_km(km,prop(1,etype(iel)),ell(iel));    
    g = g_g(:,iel);
    kv = fsparv(kv,km,g,kdiag);
end
%% !----------------------- loads and/or displacements---------------------
loads(:) = [-0.625 -1.25 -1.25 -1.25 -0.625];
if fixed_freedoms ~= 0
    % fixed nodes vector
    node = zeros(fixed_freedoms,1);
    % ﬁxed freedom numbers vector
    no = zeros(fixed_freedoms,1);
    value = zeros(fixed_freedoms,1);
    for i = 1:fixed_freedoms
        no(i) = nf(1,node(i));
    end
    kv(kdiag(no)) = kv(kdiag(no)) + penalty;
    loads(no) = kv(kdiag(no)) * value;
end
%% !----------------------equation solution ------------------------------- 
kv = sparin(kv,kdiag);
loads = spabac(kv,loads,kdiag);
fprintf("  Node       Disp \n")
for k = 1:nn
    if nf(:,k) == 0
        nf(:,k) = nn;
        loads(nf(:,k)) = 0;
    end
    fprintf("   %d   %13.4e  \n",k,loads(nf(:,k)))
end
%% !----------------------retrieve element end actions---------------------
fprintf(" Element Actions \n")
for iel = 1:nels
    km = rod_km(km,prop(1,etype(iel)),ell(iel)); 
    g = g_g(:,iel);
    g(g==0) = 5;
    eld = loads(g);
    action = km*eld;
    fprintf(" %d  %13.4e  %13.4e \n",iel,action(1),action(2))
end


















