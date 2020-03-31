% PROGRAM p43
% !-------------------------------------------------------------------------
% ! Program 4.3 Analysis of elastic beams using 2-node beam elements
% !             (elastic foundation optional).
% !-------------------------------------------------------------------------
%% ---------------------------initialisation-------------------------------
fixed_freedoms = 4;
loaded_nodes = 4;
ndof = 4;
nod = 2;
nodof = 2;
nels = 4;
nprops = 1;
np_types = 2;
nr = 0;
penalty = 1e20;

nn = nels + 1;
g = zeros(ndof,1);
num = zeros(nod,1);
nf = zeros(nodof,nn);
etype = zeros(nels,1);
ell = zeros(nels,1);
eld = zeros(ndof,1);
km = zeros(ndof,ndof);
mm = zeros(ndof,ndof);          % element ‘mass’ matrix
action = zeros(ndof,1);
g_g = zeros(ndof,nels);
prop = zeros(nprops,np_types);

prop(:) = [4e4,2e4];
etype(:) = [1,1,2,2];
ell(:) = [2.5,2.5,3,2];
nf(:,:) = 1;
nr = 0;
nf = formnf(nf);
neq = max(max(nf));
kdiag = zeros(neq,1);
loads = zeros(neq+1,1);
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
kv = zeros(kdiag(neq),1);         
for iel = 1:nels
    km = beam_km(prop(1,etype(iel)),ell(iel));
    if nprops > 1
        mm = beam_mm(prop(2,etype(iel)),ell(iel));
    end
    g = g_g(:,iel);
    kv = fsparv(kv,km+mm,g,kdiag);
end
%% !----------------------- loads and/or displacements---------------------
k = [2,3,4,5];
loads(nf(:,k)) = [-20,-6,-8.8,-1.2;0,-3,2.2,0.5333];
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    no = zeros(fixed_freedoms,1);
    sense = zeros(fixed_freedoms,1);
    value = zeros(fixed_freedoms,1);
    node(:) = [1,1,3,4];
    sense(:) = [1,2,1,1];
    value(:) = [0,-0.001,-0.005,0];
    for i = 1:fixed_freedoms
        no(i) = nf(sense(i),node(i));
    end
    kv(kdiag(no)) = kv(kdiag(no)) + penalty;
    loads(no) = kv(kdiag(no)) .* value;
end
%% !----------------------equation solution ------------------------------- 
kv = sparin(kv,kdiag);
loads = spabac(kv,loads,kdiag);
fprintf("  Node  Translation    Rotation \n")
for k = 1:nn
    if nf(:,k) == 0
        nf(:,k) = nn;
        loads(nf(:,k)) = 0;
    end
    fprintf("   %d   %13.4e  %13.4e\n",k,loads(nf(1,k)),loads(nf(2,k)))
end
%% !----------------------retrieve element end actions---------------------
fprintf(" Element Force         Moment           Force       Moment\n")
for iel = 1:nels
    km = beam_km(prop(1,etype(iel)),ell(iel));
    if nprops > 1
        mm = beam_mm(prop(2,etype(iel)),ell(iel));
    end
    g = g_g(:,iel);
    g(g==0) = nn;
    eld = loads(g);
    action = (km+mm)*eld;
    fprintf(" %d  %13.4e  %13.4e  %13.4e  %13.4e \n",...
           iel,action(1),action(2),action(3),action(4))
end







