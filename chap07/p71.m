% PROGRAM p71
% !-------------------------------------------------------------------------
% ! Program 7.1 One dimensional analysis of steady seepage using
% !             2-node "rod" elements.
% !-------------------------------------------------------------------------
%% !---------------------input and initialisation--------------------------
nod = 2;
nprops = 1;
penalty = 1e20;

nels = 6;
nn = 5;
np_types = 3;
neq = nn;

ell = zeros(nels,1);
num = zeros(nod,1);
prop = zeros(nprops,np_types);
etypes = zeros(nels,1);
kp = zeros(nod,nod);
g_num = zeros(nod,nels);
kdiag = zeros(neq,1);
loads = zeros(neq,1);
disps = zeros(neq,1);         % net nodal inﬂow/outﬂow

prop(:,:) = [3,2,1]';
etypes(:,:) = [1,3,2,2,1,3];
ell(:,:) = [1,1,1.5,1.5,1,1];
g_num(:,:) = [1,2,1,3,3,4;2,3,3,5,4,5];
%% !---------------------loop the elements to find global arrays sizes-----
for iel = 1:nels
    num = g_num(:,iel);
    kdiag = fkdiag(kdiag,num);
end
for i = 2:neq
    kdiag(i) = kdiag(i)+kdiag(i-1);
end
fprintf(" There are %d equations and the skyline storage is %d \n",...
                neq,kdiag(neq));
kv = zeros(kdiag(neq),1);
kvh = zeros(kdiag(neq),1);      % copy of kv
%% !---------------------global conductivity matrix assembly---------------
for iel = 1:nels
    kp = rod_km(kp,prop(1,etypes(iel)),ell(iel));
    num = g_num(:,iel);
    kv = fsparv(kv,kp,num,kdiag);
end
kvh = kv;
%% !---------------------specify boundary values---------------------------
loaded_nodes = 2;
k = [3,5];
loads(k) = [75,10];
fixed_freedoms = 1;
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    node(:,:) = [1];
    value = zeros(fixed_freedoms,1);
    value(:,:) = [0];
    kv(kdiag(node)) = kv(kdiag(node)) + penalty;
    loads(node) = kv(kdiag(node)) .* value;
end
%% !---------------------equation solution---------------------------------
kv = sparin(kv,kdiag);
loads = spabac(kv,loads,kdiag);
%% !---------------------retrieve nodal net flow rates---------------------
disps = linmul_sky(kvh,loads,kdiag);
fprintf("  Node Total Head  Flow rate \n")
for k = 1:nn
    fprintf("%2d %13.4e %13.4e\n",k,loads(k),disps(k))
end
fprintf("     Inflow      Outflow\n")
fprintf("%13.4e %13.4e\n",sum(disps(disps>0)),sum(disps(disps<0)))






