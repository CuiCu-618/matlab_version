% PROGRAM p73
% !-------------------------------------------------------------------------
% ! Program 7.3 Analysis of plane free-surface flow using 4-node
% !             quadrilaterals. "Analytical" form of element conductivity
% !             matrix.
% !-------------------------------------------------------------------------
%% !---------------------input and initialisation--------------------------
ndim = 2;
nod = 4;
d180 = 180;
penalty = 1e20;

nxe = 8;
nye = 7;
tol = 0.01;
limit = 20;
np_types = 1;

nels = nxe*nye;
nn = (nxe+1)*(nels/nxe+1); 
neq = nn;

g_coord = zeros(ndim,nn);
coord = zeros(nod,ndim);
bottom_width = zeros(nxe+1,1);  % x-coordinates of nodes at base of mesh
top_width = zeros(nxe+1,1);     % x-coordinates of initial nodes at top of mesh
surf = zeros(nxe+1,1);          % holds current total head values of free surface
angs = zeros(nxe+1,1);          % angles made by sloping mesh lines to horizontal
kp = zeros(nod,nod);
num = zeros(nod,1);
g_num = zeros(nod,nels);
prop = zeros(ndim,np_types);
kdiag = zeros(neq,1);
kay = zeros(ndim,ndim);
etypes = zeros(nels,1);
loads = zeros(neq,1);
disps = zeros(neq,1);
oldpot = zeros(neq,1);          % nodal total head values from previous iteration

prop(:,:) = [0.001,0.001];
etypes(:,:) = 1;
bottom_width(:,:) = 0:8;
top_width(:,:) = 0:8;
initial_height = 7;
surf(:,:) = initial_height;
angs(:,:) = atan(surf./(top_width-bottom_width))*d180/pi;
hup = 7;                        % ﬁxed total head on upstream side
fixed_up = 8;                   % number of nodes on upstream side
node_up = [9,18,27,36,45,54,63,72];
hdown = 2;                      % ﬁxed total head on downstream side
fixed_down = 3;
node_down = [46,55,64];
fixed_seep = nels/nxe-fixed_down; % number of nodes on seepage surface
node_seep = zeros(fixed_seep,1);% nodes ﬁxed on downstream seepage surface
for i = 1:fixed_seep
    node_seep(i) = i*(nxe+1)+1;
end
%% !---------------------loop the elements to find global arrays sizes-----
for iel = 1:nels
    [coord,num] = geom_freesurf(iel,nxe,fixed_seep,fixed_down,hdown,bottom_width,angs,surf,coord,num);
    kdiag = fkdiag(kdiag,num);
end
for i = 2:neq
    kdiag(i) = kdiag(i)+kdiag(i-1);
end
fprintf(" There are %d equations and the skyline storage is %d \n",...
                neq,kdiag(neq));
kv = zeros(kdiag(neq),1);
kvh = zeros(kdiag(neq),1);
%% !---------------------global conductivity matrix assembly---------------
for iters = 1:limit
    kv(:,:) = 0;
    for iel = 1:nels
        kay(:,:) = 0;
        for i = 1:ndim
            kay(i,i) = prop(i,etypes(iel));
        end
        [coord,num] = geom_freesurf(iel,nxe,fixed_seep,fixed_down,hdown,bottom_width,angs,surf,coord,num);
        g_num(:,iel) = num;
        g_coord(:,num) = coord';
        kp = seep4(kp,coord,kay);
        kv = fsparv(kv,kp,num,kdiag);
    end
    kvh = kv;

%% !---------------------specify boundary values---------------------------
    loads(:,:) = 0;
    kv(kdiag(node_up))=kv(kdiag(node_up))+penalty;
    loads(node_up)=kv(kdiag(node_up))*hup;
    kv(kdiag(node_down))=kv(kdiag(node_down))+penalty;
    loads(node_down)=kv(kdiag(node_down))*hdown;
    kv(kdiag(node_seep))=kv(kdiag(node_seep))+penalty;
    for i = 1:fixed_seep
        loads(node_seep(i))=kv(kdiag(node_seep(i)))* ...
       (hdown+(surf(1)-hdown)*(fixed_seep+1-i)/(fixed_seep+1));
    end
%% !---------------------equation solution---------------------------------
    kv = sparin(kv,kdiag);
    loads = spabac(kv,loads,kdiag);
    surf(1:nxe) = loads(1:nxe);
%% !---------------------check convergence---------------------------------
    [converged,oldpot] = checon(loads,oldpot,tol);
    if converged
        break
    end
end
disps = linmul_sky(kvh,loads,kdiag);
fprintf("  Node Total Head  Flow rate \n")
for k = 1:nn
    fprintf("%2d %13.4e %13.4e\n",k,loads(k),disps(k))
end
fprintf("     Inflow      Outflow\n")
fprintf("%13.4e %13.4e\n",sum(disps(disps>0)),sum(disps(disps<0)))
fprintf(" Converged in %2d iterations\n",iters)






