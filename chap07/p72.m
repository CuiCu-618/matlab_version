% PROGRAM p72             
% !-------------------------------------------------------------------------
% ! Program 7.2 Plane or axisymmetric analysis of steady seepage using
% !             4-node rectangular quadrilaterals. Mesh numbered
% !             in x(r)- or y(z)- direction.
% !-------------------------------------------------------------------------
%% !---------------------input and initialisation--------------------------
ndim = 2;
nip = 4;
nod = 4;
penalty = 1e20;
element = "quadrilateral";

type_2d = "plane";
dir = "x";
nxe = 6;
nye = 6;
np_types = 1;
[nels,nn] = mesh_size(element,nod,nxe,nye);
neq = nn;

points = zeros(nip,ndim);
g_coord = zeros(ndim,nn);
coord = zeros(nod,ndim);
jac = zeros(ndim,ndim);
weights = zeros(nip,1);
der = zeros(ndim,nod);
deriv = zeros(ndim,nod);
kp = zeros(nod,nod);
num = zeros(nod,1);
g_num = zeros(nod,nels);
kay = zeros(ndim,ndim);         % permeability matrix
etypes = zeros(nels,1);
x_coords = zeros(nxe+1,1);
y_coords = zeros(nye+1,1);
prop = zeros(ndim,np_types);
gc = zeros(ndim,1);
fun = zeros(nod,1);
kdiag = zeros(neq,1);
loads = zeros(neq,1);
disps = zeros(neq,1);

prop(:,:) = [1,1];
etypes(:,:) = 1;
x_coords(:,:) = [0,1,2,3,4,5,6];
y_coords(:,:) = -[0,1,2,3,4,5,6];
%% !---------------------loop the elements to find global arrays sizes-----
for iel = 1:nels
    [coord,num] = geom_rect(element,iel,x_coords,num,y_coords,dir);
    g_num(:,iel) = num;
    g_coord(:,num) = coord';
    kdiag = fkdiag(kdiag,num);
end
for i = 2:neq
    kdiag(i) = kdiag(i)+kdiag(i-1);
end
fprintf(" There are %d equations and the skyline storage is %d \n",...
                neq,kdiag(neq));
kv = zeros(kdiag(neq),1);
kvh = zeros(kdiag(neq),1);
[points,weights] = sample(element,points,weights);
gc(:) = 1;
%% !---------------------global conductivity matrix assembly---------------
for iel = 1:nels
    kay(:,:) = 0;
    for i = 1:ndim
        kay(i,i) = prop(i,etypes(iel));
    end
    num(:,1) = g_num(:,iel);
    coord = g_coord(:,num)';
    kp(:,:) = 0;
    for i = 1:nip
        der(:,:) = shape_der(der,points,i);
        fun(:,1) = shape_fun(fun,points,i);
        jac = der*coord;
        dete = det(jac);
        deriv = jac\der;
        if type_2d == "axisymmetric"
            gc = fun*coord;
        end
        kp = kp + deriv'*kay*deriv*dete*weights(i)*gc(1);
    end
    kv = fsparv(kv,kp,num,kdiag);
end
kvh = kv;
%% !---------------------specify boundary values---------------------------
loaded_nodes = 0;
fixed_freedoms = 11;
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    node(:,:) = [1,2,3,4,5,6,7,28,35,42,49];
    value = zeros(fixed_freedoms,1);
    value(:,:) = [0,0,0,0,0,0,0,50,50,50,50];
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











