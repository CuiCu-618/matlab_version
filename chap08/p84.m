% PROGRAM p84
% !-------------------------------------------------------------------------
% ! Program 8.4 Plane or axisymmetric transient analysis using 4-node
% !             rectangular quadrilaterals. Mesh numbered in x(r)- or y(z)-
% !             direction. Implicit time integration using the "theta"
% !             method.
% !-------------------------------------------------------------------------
%% !-----------------------input and initialisation------------------------
ndim = 2;
nip = 4;
nod = 4;
penalty = 1e20;
element = "quadrilateral";

type_2d = "plane";
dir = "x";
nxe = 5;
nye = 5;
np_types = 1;
[nels,nn] = mesh_size(element,nod,nxe,nye);
neq = nn;

points = zeros(nip,ndim);
weights = zeros(nip,1);
kay = zeros(ndim,ndim);
coord = zeros(nod,ndim);
fun = zeros(nod,1);
jac = zeros(ndim,ndim);
g_coord = zeros(ndim,nn);
der = zeros(ndim,nod);
deriv = zeros(ndim,nod);
mm = zeros(nod,nod);
g_num = zeros(nod,nels);
kc = zeros(nod,nod);
ntn = zeros(nod,nod);           % cross product of shape functions
num = zeros(nod,1);
etype = zeros(nels,1);
kdiag = zeros(neq,1);
loads = zeros(neq,1);
newlo = zeros(neq,1);
x_coord = zeros(nxe+1,1);
y_coord = zeros(nye+1,1);
prop = zeros(ndim,np_types);
gc = zeros(ndim,1);

prop(:,:) = [1,1];
etype(:,:) = 1;
x_coord(:,:) = [0.0 0.2 0.4 0.6 0.8 1.0];
y_coord(:,:) = [0.0 -0.2 -0.4 -0.6 -0.8 -1.0];
dtim = 0.01;
nstep = 150;
theta = 0.5;
npri = 10;
nres = 31;
ntime = 100;
%% ! --------loop the elements to set up global geometry and kdiag --------
for iel = 1:nels
    [coord,num] = geom_rect(element,iel,x_coord,num,y_coord,dir);
    g_num(:,iel) = num;
    g_coord(:,num) = coord';
    kdiag = fkdiag(kdiag,num);
end
for i = 2:neq
    kdiag(i) = kdiag(i) + kdiag(i-1);
end
kv = zeros(kdiag(neq),1);
bp = zeros(kdiag(neq),1);
fprintf(" There are %d equations and the skyline storage is %d \n",...
                neq,kdiag(neq));
[points,weights] = sample(element,points,weights);            
gc(:,:) = 1;
%% !---------------------global conductivity and "mass" matrix assembly----
for iel = 1:nels
    kay(:,:) = 0;
    for i = 1:ndim
        kay(i,i) = prop(i,etype(iel));
    end
    num(:,1) = g_num(:,iel);
    coord(:,:) = g_coord(:,num)';
    kc(:,:) = 0;
    mm(:,:) = 0;
    for i = 1:nip
        der(:,:) = shape_der(der,points,i);
        fun(:,1) = shape_fun(fun,points,i);
        jac = der*coord;
        dete = det(jac);
        deriv = jac\der;
        if type_2d == "axisymmetric"
            gc = fun'*coord;
        end
        kc = kc + deriv'*kay*deriv*dete*weights(i)*gc(1);
        ntn = cross_product(fun,fun);
        mm = mm + ntn*dete*weights(i)*gc(1);
    end
    kv = fsparv(kv,kc,num,kdiag);
    bp = fsparv(bp,mm,num,kdiag);
end
kv = kv*theta*dtim;
bp = bp+kv;
kv = bp-kv/theta;
%% !---------------------specify initial and boundary values---------------
loads(:,1) = [0.0,0.0,0.0,0.0,0.0,0.0,100.0,100.0,100.0,100.0,100.0,0.0,...
              100.0,100.0,100.0,100.0,100.0,0.0,100.0,100.0,100.0,100.0,...
              100.0,0.0,100.0,100.0,100.0,100.0,100.0,0.0,100.0,100.0,...
              100.0,100.0,100.0,0.0];
fixed_freedoms = 11;
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    node(:,1) = [1,2,3,4,5,6,12,18,24,30,36];
    value = zeros(fixed_freedoms,1);
    value(:,1) = zeros(fixed_freedoms,1);
    storbp  = zeros(fixed_freedoms,1); 
    bp(kdiag(node)) = bp(kdiag(node)) + penalty;
    storbp(:,1) = bp(kdiag(node));
end
%% !---------------------factorise equations-------------------------------
bp = sparin(bp,kdiag);
%% !---------------------time stepping loop--------------------------------
fprintf("    Time         Pressure (node %2d)\n",nres)
fprintf("%13.4e %13.4e\n",0,loads(nres))
for j = 1:nstep
    time = j*dtim;
    newlo(:,1) = linmul_sky(kv,loads,kdiag);
    if fixed_freedoms ~= 0
        newlo(node) = storbp.*value;
    end
    newlo = spabac(bp,newlo,kdiag);
    loads = newlo;

    if mod(j,npri) == 0
        fprintf("%13.4e %13.4e\n",time,loads(nres))
    end
end








