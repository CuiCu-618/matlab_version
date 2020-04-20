% PROGRAM p87
% !-------------------------------------------------------------------------
% ! Program 8.7 Plane or axisymmetric analysis of the consolidation equation
% !             using 4-node rectangular quadrilaterals. Mesh numbered in
% !             x(r)- or y(z)- direction. "theta" method using an
% !             "element-by-element" (ebe) product algorithm.
% !-------------------------------------------------------------------------
%% !-----------------------input and initialisation------------------------
ndim = 2;
nip = 4;
nod = 4;
penalty = 1e20;
pt5 = 0.5;
one = 1;
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
globma = zeros(nn,1);
store_kc = zeros(nod,nod,nels);
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
[points,weights] = sample(element,points,weights);            
gc(:,:) = 1;
%%  !-------create and store element and global lumped mass matrices--------
for iel = 1:nels
    [coord,num] = geom_rect(element,iel,x_coord,num,y_coord,dir);
    kay(:,:) = 0;
    for i = 1:ndim
        kay(i,i) = prop(i,etype(iel));
    end
    g_num(:,iel) = num;
    g_coord(:,num) = coord';
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
    store_kc(:,:,iel) = kc;
    for i = 1:nod
        globma(num(i)) = globma(num(i)) + sum(mm(i,:));
    end
end
%% !---------------------recover element a and b matrices------------------
for iel = 1:nels
    num(:,:) = g_num(:,iel);
    kc = -store_kc(:,:,iel)*(one-theta)*dtim*pt5;
    mm = store_kc(:,:,iel)*theta*dtim*pt5;
    for i = 1:nod
        mm(i,i) = mm(i,i) + globma(num(i));
        kc(i,i) = kc(i,i) + globma(num(i));
    end
    mm = mm\kc;
    store_kc(:,:,iel) = mm;
end
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
end
%% !---------------------time stepping loop--------------------------------
fprintf("    Time         Pressure (node %2d)\n",nres)
fprintf("%13.4e %13.4e\n",0,loads(nres))
for j = 1:nstep
    time = j*dtim;
%% !---------------------first pass (1 to nels)----------------------------
    for iel = 1:nels
        num(:,:) = g_num(:,iel);
        mm = store_kc(:,:,iel);
        loads(num) = mm*loads(num);
        loads(node) = value;
    end
%% !-----------------------second pass (nels to 1)---------------------------
    for iel = nels:-1:1
        num(:,:) = g_num(:,iel);
        mm = store_kc(:,:,iel);
        loads(num) = mm*loads(num);
        loads(node) = value;
    end
    if mod(j,npri) == 0
        fprintf("%13.4e %13.4e\n",time,loads(nres))
    end   

end


















