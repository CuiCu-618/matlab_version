% PROGRAM p89      
% !-------------------------------------------------------------------------
% ! Program 8.9 Plane analysis of the diffusion-convection equation
% !             using 4-node rectangular quadrilaterals. Implicit time
% !             integration using the "theta" method.
% !             Self-adjoint transformation.
% !-------------------------------------------------------------------------
%% !-----------------------input and initialisation------------------------
ndim = 2;
nip = 4;
nod = 4;
penalty = 1e20;
d6 = 6;
d12 = 12;
pt25 = 0.25;
two = 2;
element = "quadrilateral";

nxe = 1;
nye = 40;
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
ans_ = zeros(neq,1);
x_coords = zeros(nxe+1,1);
y_coords = zeros(nye+1,1);
prop = zeros(ndim,np_types);
gc = zeros(ndim,1);

prop(:,:) = [1,1];
etype(:,:) = 1;

x_coords(:,:) = [0,1.4];
y_coords(:,:) = linspace(56,0,41);
dtim = 300;
nstep = 20;
theta = 1;
npri = 1;
nres = 82;
ntime = 5;
ux = 0;
uy = 0.0135;
%% ! --------loop the elements to set up global geometry and kdiag --------
for iel = 1:nels
    [coord,num] = geom_rect(element,iel,x_coords,num,y_coords,"x");
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
        kc = kc + deriv'*kay*deriv*dete*weights(i);
        ntn = cross_product(fun,fun);
        mm = mm + ntn*dete*weights(i);
    end
    kc = kc+mm*(ux*ux/kay(1,1)+uy*uy/kay(2,2))*pt25;
    mm = mm/(theta*dtim);
%% !---------------------derivative boundary conditions--------------------
    if iel == 1
        D = x_coords(2)-x_coords(1);
        kc(2,2) = kc(2,2)+uy*D/d6;
        kc(2,3) = kc(2,3)+uy*D/d12;
        kc(3,2) = kc(3,2)+uy*D/d12;
        kc(3,3) = kc(3,3)+uy*D/d6;
    elseif iel == nels
        D = x_coords(2)-x_coords(1);
        kc(1,1) = kc(1,1)+uy*D/d6;
        kc(1,4) = kc(1,4)+uy*D/d12;
        kc(4,1) = kc(4,1)+uy*D/d12;
        kc(4,4) = kc(4,4)+uy*D/d6;        
    end
    kv = fsparv(kv,kc,num,kdiag);
    bp = fsparv(bp,mm,num,kdiag);
end
f1 = uy*D/(two*theta);
f2 = f1;
bp = bp+kv;
kv = bp-kv/theta;
%% !---------------------factorise equations-------------------------------
bp = sparin(bp,kdiag);
%% !---------------------time stepping loop--------------------------------
fprintf("    Time         Concentration (node %2d)\n",nres)
fprintf("%13.4e %13.4e\n",0,loads(nres))
for j = 1:nstep
    time = j*dtim;
    ans_(:,1) = linmul_sky(kv,loads,kdiag);
    ans_(neq) = ans_(neq) + f1;
    ans_(neq-1) = ans_(neq-1) + f1;
    ans_ = spabac(bp,ans_,kdiag);
    loads = ans_;

    if mod(j,npri) == 0
        fprintf("%13.4e %13.4e\n",time,loads(nres)*exp(-ux*g_coord(1,nres)/two/kay(1,1))*...
                    exp(-uy*g_coord(2,nres)/two/kay(2,2)));
    end
end



