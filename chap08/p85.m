% PROGRAM p85
% !-------------------------------------------------------------------------
% ! Program 8.5 Plane or axisymmetric consolidation analysis using 4-node
% !             rectangular quadrilaterals. Mesh numbered in x(r)- or y(z)-
% !             direction. Implicit time integration using the "theta"
% !             method. No global matrix assembly. Diagonal
% !             preconditioner conjugate gradient solver
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
cg_tol = 0.0001;
cg_limit = 100;
np_types = 1;
[nels,nn] = mesh_size(element,nod,nxe,nye);
neq = nn;
fprintf(" There are %d equations\n",neq);

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
storka = zeros(nod,nod,nels);
storkb = zeros(nod,nod,nels);
etype = zeros(nels,1);
loads = zeros(neq,1);
x_coord = zeros(nxe+1,1);
y_coord = zeros(nye+1,1);
prop = zeros(ndim,np_types);
gc = zeros(ndim,1);
diag_precon = zeros(neq,1);
u = zeros(neq,1);
d = zeros(neq,1);
p = zeros(neq,1);
x = zeros(neq,1);
r = zeros(neq,1);               % holds ﬁxed rhs terms in pcg solver
xnew = zeros(neq,1);

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
end
[points,weights] = sample(element,points,weights);            
gc(:,:) = 1;
%% !----------element matrix integration, storage and preconditioner------ 
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
    storka(:,:,iel) = mm+kc*theta*dtim;
    storkb(:,:,iel) = mm-kc*(1-theta)*dtim;
    for k = 1:nod
        diag_precon(num(k)) = diag_precon(num(k)) + storkb(k,k,iel);
    end
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
    
    store  = zeros(fixed_freedoms,1); 
    diag_precon(node) = diag_precon(node)+penalty;
    store(:,:) = diag_precon(node);
end
diag_precon = 1./diag_precon;
%% !-----------------------time stepping loop--------------------------------
fprintf("    Time         Pressure (node %2d)    cg iters\n",nres)
fprintf("%13.4e %13.4e\n",0,loads(nres))
for j = 1:nstep
    time = j*dtim;
    u(:,:) = 0;
    for iel = 1:nels
        num = g_num(:,iel);
        kc = storkb(:,:,iel);
        u(num) = u(num) + kc*loads(num);
    end
    r = u;
    if fixed_freedoms ~= 0
        r(node) = store.*value;
    end
    d = diag_precon.*r;
    p = d;
    x(:,:) = 0;
%% !-----------------------pcg equation solution-----------------------------
    for cg_iters = 1:cg_limit
        u(:,:) = 0;             % Q
        for iel = 1:nels
            num = g_num(:,iel);
            kc = storka(:,:,iel);
            u(num) = u(num) + kc*p(num);
        end
        if fixed_freedoms ~= 0
            u(node) = store.*p(node);
        end
        up = r'*d;
        alpha = up/(p'*u);
        xnew = x + p*alpha;
        r = r - u*alpha;
        d = diag_precon.*r;
        beta = r'*d/up;
        p = d + p*beta;
        [cg_converged,x] = checon(xnew,x,cg_tol);
        if cg_converged
            break
        end
    end
    loads = xnew;
    if mod(j,npri) == 0
        fprintf("%13.4e %13.4e      %10d\n",time,loads(nres),cg_iters)
    end
end




