% PROGRAM p92
% !--------------------------------------------------------------------------
% ! Program 9.2 Analysis of the plane steady state Navier-Stokes equation
% !             using 8-node rectangular quadrilaterals for velocities
% !             coupled to 4-node rectangular quadrilaterals for pressures.
% !             Mesh numbered in x- or y-direction. Freedoms numbered in the
% !             order u-p-v. Element by element solution using BiCGStab(l).
% !             with no preconditioning. No global matrix assembly,
% !--------------------------------------------------------------------------
%% !---------------------input and initialisation--------------------------
ndim = 2;
nip = 4;
nod = 8;
nodf = 4;                       % number of nodes per ﬂuid element
nodof = 3;
ntot = 20;                      % total number of degrees of freedom per element
one = 1;
penalty = 1e20;
pt5 = 0.5;
element = "quadrilateral";

nxe = 5;
nye = 5;
tol = 0.001;
limit = 30;
visc = 0.01;
rho = 1;
cg_tol = 1e-5;
cg_limit = 200;
x0 = 1;                         % initialisation value
ell = 4;                        % BiCGStab parameter; taken as 4 in this example
kappa = 0;                      % BiCGStab parameter; taken as zero in this example
[nels,nn] = mesh_size(element,nod,nxe,nye);

points = zeros(nip,ndim);
coord = zeros(nod,ndim);
derivf = zeros(ndim,nodf);      % ﬂuid shape function derivatives wrt global coordinates
uvel = zeros(nod,1);            % element nodal x-velocity
jac = zeros(ndim,ndim);
kay = zeros(ndim,ndim);
der = zeros(ndim,nod);
deriv = zeros(ndim,nod);
vvel = zeros(nod,1);            % element nodal y-velocity
derf = zeros(ndim,nodf);        % ﬂuid shape function derivatives wrt local coordinates
funf = zeros(nodf,1);           % ﬂuid shape functions
coordf = zeros(nodf,ndim);      % ﬂuid element nodal coordinates
g_g = zeros(ntot,nels);
c11 = zeros(nod,nod);
c12 = zeros(nod,nodf);
c21 = zeros(nodf,nod);
c23 = zeros(nodf,nod);
g = zeros(ntot,1);
c32 = zeros(nod,nodf);
ke = zeros(ntot,ntot);          % element ‘stiffness’ matrix
fun = zeros(nod,1);
x_coords = zeros(nxe+1,1);
y_coords = zeros(nye+1,1);
nf = zeros(nodof,nn);
g_coord = zeros(ndim,nn);
g_num = zeros(nod,nels);
num = zeros(nod,1);
weights = zeros(nip,1);
nd1 = zeros(nod,nod);           % product [fun] T [deriv(1,:)]
nd2 = zeros(nod,nod);           % product [fun] T [deriv(2,:)]
ndf1 = zeros(nod,nodf);         % product [fun] T [derivf(1,:)]
ndf2 = zeros(nod,nodf);         % product [fun] T [derivf(2,:)]
nfd1 = zeros(nodf,nod);         % product [funf] T [deriv(1,:)]
nfd2 = zeros(nodf,nod);         % product [funf] T [deriv(2,:)]
storke = zeros(ntot,ntot,nels);
s = zeros(ell+1,1);             % small local vector
gg = zeros(ell+1,ell+1);        % small local array
gamma = zeros(ell+1,1);         % small local array

x_coords(:,:) = 0:0.2:1;
y_coords(:,:) = 0:-0.2:-1;

kay(1,1) = visc/rho;
kay(2,2) = visc/rho;
nf(:,:) = 1;
nr = 81;
k = [1:19,21:2:27,28:36,38:2:44,45:55,57:2:61,62:70,72:2:78,79:96];
Temp = [1 1 0 0   2 1 0 0   3 1 1 0   4 1 0 0   5 1 1 0   6 1 0 0 ...
        7 1 1 0   8 1 0 0   9 1 1 0  10 1 0 0  11 1 1 0  12 0 0 0 ...
       13 1 0 1  14 1 0 1  15 1 0 1  16 1 0 1  17 0 0 0  18 0 1 0 ...
       19 1 0 1  21 1 0 1  23 1 0 1  25 1 0 1  27 1 0 1  28 0 1 0 ...
       29 0 0 0  30 1 0 1  31 1 0 1  32 1 0 1  33 1 0 1  34 0 0 0 ...
       35 0 1 0  36 1 0 1  38 1 0 1  40 1 0 1  42 1 0 1  44 1 0 1 ...
       45 0 1 0  46 0 0 0  47 1 0 1  48 1 0 1  49 1 0 1  50 1 0 1 ...
       51 0 0 0  52 0 1 0  53 1 0 1  54 1 1 1  55 1 0 1  57 1 0 1 ...
       59 1 0 1  61 1 0 1  62 0 1 0  63 0 0 0  64 1 0 1  65 1 0 1 ...
       66 1 0 1  67 1 0 1  68 0 0 0  69 0 1 0  70 1 0 1  72 1 0 1 ...        
       74 1 0 1  76 1 0 1  78 1 0 1  79 0 1 0  80 0 0 0  81 1 0 1 ...
       82 1 0 1  83 1 0 1  84 1 0 1  85 0 0 0  86 0 1 0  87 0 0 0 ...
       88 0 1 0  89 0 0 0  90 0 1 0  91 0 0 0  92 0 1 0  93 0 0 0 ...
       94 0 1 0  95 0 0 0  96 0 1 0];
for i = 1:nr
    ik = 4*i-3;
    nf(:,k(i)) = [Temp(ik+1);Temp(ik+2);Temp(ik+3)];
end
nf = formnf(nf);
neq = max(max(nf));
[points,weights] = sample(element,points,weights);
%% !---------------------loop the elements to find global arrays sizes-----
for iel = 1:nels
    [coord,num] = geom_rect(element,iel,x_coords,num,y_coords,"x");
    g(1:8) = nf(1,num(1:8));
    g(9:12) = nf(2,num(1:2:7));
    g(13:20) = nf(3,num(1:8));
    g_num(:,iel) = num;
    g_coord(:,num) = coord';
    g_g(:,iel) = g;
end
fprintf(" There are %d equations\n",neq)
loads = zeros(neq+1,1);
oldlds = zeros(neq+1,1);
rt = zeros(neq+1,1);            % initial residual vector
r = zeros(neq+1,ell+1);         % residual vector
u = zeros(neq+1,ell+1);         % gather/scatter arrays
b = zeros(neq+1,1);             % right-hand-side vector
diag = zeros(neq+1,1);          % diagonal of left-hand-side matrix
y = zeros(neq+1,1);             % gather/scatter arrays
y1 = zeros(neq+1,1);            % gather/scatter arrays


fixed_freedoms = 11;
node = zeros(fixed_freedoms,1);
sense = zeros(fixed_freedoms,1);
value = zeros(fixed_freedoms,1);
no = zeros(fixed_freedoms,1);
node(:,:) = [1,2,3,4,5,6,7,8,9,10,11];
value(:,:) = 1;
sense(:,:) = 1;
cg_tot = 0;
%% !---------------------iteration loop------------------------------------
for iters = 1:limit
    converged = false;
    ke(:,:) = 0;
    diag(:) = 0;
    b(:) = 0;
%% !---------------------element stiffness integration and storage --------
    for iel = 1:nels
        num(:,:) = g_num(:,iel);
        coord(:,:) = g_coord(:,num)';
        g = g_g(:,iel);
        g(g==0) = neq+1;
        coordf = coord(1:2:7,:);
        uvel(:,:) = (loads(g(1:nod)) + oldlds(g(1:nod)))*pt5;
        for i = nod+nodf+1:ntot
            vvel(i-nod-nodf) = (loads(g(i)) + oldlds(g(i)))*pt5;
        end
        c11(:,:) = 0;
        c12(:,:) = 0;
        c21(:,:) = 0;
        c23(:,:) = 0;
        c32(:,:) = 0;
%% !---------------------velocity contribution-----------------------------
        for i = 1:nip
            fun(:,:) = shape_fun(fun,points,i);
            der(:,:) = shape_der(der,points,i);
            jac = der*coord;
            dete = det(jac);
            deriv = jac\der;
            ubar = fun'*uvel;
            vbar = fun'*vvel;
            if iters == 1
                ubar = 1;
                vbar = 0;
            end
            nd1(:,:) = cross_product(fun,deriv(1,:)');
            nd2(:,:) = cross_product(fun,deriv(2,:)');
            c11 = c11 + dete*weights(i)*(deriv'*kay*deriv + nd1*ubar+nd2*vbar);
%% !---------------------pressure contribution-----------------------------
            funf(:,:) = shape_fun(funf,points,i);
            derf(:,:) = shape_der(derf,points,i);
            jac = derf*coordf;
            dete = det(jac);
            derivf = jac\derf;
            ndf1(:,:) = cross_product(fun,derivf(1,:)');
            ndf2(:,:) = cross_product(fun,derivf(2,:)');
            nfd1(:,:) = cross_product(funf,deriv(1,:)');
            nfd2(:,:) = cross_product(funf,deriv(2,:)');
            c12 = c12+ndf1*dete*weights(i)/rho;
            c32 = c32+ndf2*dete*weights(i)/rho;
            c21 = c21+nfd1*dete*weights(i);
            c23 = c23+nfd2*dete*weights(i);
        end
        ke = formupv(c11,c12,c21,c23,c32);
        storke(:,:,iel) = ke;
        for k = 1:ntot
            diag(g(k)) = diag(g(k)) + ke(k,k);
        end
    end
%% !---------------------specify pressure and velocity boundary values-----
    for i = 1:fixed_freedoms
        no(i) = nf(sense(i),node(i));
    end
    diag(no) = diag(no)+penalty;
    b(no) = diag(no).*value;
    store = diag(no);
%% !-solve the simultaneous equations element by element using BiCGStab(l)-
%% !---------------------initialisation phase------------------------------
    if iters == 1
        loads(1:end-1) = x0;
    end
    y = loads;
    y1(:) = 0;
    for iel = 1:nels
        g = g_g(:,iel);
        g(g==0) = neq+1;
        ke = storke(:,:,iel);
        y1(g) = y1(g) + ke*y(g);
    end
    y1(end) = 0;
    y1(no) = y(no).*store;
    y = y1;
    rt = b-y;
    r(:,:) = 0;
    r(:,1) = rt;
    u(:) = 0;
    gama = 1;
    omega = 1;
    k = 0;
    norm_r = norm(rt);
    r0_norm = norm_r;
    error = 1;
%% !---------------------BiCGStab(l) iterations----------------------------
    for cg_iters = 1:cg_limit
        cg_converged = (error < cg_tol);
        if cg_converged 
            break;
        end
        gama = -omega*gama;
        y = r(:,1);
        for j = 1:ell
            rho1 = rt'*y;
            beta = rho1/gama;
            u(:,1:j) = r(:,1:j) - beta*u(:,1:j);
            y = u(:,j);
            y1(:) = 0;
            for iel = 1:nels
                g = g_g(:,iel);
                g(g==0) = neq+1;
                ke = storke(:,:,iel);
                y1(g) = y1(g) + ke*y(g);
            end
            y1(end) = 0;
            y1(no) = y(no).*store;
            y = y1; 
            u(:,j+1) = y;
            gama = rt'*y;
            alpha = rho1/gama;
            loads = loads + alpha*u(:,1);
            r(:,1:j) = r(:,1:j)-alpha*u(:,2:j+1);
            y = r(:,j);
            y1(:) = 0;
            for iel = 1:nels
                g = g_g(:,iel);
                g(g==0) = neq+1;
                ke = storke(:,:,iel);
                y1(g) = y1(g) + ke*y(g);
            end
            y1(end) = 0;
            y1(no) = y(no).*store;
            y = y1; 
            r(:,j+1) = y;
        end
        gg = r'*r;
        [omega,gamma,s] = form_s(gg,ell,kappa,omega,gamma,s);
        loads = loads - r*s;
        r(:,1) = r*gamma;
        u(:,1) = u*gamma;
        norm_r = norm(r(:,1));
        error = norm_r/r0_norm;
        k = k + 1;
    end
    cg_tot = cg_tot + cg_iters;
%% !----------------------end of BiCGStab(l) process-----------------------
    [converged,oldlds] = checon(loads,oldlds,tol);
    if converged
        break
    end 
end
u = zeros(nn,1);
v = zeros(nn,1);
fprintf("  Node     u-velocity  pressure    v-velocity\n")
for k = 1:nn
    nf(nf==0) = neq+1;
    u(k) = loads(nf(1,k));
    v(k) = loads(nf(3,k));
    fprintf("%2d %13.4e %13.4e %13.4e\n",k,u(k),loads(nf(2,k)),v(k))
end
fprintf(" Converged in %d  iterations\n with an average of %2d BiCGStab(l) iterations\n",...
            iters,floor(cg_tot/iters))
quiver(g_coord(1,12:end),g_coord(2,12:end),u(12:end)',v(12:end)')




