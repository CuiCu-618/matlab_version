% PROGRAM p56         
% !-------------------------------------------------------------------------
% ! Program 5.6 Three-dimensional strain of an elastic solid using
% !             8-, 14- or 20-node brick hexahedra. Mesh numbered in x-z
% !             planes then in the y-direction. No global stiffness matrix
% !             assembly. Diagonally preconditioned conjugate gradient solver.
% !-------------------------------------------------------------------------
%% ---------------------------initialisation-------------------------------
cg_limit = 200;                 % pcg iteration ceiling
fixed_freedoms = 0;
loaded_nodes = 8;
ndim = 3;
nip = 8;
nprops = 2;
np_types = 2;
nod = 20;
nodof = 3;
nr = 46;
nst = 6;
nxe = 1;
nye = 3;
nze = 2;
penalty=1e20;
element = "hexahedron";
cg_converged = false;

[nels,nn] = mesh_size(element,nod,nxe,nye,nze);
ndof = nod*nodof;

nf = zeros(nodof,nn);
points = zeros(nip,ndim);
g = zeros(ndof,1);
g_coord = zeros(ndim,nn);
fun = zeros(nod,1);             % shape functions
coord = zeros(nod,ndim);
jac = zeros(ndim,ndim);         % Jacobian matrix
g_num = zeros(nod,nels);
der = zeros(ndim,nod);          % shape function derivatives with respect to local coordinates
deriv = zeros(ndim,nod);        % shape function derivatives with respect to global coordinates
bee = zeros(nst,ndof);          % strain–displacement matrix
km = zeros(ndof,ndof);
eld = zeros(ndof,1);
weights = zeros(nip,1);         % weighting coefﬁcients
g_g = zeros(ndof,nels);
prop = zeros(nprops,np_types);
num = zeros(nod,1);
x_coords = zeros(nxe+1,1);      % x-coordinates of mesh layout
y_coords = zeros(nye+1,1);      % y-coordinates of mesh layout
z_coords = zeros(nze+1,1);      % z-coordinates of mesh layout
etype = zeros(nels,1);
gc = zeros(ndim,1);             % integrating point coordinates
dee = zeros(nst,nst);           % stress strain matrix
sigma = zeros(nst,1);           % stress terms
storkm = zeros(ndof,ndof,nels); % holds element stiffness matrices

prop(:,:) = [100,50;0.3,0.3];
etype(:,:) = [1,2,1,2,1,2];
x_coords(:,:) = [0,0.5];
y_coords(:,:) = [0,1,2,3];
z_coords(:,:) = [0,-1,-2];
nf(:,:) = 1;
k = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,19,20,23,25,28,30,31,32,...
     33,35,37,38,39,42,44,47,49,50,51,52,54,56,57,58,61,63,66,68,69,70];
nf(:,k) = [0,1,1,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
           0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,0,0;...
           1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,0,0];
nf = formnf(nf);
neq = max(max(nf));
loads = zeros(neq+1,1);
kdiag = zeros(neq,1);
p = zeros(neq+1,1);             % ‘descent’ vector used in equations (3.22)
x = zeros(neq+1,1);             % ‘old’ solution vectors
xnew = zeros(neq+1,1);          % ‘new’ solution vectors
u = zeros(neq+1,1);             % vector used in equation (3.22)
diag_precon = zeros(neq+1,1);   % diagonal preconditioner vector
d = zeros(neq+1,1);             % preconditioned rhs vector
fprintf(" There are %d equations\n",neq);
[points,weights] = sample(element,points,weights);
%% !--------element stiffness integration, storage and preconditioner------ 
for iel = 1:nels
    [coord,num] = hexahedron_xz(iel,x_coords,y_coords,z_coords,num);
    g(:,1) = num_to_g(num,nf);
    g_num(:,iel) = num;
    g_coord(:,num) = coord';
    g_g(:,iel) = g;
    kdiag = fkdiag(kdiag,g);
    dee = deemat(dee,prop(1,etype(iel)),prop(2,etype(iel)));
    num(:,1) = g_num(:,iel);
    coord(:,:) = g_coord(:,num)';
    g(:,1) = g_g(:,iel);
    km(:,:) = 0;
    for i = 1:nip
        der(:) = shape_der(der,points,i);
        jac = der*coord;
        dete = det(jac);
        deriv = jac\der;
        bee = beemat(bee,deriv);
        km = km + bee'*dee*bee*dete*weights(i);
    end
    storkm(:,:,iel) = km;
    for k = 1:ndof
        diag_precon(g(k)) = diag_precon(g(k)) + km(k,k);
    end
end
%% !---------------------invert the preconditioner and get starting loads--
k = [1,2,3,14,15,20,21,22];
nf(nf == 0) = neq+1;
loads(nf(:,k)) = [zeros(2,8);0.0417,-0.1667,0.0417,-0.1667,-0.1667,0.0417,-0.1667,0.0417];
loads(neq+1) = 0;
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    no = zeros(fixed_freedoms,1);
    sense = zeros(fixed_freedoms,1);
    value = zeros(fixed_freedoms,1);
    store = zeros(fixed_freedoms,1);    % stores augmented diagonal terms
    for i = 1:fixed_freedoms
        no(i) = nf(sense(i),node(i));
    end
    diag_precon(no) = diag_precon(no) + penalty;
    loads(no) = diag_precon(no) .* value;
    store(no) = diag_precon(no);        % ?
end
diag_precon(1:end-1) = 1./diag_precon(1:end-1);
diag_precon(end) = 0;
d = diag_precon*loads;
p = d;
cg_iters = 0;                   % pcg iteration counter
%% !---------------------pcg equation solution-----------------------------
for cg_iters = 1:cg_limit
    u(:,:) = 0;
    for iel = 1:nels
        g = g_g(:,iel);
        km = storkm(:,:,iel);
        u(g) = u(g) + km*p(g);
    end
end




