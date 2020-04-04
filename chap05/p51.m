% PROGRAM p51
% !------------------------------------------------------------------------
% ! Program 5.1 Plane or axisymmetric strain analysis of an elastic solid
% !             using 3-, 6-, 10- or 15-node right-angled triangles or
% !             4-, 8- or 9-node rectangular quadrilaterals. Mesh numbered
% !             in x(r)- or y(z)- direction.
% !------------------------------------------------------------------------
%% ---------------------------initialisation-------------------------------
fixed_freedoms = 0;
loaded_nodes = 3;
ndim = 2;
nip = 1;                        % number of integration points per element
nod = 3;
nodof = 2;
nprops = 2;
np_types = 1;
nr = 5;
nst = 3;                        % number of stress (strain) terms
nxe = 2;
nye = 2;

type_2d = "plane";
element = "triangle";
dir = "x";
[nels,nn] = mesh_size(element,nod,nxe,nye);
ndof = nod*nodof;
if type_2d == "'axisymmetric"
    nst = 4;
end

nf = zeros(nodof,nn);
points = zeros(nip,ndim);       % integrating point local coordinates
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
y_coords = zeros(nxe+1,1);      % y-coordinates of mesh layout
etype = zeros(nels,1);
gc = zeros(ndim,1);             % integrating point coordinates
dee = zeros(nst,nst);           % stress strain matrix
sigma = zeros(nst,1);           % stress terms

prop(:,:) = [1e6;0.3];
etype(:) = 1;
x_coords(:,1) = [0,0.5,1]';
y_coords(:,1) = [0,-0.5,-1]';
nf(:,:) = 1;
k = [1,4,7,8,9];
nf(:,k) = [0,0,0,1,1;1,1,0,0,0];
nf = formnf(nf);









