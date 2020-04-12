% PROGRAM p75
% !-------------------------------------------------------------------------
% ! Program 7.5 General two- (plane) or three-dimensional analysis of steady
% !     seepage. No global conductivity matrix assembly.
% !     Diagonally preconditioned conjugate gradient solver.
% !-------------------------------------------------------------------------
%% !---------------------input and initialisation--------------------------
one = 1;
penalty = 1e20;

element = "hexahedron";
nod = 8;
nels = 64;
nn = 125;
nip = 8;
ndim = 3;
cg_tol = 1e-7;
cg_limit = 200;
np_types = 2;
neq = nn;
fprintf(" There are %d equations\n",neq);

points = zeros(nip,ndim);
g_coord = zeros(ndim,nn);
coord = zeros(nod,ndim);
etype = zeros(nels,1);
jac = zeros(ndim,ndim);
weights = zeros(nip,1);
num = zeros(nod,1);
g_num = zeros(nod,nels);
der = zeros(ndim,nod);
deriv = zeros(ndim,nod);
kp = zeros(nod,nod);
kay = zeros(ndim,ndim);
prop = zeros(ndim,np_types);
p = zeros(neq,1);
loads = zeros(neq,1);
x = zeros(neq,1);
xnew = zeros(neq,1);
u = zeros(neq,1);
diag_precon = zeros(neq,1);
d = zeros(neq,1);
disps = zeros(neq,1);
storkp = zeros(nod,nod,nels);

prop(:,:) = [2,2,2;1,1,5]';
etype(:,:) = [ones(1,32),2*ones(1,32)];
g_coord(:,:) = [0.0  0.0  0.0;1.0  0.0  0.0;2.0  0.0  0.0;3.0  0.0  0.0
                4.0  0.0  0.0;0.0  0.0 -1.0;1.0  0.0 -1.0;2.0  0.0 -1.0
                3.0  0.0 -1.0;4.0  0.0 -1.0;0.0  0.0 -2.0;1.0  0.0 -2.0
                2.0  0.0 -2.0;3.0  0.0 -2.0;4.0  0.0 -2.0;0.0  0.0 -3.0
                1.0  0.0 -3.0;2.0  0.0 -3.0;3.0  0.0 -3.0;4.0  0.0 -3.0
                0.0  0.0 -4.0;1.0  0.0 -4.0;2.0  0.0 -4.0;3.0  0.0 -4.0
                4.0  0.0 -4.0;0.0  1.0  0.0;1.0  1.0  0.0;2.0  1.0  0.0
                3.0  1.0  0.0;4.0  1.0  0.0;0.0  1.0 -1.0;1.0  1.0 -1.0
                2.0  1.0 -1.0;3.0  1.0 -1.0;4.0  1.0 -1.0;0.0  1.0 -2.0
                1.0  1.0 -2.0;2.0  1.0 -2.0;3.0  1.0 -2.0;4.0  1.0 -2.0
                0.0  1.0 -3.0;1.0  1.0 -3.0;2.0  1.0 -3.0;3.0  1.0 -3.0
                4.0  1.0 -3.0;0.0  1.0 -4.0;1.0  1.0 -4.0;2.0  1.0 -4.0
                3.0  1.0 -4.0;4.0  1.0 -4.0;0.0  2.0  0.0;1.0  2.0  0.0
                2.0  2.0  0.0;3.0  2.0  0.0;4.0  2.0  0.0;0.0  2.0 -1.0
                1.0  2.0 -1.0;2.0  2.0 -1.0;3.0  2.0 -1.0;4.0  2.0 -1.0
                0.0  2.0 -2.0;1.0  2.0 -2.0;2.0  2.0 -2.0;3.0  2.0 -2.0
                4.0  2.0 -2.0;0.0  2.0 -3.0;1.0  2.0 -3.0;2.0  2.0 -3.0
                3.0  2.0 -3.0;4.0  2.0 -3.0;0.0  2.0 -4.0;1.0  2.0 -4.0
                2.0  2.0 -4.0;3.0  2.0 -4.0;4.0  2.0 -4.0;0.0  3.0  0.0
                1.0  3.0  0.0;2.0  3.0  0.0;3.0  3.0  0.0;4.0  3.0  0.0
                0.0  3.0 -1.0;1.0  3.0 -1.0;2.0  3.0 -1.0;3.0  3.0 -1.0
                4.0  3.0 -1.0;0.0  3.0 -2.0;1.0  3.0 -2.0;2.0  3.0 -2.0
                3.0  3.0 -2.0;4.0  3.0 -2.0;0.0  3.0 -3.0;1.0  3.0 -3.0
                2.0  3.0 -3.0;3.0  3.0 -3.0;4.0  3.0 -3.0;0.0  3.0 -4.0
                1.0  3.0 -4.0;2.0  3.0 -4.0;3.0  3.0 -4.0;4.0  3.0 -4.0
                0.0  4.0  0.0;1.0  4.0  0.0;2.0  4.0  0.0;3.0  4.0  0.0
                4.0  4.0  0.0;0.0  4.0 -1.0;1.0  4.0 -1.0;2.0  4.0 -1.0
                3.0  4.0 -1.0;4.0  4.0 -1.0;0.0  4.0 -2.0;1.0  4.0 -2.0
                2.0  4.0 -2.0;3.0  4.0 -2.0;4.0  4.0 -2.0;0.0  4.0 -3.0
                1.0  4.0 -3.0;2.0  4.0 -3.0;3.0  4.0 -3.0;4.0  4.0 -3.0
                0.0  4.0 -4.0;1.0  4.0 -4.0;2.0  4.0 -4.0;3.0  4.0 -4.0
                4.0  4.0 -4.0]';
g_num(:,:) = [6   1   2   7  31  26  27  32;7   2   3   8  32  27  28  33
              8   3   4   9  33  28  29  34;9   4   5  10  34  29  30  35
             11   6   7  12  36  31  32  37;12   7   8  13  37  32  33  38
             13   8   9  14  38  33  34  39;14   9  10  15  39  34  35  40
             16  11  12  17  41  36  37  42;17  12  13  18  42  37  38  43
             18  13  14  19  43  38  39  44;19  14  15  20  44  39  40  45
             21  16  17  22  46  41  42  47;22  17  18  23  47  42  43  48
             23  18  19  24  48  43  44  49;24  19  20  25  49  44  45  50
             31  26  27  32  56  51  52  57;32  27  28  33  57  52  53  58
             33  28  29  34  58  53  54  59;34  29  30  35  59  54  55  60
             36  31  32  37  61  56  57  62;37  32  33  38  62  57  58  63
             38  33  34  39  63  58  59  64;39  34  35  40  64  59  60  65
             41  36  37  42  66  61  62  67;42  37  38  43  67  62  63  68
             43  38  39  44  68  63  64  69;44  39  40  45  69  64  65  70
             46  41  42  47  71  66  67  72;47  42  43  48  72  67  68  73
             48  43  44  49  73  68  69  74;49  44  45  50  74  69  70  75
             56  51  52  57  81  76  77  82;57  52  53  58  82  77  78  83
             58  53  54  59  83  78  79  84;59  54  55  60  84  79  80  85
             61  56  57  62  86  81  82  87;62  57  58  63  87  82  83  88
             63  58  59  64  88  83  84  89;64  59  60  65  89  84  85  90
             66  61  62  67  91  86  87  92;67  62  63  68  92  87  88  93
             68  63  64  69  93  88  89  94;69  64  65  70  94  89  90  95
             71  66  67  72  96  91  92  97;72  67  68  73  97  92  93  98
             73  68  69  74  98  93  94  99;74  69  70  75  99  94  95 100
             81  76  77  82 106 101 102 107;82  77  78  83 107 102 103 108
             83  78  79  84 108 103 104 109;84  79  80  85 109 104 105 110
             86  81  82  87 111 106 107 112;87  82  83  88 112 107 108 113
             88  83  84  89 113 108 109 114;89  84  85  90 114 109 110 115
             91  86  87  92 116 111 112 117;92  87  88  93 117 112 113 118
             93  88  89  94 118 113 114 119;94  89  90  95 119 114 115 120
             96  91  92  97 121 116 117 122;97  92  93  98 122 117 118 123
             98  93  94  99 123 118 119 124;99  94  95 100 124 119 120 125]';
[points,weights] = sample(element,points,weights);
%% !--------element conductivity integration, storage and preconditioner--- 
for iel = 1:nels
    kay(:,:) = 0;
    for i = 1:ndim
        kay(i,i) = prop(i,etype(iel));
    end
    num = g_num(:,iel);
    coord = g_coord(:,num)';
    kp(:,:) = 0;
    for i = 1:nip
        der = shape_der(der,points,i);
        jac = der*coord;
        dete = det(jac);
        deriv = jac\der;
        kp = kp + deriv'*kay*deriv*dete*weights(i);
    end
    storkp(:,:,iel) = kp;
    for k = 1:nod
        diag_precon(num(k)) = diag_precon(num(k)) + kp(k,k);
    end
end
%% !---------------------invert the preconditioner and get starting loads--
loaded_nodes = 1;
k = 21;
loads(k) = 100;
fixed_freedoms = 61;
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    node(:,:) = [  1 ,2 ,3 ,4 ,5 ,10 ,15 ,20 ,25 ,26 ,27 ,28 ,29 ,30 ,35 ,40 ,...
                  45 ,50 ,51 ,52 ,53 ,54 ,55 ,60 ,65 ,70 ,75 ,76 ,77 ,78 ,79 ,...
                  80 ,85 ,90 ,95 ,100 ,101 ,102 ,103 ,104 ,105 ,106 ,107 ,108 ,...
                  109 ,110 ,111 ,112 ,113 ,114 ,115 ,116 ,117 ,118 ,119 ,120 ,...
                  121 ,122 ,123 ,124 ,125];
    value = zeros(fixed_freedoms,1);
    value(:,:) = 0;
    store = zeros(fixed_freedoms,1);
    diag_precon(node)=diag_precon(node)+penalty;
    loads(node)=diag_precon(node).*value; 
    store(:,1)=diag_precon(node);
end
diag_precon = one./diag_precon;
d = diag_precon.*loads;
p = d;
x(:,:) = 0;
%%  !---------------------pcg equation solution-----------------------------
for cg_iters = 1:cg_limit
    u(:,:) = 0;                 % Q
    for iel = 1:nels
        num = g_num(:,iel);
        kp = storkp(:,:,iel);
        u(num) = u(num) + kp*p(num);
    end
    if fixed_freedoms ~= 0
        u(node) = p(node).*store;
    end
    up = loads'*d;
    alpha = up/(p'*u);
    xnew(:,:) = x(:,:) + alpha*p;
    loads = loads - u*alpha;    % R
    d = diag_precon.*loads;      % R
    beta = loads'*d/up; 
    p = d+p*beta;
    [cg_converged,x] = checon(xnew,x,cg_tol);
    if cg_converged
        break;
    end
end
fprintf(" Number of cg iterations to convergence was %2d\n",cg_iters)
%% !---------------------retrieve nodal net flow rates---------------------
loads(:,:) = xnew;
for iel = 1:nels
    num = g_num(:,iel);
    kp = storkp(:,:,iel);
    disps(num) = disps(num) + kp*loads(num);
end
fprintf("  Node Total Head  Flow rate \n")
for k = 1:nn
    fprintf("%3d %13.4e %13.4e\n",k,loads(k),disps(k))
end
fprintf("     Inflow      Outflow\n")
fprintf("%13.4e %13.4e\n",sum(disps(disps>0)),sum(disps(disps<0)))















