% PROGRAM p83
% !-------------------------------------------------------------------------
% ! Program 8.3 One dimensional consolidation analysis using 2-node "rod"
% !             elements. Explicit time integration. Element-by-element.
% !             Lumped mass.
% !-------------------------------------------------------------------------
%% !---------------------input and initialisation--------------------------
nod = 2;
nprops = 2;
penalty = 1e20;
pt5 = 0.5;

nels = 10;
np_types = 1;
neq = nels+1;
num = zeros(nod,1);
etypes = zeros(nels,1);
kc = zeros(nod,nod);
mm = zeros(nod,nod);
press = zeros(neq);
prop = zeros(nprops,np_types);
ell = zeros(nels,1);
loads = zeros(neq,1);
newlo = zeros(neq,1);           % new excess pore pressure values
mass = zeros(nod,1);            % element lumped mass vector
globma = zeros(neq,1);          % global lumped mass matrix (stored as a vector)
store_mm = zeros(nod,nod,nels); % stores lhs element matrices

prop(:,:) = [1,10];
etypes(:,:) = 1;
ell(:,:) = 0.1*ones(1,nels);
dtim = 0.001;                   % calculation time step
nstep = 2000;                   % number of time steps required
npri = 100;                     % output printed every npri time steps
nres = 11;                      % node number at which time history is to be printed
ntime = 1000;                   % time step number at which spatial distribution is to be printed or contoured
fprintf(" There are %d equations\n", neq);
%% !---------------------global conductivity and "mass" matrix assembly----
for iel = 1:nels
    num(:,1) = [iel,iel+1];
    kc = rod_km(kc,prop(1,etypes(iel)),ell(iel));
    mm(:,:) = 0;
    for i = 1:nod
        mm(i,i) = ell(iel)/2;
        mass(i) = ell(iel)/2;
    end
    store_mm(:,:,iel) = mm-kc*dtim;
    globma(num) = globma(num)+mass;
end
%% !---------------------specify initial and boundary values---------------
loads(:,1) = 100*ones(neq,1);
% loads(:,1) = 10*[10:-1:0];
fixed_freedoms = 1;
globma = 1./globma;
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    node(:,1) = 1;
    value = zeros(fixed_freedoms,1);
    value(:,1) = 0;
end
%% !---------------------time stepping loop--------------------------------
fprintf("    Time           Uav          Pressure (node %2d)\n",nres)
fprintf("%13.4e %13.4e %13.4e\n",0,0,loads(nres))
a0 = 0;                         % holds area beneath isochrone by trapezoid rule at time t =0
for iel = 1:nels
    a0 = a0+pt5*ell(iel)*(loads(iel)+loads(iel+1));
end
for j = 1:nstep
    time = j*dtim;
    newlo(:) = 0;
    for iel = 1:nels
        num(:,1) = [iel,iel+1];
        mm = store_mm(:,:,iel);
        newlo(num) = newlo(num) + mm*loads(num);
    end
    loads = newlo.*globma;
    if fixed_freedoms ~= 0
        loads(node) = value;
    end
    
    at = 0;
    for iel = 1:nels
        at = at+pt5*ell(iel)*(loads(iel)+loads(iel+1));
    end
    if j == ntime
        press(:,1) = loads(:,1);
    end
%     plot([0:10],loads,'LineWidth',2)
%         axis([0 10 0 100])
%         drawnow
    if mod(j,npri) == 0
        fprintf("%13.4e %13.4e %13.4e\n",time,(a0-at)/a0,loads(nres))
    end
end
fprintf("    Depth       Pressure (time= %.4e)\n",ntime*dtim)
fprintf("%13.4e %13.4e\n",0,press(1))
for i = 1:nels
    fprintf("%13.4e %13.4e\n",sum(ell(1:i)),press(i+1))
end








