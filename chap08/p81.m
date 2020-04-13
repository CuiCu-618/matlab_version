% PROGRAM p81
% !-------------------------------------------------------------------------
% ! Program 8.1 One dimensional consolidation analysis using 2-node "rod"
% !             elements. Implicit time integration using the "theta" method.
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
kdiag = zeros(neq,1);
loads = zeros(neq,1);
newlo = zeros(neq,1);           % new excess pore pressure values

prop(:,:) = [1,10];
etypes(:,:) = 1;
etypes([4,5,6],1) = 2;
ell(:,:) = 0.1*ones(1,nels);
dtim = 0.001;                   % calculation time step
nstep = 2000;                   % number of time steps required
theta = 0.5;                    % time-integration weighting parameter
npri = 100;                     % output printed every npri time steps
nres = 11;                      % node number at which time history is to be printed
ntime = 1000;                   % time step number at which spatial distribution is to be printed or contoured
%% !---------------------loop the elements to find global arrays sizes-----
for iel = 1:nels
    num(:,1) = [iel,iel+1];
    kdiag = fkdiag(kdiag,num);
end
for i = 2:neq
    kdiag(i) = kdiag(i) + kdiag(i-1);
end
kv = zeros(kdiag(neq),1);       % global conductivity matrix
bp = zeros(kdiag(neq),1);       % global ‘mass’ matrix
fprintf(" There are %d equations and the skyline storage is %d \n",...
                neq,kdiag(neq));
%% !---------------------global conductivity and "mass" matrix assembly----
for iel = 1:nels
    num(:,1) = [iel,iel+1];
    kc = rod_km(kc,prop(etypes(iel),1),ell(iel));
    mm = rod_mm(mm,ell(iel));
    kv = fsparv(kv,kc,num,kdiag);
    bp = fsparv(bp,mm,num,kdiag);
end
kv = kv*theta*dtim;
bp = bp+kv;
kv = bp-kv/theta;
%% !---------------------specify initial and boundary values---------------
loads(:,1) = 100*ones(neq,1);
% loads(:,1) = 10*[10:-1:0];
a0 = 0;                         % holds area beneath isochrone by trapezoid rule at time t =0
for iel = 1:nels
    a0 = a0+pt5*ell(iel)*(loads(iel)+loads(iel+1));
end
fixed_freedoms = 1;
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    node(:,1) = 1;
    value = zeros(fixed_freedoms,1);
    value(:,1) = 0;
    storbp  = zeros(fixed_freedoms,1); % stores augmented diagonal terms
    bp(kdiag(node)) = bp(kdiag(node)) + penalty;
    storbp(:,1) = bp(kdiag(node));
end
%% !---------------------factorise equations-------------------------------
bp = sparin(bp,kdiag);
%% !---------------------time stepping loop--------------------------------
fprintf("    Time           Uav          Pressure (node %2d)\n",nres)
fprintf("%13.4e %13.4e %13.4e\n",0,0,loads(nres))

for j = 1:nstep
    time = j*dtim;
    newlo(:,1) = linmul_sky(kv,loads,kdiag);
    if fixed_freedoms ~= 0
        newlo(node) = storbp.*value;
    end
    newlo = spabac(bp,newlo,kdiag);
    loads = newlo;
    
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

