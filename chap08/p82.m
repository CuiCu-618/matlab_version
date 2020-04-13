% PROGRAM p82
% !-------------------------------------------------------------------------
% ! Program 8.2 One dimensional consolidation analysis
% ! (settlement and excess pore pressure) using 2-node "line" elements. 
% ! Implicit time integration using the "theta" method.
% !-------------------------------------------------------------------------
%% !---------------------input and initialisation--------------------------
nod = 2;
nprops = 2;
penalty = 1e20;
pt5 = 0.5;

nels = 40;
np_types = 4;

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
newlo = zeros(neq,1);          

prop(:,:) = [2.78e-11 6.41e-5;8.25e-11 4.08e-5;1.17e-11 2.04e-5;2.94e-11 4.08e-5]';
gamw = 9.81;
etypes(:,:) = 1;
etypes(:,1) = [ones(1,10),2*ones(1,10),3*ones(1,10),4*ones(1,10)];

ell(:,1) = [0.305*ones(1,10),0.61*ones(1,10),0.914*ones(1,10),0.61*ones(1,10)];
dtim = 86400;
nstep = 7200;
theta = 0.5;
npri = 400;
nres = 21;
ntime = 740;
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
    kc = rod_km(kc,prop(1,etypes(iel))/gamw,ell(iel));
    mm = rod_mm(mm,ell(iel)*prop(2,etypes(iel)));
    kv = fsparv(kv,kc,num,kdiag);
    bp = fsparv(bp,mm,num,kdiag);
end
kv = kv*theta*dtim;
bp = bp+kv;
kv = bp-kv/theta;
%% !---------------------specify initial and boundary values---------------
loads(:,1) = 100*ones(neq,1);
a0 = 0;                         
sc = 0;                         % ultimate consolidation settlement
for iel = 1:nels
    a0 = a0+pt5*ell(iel)*(loads(iel)+loads(iel+1));
    sc = sc+pt5*ell(iel)*prop(2,etypes(iel))*(loads(iel)+loads(iel+1));
end
fixed_freedoms = 2;
if fixed_freedoms ~= 0
    node = zeros(fixed_freedoms,1);
    node(:,1) = [1,41];
    value = zeros(fixed_freedoms,1);
    value(:,1) = [0,0];
    storbp  = zeros(fixed_freedoms,1); 
    bp(kdiag(node)) = bp(kdiag(node)) + penalty;
    storbp(:,1) = bp(kdiag(node));
end
%% !---------------------factorise equations-------------------------------
bp = sparin(bp,kdiag);
%% !---------------------time stepping loop--------------------------------
fprintf("    Time           Uav          Uavs        Settlement     Pressure (node %2d)\n",nres)
fprintf("%13.4e %13.4e %13.4e %13.4e %13.4e\n",0,0,0,0,loads(nres))

for j = 1:nstep
    time = j*dtim;
    newlo(:,1) = linmul_sky(kv,loads,kdiag);
    if fixed_freedoms ~= 0
        newlo(node) = storbp.*value;
    end
    newlo = spabac(bp,newlo,kdiag);
    loads = newlo;
    
    at = 0;
    sl = 0;
    for iel = 1:nels
        at = at+pt5*ell(iel)*(loads(iel)+loads(iel+1));
        sl = sl+pt5*ell(iel)*prop(2,etypes(iel))*(loads(iel)+loads(iel+1));
    end
    uav = (a0-at)/a0;
    uavs = (sc-sl)/sc;
    if j == ntime
        press(:,1) = loads(:,1);
    end
    
    if mod(j,npri) == 0
        plot(loads,'LineWidth',2)
        axis([1 41 0 100])
        drawnow
        fprintf("%13.4e %13.4e %13.4e %13.4e %13.4e\n",time,uav,uavs,sc-sl,loads(nres))
    end
end
fprintf("    Depth       Pressure (time= %.4e)\n",ntime*dtim)
fprintf("%13.4e %13.4e\n",0,press(1))
for i = 1:nels
    fprintf("%13.4e %13.4e\n",sum(ell(1:i)),press(i+1))
end





