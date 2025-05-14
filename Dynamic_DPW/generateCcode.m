function [] = generateCcode()

%================Geometric Data============================================
% FINGER PARAMETERS
rb = 5;

Hv = [4.0];
Ch_L = 6.5*ones(size(Hv));
Unit_sep = 1.0*ones(size(Hv));

t = 0.8;
t_lim = 1.25;

% PROBLEM PARAMETER
t_ch = 1.0;
number_of_segments = length(Hv);

E =  26.0;
nu = 0.35;

base = 2*rb;
UC = 1.5*base;
R = (Hv.^2 + rb^2)./(2*Hv);

start_angle = 0;
BC_nodes = [1 2 3];

segment_length = Ch_L + t_ch + Unit_sep;
segment_height = UC/2 + t_lim + t;

geo_matrix =  DPG_Geometry(segment_height,segment_length,number_of_segments);
coords0 = geo_matrix.coords0;
connecRigid = geo_matrix.connecRigid;
connecBistable = geo_matrix.connecBistable;
connecTorsion = geo_matrix.connecTorsion;
Nnodes = geo_matrix.Nnodes;
connecNodes = geo_matrix.connecNodes;

% ================Spring Data==============================================
k_lim = E*(t_lim*UC)./(Ch_L+t+t_ch); % Linear Spring
k_s = (E*(t*UC)./segment_height)*ones(1,number_of_segments); % Rigid Links

D = (E*(t_lim)^3)/(12*(1-nu^2));
kt = D*(UC./(Unit_sep+t+t_ch));                           % Torsional Spring

Ks = [k_s',k_lim',kt'];

% Bistable Spring
kBS = (E ./ R.^2) .* ( ...
    - (1.97 .* R .* t.^5) ./ Hv.^3 + ...
    (7.4 .* Hv.^2 .* t.^3) ./ R.^2 + ...
    (3.5 .* R .* t.^4) ./ Hv.^2 + ...
    (0.37 .* Hv.^2 .* t.^2) ./ R + ...
    (42.2 .* t.^5) ./ Hv.^2 - ...
    (35.8 .* Hv .* t.^4) ./ R.^2 + ...
    (71.8 .* t.^5) ./ (Hv .* R) - ...
    (3.4 .* R .* t.^3) ./ Hv - ...
    (67.7 .* t.^4) ./ Hv + ...
    4.2 .* R .* t.^2 + ...
    11.1 .* t.^3);

al = (((0.2941 .* Hv.^3) ./ R.^3 - ...
    (0.4412 .* t.^3) ./ Hv.^3 - ...
    (2.2059 .* Hv.^2 .* t) ./ R.^3 - ...
    (0.7353 .* Hv.^2) ./ R.^2 - ...
    (6.9118 .* t.^3) ./ (Hv.^2 .* R) + ...
    (1.1765 .* t.^2) ./ Hv.^2 - ...
    (3.0882 .* t.^3) ./ (Hv .* R.^2) + ...
    (7.0588 .* Hv .* t) ./ R.^2 + ...
    (13.3824 .* t.^2) ./ (Hv .* R) + ...
    (0.4412 .* Hv) ./ R + ...
    (13.6765 .* t.^3) ./ R.^3 - ...
    (15.0 .* t.^2) ./ R.^2 - ...
    (4.7059 .* t) ./ R));

tip_dis =  0.2565 * t + 2.1425 * Hv;

dd = tip_dis - Hv - Unit_sep - t_ch/2 - t/2;

for i = 1:length(dd)
    if dd(i)<0
        dd(i) = 0.01;
    end
end


Bis_S = [kBS',al',dd'];

% ================Dynamic Data ============================================
eta_iso = 0.01;
eta = 0.05;

% =========================== Spatial Symbolic Setup=======================
NoBcNodes = Nnodes-length(BC_nodes);
numVar = NoBcNodes*2;

%===========================Dynamic Setup=========================
% PRESSURE PARAMETERS
P0 = 0.325; % Presure [MPa]
k_dome = 1.0;

tau_1 = 1.0;
tau_2 = 4.0;
tau_3 = 5.0;

% MASS PARAMETERS
DPG_mass = (7*10^-6)*(number_of_segments/5); %[Ton]
mass_V = (DPG_mass/Nnodes)*ones(numVar,1);

% =========================================================================
% Integration scheme for time (RK45)
% =========================================================================
x0s = zeros(numVar,1);
v0 = zeros(numVar,1);
x0 = [x0s;v0];

TT = 6.0;
tspan = [0, TT];


% C++ for Solve system
x0sType = coder.typeof(x0s, [Inf,1], [1,0]);
v0Type = coder.typeof(v0, [Inf,1], [1,0]);
x0Type = coder.typeof(x0, [Inf,1], [1,0]);
mass_VType = coder.typeof(mass_V, [Inf,1], [1,0]);
KsType = coder.typeof(Ks, [Inf,3], [1,0]);
Bis_SType = coder.typeof(Bis_S, [Inf,3], [1,0]);
connecNodesType = coder.typeof(connecNodes, [Inf,2], [1,0]);
connecBistableType = coder.typeof(connecBistable, [Inf,2], [1,0]);
connecRigidType = coder.typeof(connecRigid, [Inf,2], [1,0]);
connecTorsionType = coder.typeof(connecNodes, [Inf,3], [1,0]);
BC_nodesType = coder.typeof(BC_nodes, [1,Inf], [0,1]);
coords0Type = coder.typeof(coords0, [2,Inf], [0,1]);

codegen -config:mex solve_system.m -args {x0Type,tspan,numVar,mass_VType,eta,eta_iso,k_dome,P0,tau_1,tau_2,tau_3,KsType,Bis_SType,connecNodesType,connecRigidType,connecBistableType,connecTorsionType,BC_nodesType,coords0Type,Nnodes}

% C++ for Gradient and term
codegen -config:mex sym_elasticenergy.m -args {x0sType,KsType,Bis_SType,coords0Type,connecRigidType,connecBistableType,connecTorsionType,BC_nodesType,Nnodes}
codegen -config:mex Grad_sym_elasticenergy.m -args {x0sType,KsType,Bis_SType,coords0Type,connecRigidType,connecBistableType,connecTorsionType,BC_nodesType,Nnodes}
codegen -config:mex damping_F_in.m -args {x0sType,v0Type,eta,coords0Type,connecBistableType,BC_nodesType,Nnodes}

end