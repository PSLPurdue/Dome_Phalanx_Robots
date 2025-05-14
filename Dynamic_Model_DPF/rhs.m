function dx_dt = rhs(t, x, numVar, mass_V, eta, eta_iso, k_dome, P0, tau_1, tau_2, tau_3, tau_4, Ks, Bis_S, connecNodes, connecRigid, connecBistable, connecTorsion, BC_nodes, coords0, Nnodes)
% RHS  Right-hand side of ODE system for dynamic simulation
%   Assembles velocities and accelerations from external, internal,
%   damping, and inertial forces to form state derivatives.
%
%   Inputs:
%     t               – current time
%     x               – state vector [displacements; velocities]
%     numVar          – number of displacement DOFs
%     mass_V          – vector of nodal masses (size numVar×1)
%     eta             – damping coefficient for viscous damping
%     eta_iso         – isotropic damping coefficient
%     k_dome          – stiffness multiplier for dome forces
%     P0, tau_1-4     – parameters for time-dependent pressure P_t
%     Ks, Bis_S       – spring constants and bistable params
%     connecNodes     – connectivity for damping (bistable) springs
%     connecRigid     – rigid & limit-layer spring connectivity
%     connecBistable  – bistable spring connectivity
%     connecTorsion   – torsional spring connectivity
%     BC_nodes        – indices of fixed (boundary) nodes
%     coords0         – original undeformed node coordinates
%     Nnodes          – total number of nodes
%
%   Output:
%     dx_dt           – time derivative [velocities; accelerations]

%% 1. Split state vector into displacement and velocity
x_1 = x(1:numVar);           % current displacements
x_2 = x(numVar+1:2*numVar);  % current velocities

%% 2. Compute external force vector at time t
F_t = F_ext(t, x_1, P0, tau_1, tau_2, tau_3, tau_4, k_dome, ...
    connecRigid, connecBistable, BC_nodes, coords0, Nnodes);

%% 3. Compute internal elastic forces (negative gradient of energy)
% Grad_sym_elasticenergy returns row vector; transpose to column
F_in = Grad_sym_elasticenergy(x_1, Ks, Bis_S, coords0, ...
    connecRigid, connecBistable, connecTorsion, BC_nodes, Nnodes)';

%% 4. Compute damping forces
% 4.1 Viscous damping from bistable springs
F_in_vis = damping_F_in(x_1, x_2, eta, coords0, connecNodes, BC_nodes, Nnodes);
% 4.2 Isotropic (Rayleigh) damping proportional to velocity
F_iso    = eta_iso * x_2;

%% 5. Formulate equations of motion
% 5.1 Velocities derivative is current velocities
xp_1 = x_2;
% 5.2 Accelerations from net force divided by nodal mass
%       net force = external - internal - viscous - isotropic
dxp = (F_t - F_in - F_in_vis - F_iso) ./ mass_V;

%% 6. Return time derivative of state
dx_dt = [xp_1; dxp];
end
