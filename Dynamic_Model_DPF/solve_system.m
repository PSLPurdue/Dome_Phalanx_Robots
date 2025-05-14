function [time, x] = solve_system(x0, tspan, numVar, mass_V, eta, eta_iso, k_dome, P0, tau_1, tau_2, tau_3, tau_4, Ks, Bis_S, connecNodes, connecRigid, connecBistable, connecTorsion, BC_nodes, coords0, Nnodes)
% SOLVE_SYSTEM  Solve dynamic ODE system using MATLAB's stiff solver
%   Sets up and integrates the coupled equations of motion defined in 'rhs'.
%
%   Inputs:
%     x0             – initial state [displacements; velocities]
%     tspan          – [t0 tf] integration interval
%     numVar         – number of displacement DOFs
%     mass_V         – vector of nodal masses (size numVar×1)
%     eta, eta_iso   – damping coefficients
%     k_dome         – dome stiffness multiplier
%     P0, tau_1-4    – parameters for time-dependent pressure profile
%     Ks, Bis_S      – spring constants and bistable parameters
%     connecNodes    – connectivity for damping springs
%     connecRigid    – rigid spring connectivity pairs
%     connecBistable – bistable spring connectivity pairs
%     connecTorsion  – torsional spring triplets
%     BC_nodes       – fixed node indices
%     coords0        – undeformed node coordinates (2×Nnodes)
%     Nnodes         – total number of nodes
%
%   Outputs:
%     time           – vector of time points returned by ode15s
%     x              – solution array [time steps×2*numVar]

%% 1. Create ODE function handle
% Capture all parameters needed by 'rhs'
ode_fun = @(t, x) rhs(t, x, numVar, mass_V, eta, eta_iso, ...
    k_dome, P0, tau_1, tau_2, tau_3, tau_4, Ks, Bis_S, ...
    connecNodes, connecRigid, connecBistable, connecTorsion, ...
    BC_nodes, coords0, Nnodes);

%% 2. Integrate using ode15s (suitable for stiff dynamics)
[time, x] = ode15s(ode_fun, tspan, x0);

end