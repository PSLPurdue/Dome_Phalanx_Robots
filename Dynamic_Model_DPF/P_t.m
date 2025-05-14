function P_t = P_t(P0, t, tau_1, tau_2, tau_3, tau_4)
% P_T  Time-dependent pressure profile with ramp-up, hold, and ramp-down
%   Computes pressure P_t based on piecewise linear segments defined by
%   time thresholds tau_1 to tau_4 and peak pressure P0.
%
%   Inputs:
%     P0    – peak pressure magnitude
%     t     – current time
%     tau_1 – time when pressure starts ramping up
%     tau_2 – time when pressure reaches P0
%     tau_3 – time when pressure begins ramping down
%     tau_4 – time when pressure returns to zero
%
%   Output:
%     P_t   – scalar pressure at time t

% 1. Before ramp-up: pressure is zero
if t <= tau_1
    P_t = 0;

% 2. Ramp-up: linearly increase from 0 to P0 between tau_1 and tau_2
elseif t <= tau_2
    P_t = P0 * (t - tau_1) / (tau_2 - tau_1);

% 3. Hold: maintain peak pressure P0 from tau_2 to tau_3
elseif t <= tau_3
    P_t = P0;

% 4. Ramp-down: linearly decrease from P0 back to 0 between tau_3 and tau_4
elseif t <= tau_4
    P_t = P0 * (tau_4 - t) / (tau_4 - tau_3);

% 5. After ramp-down: pressure returns to zero
else
    P_t = 0;
end

end
