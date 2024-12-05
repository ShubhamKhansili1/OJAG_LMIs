% This code finds an upper bound on Delta such that the LMIs are feasible.
a=0.1; % Diffusion coefficient.
K=1; % Global control gain.
k=1; % First Local control gain giving rise to dirichlet type boundary.
sigma=0.8; % Second local control gain giving rise to Wentzell boundary.
delta= 0.8;  % Decay rate.
L= 0.001;  % Lipschitz Constant.

% Range of Delta values to check
Delta_values = 0:0.01:1 ; 

%% Wentzell | No communication.
for i = 1:length(Delta_values)
    Delta_W_NC = Delta_values(i);
    
    % Call the feasibility function
    feas = OJAG_Wentzell_No_Com(a,K,k,sigma,L,delta,Delta_W_NC);
    
    if feas ~= 1
        fprintf('LMI is not feasible for Delta with Wentzell BC (consant case) = %.4f.\n', Delta_W_NC);
        break; % Stop the loop
    end
end
%% Wentzell | with Communication.
for i = 1:length(Delta_values)
    Delta_W_C = Delta_values(i);
    
    % Call the feasibility function
    feas = OJAG_Wentzell_Com(a,K,k,sigma,L,delta, Delta_W_C);
    
    if feas ~= 1
        fprintf('LMI is not feasible for Delta with Wentzell BC (linear case) = %.4f.\n', Delta_W_C);
        break; % Stop the loop
    end
end

