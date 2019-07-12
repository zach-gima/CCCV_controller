%%% Re-purposing of ode_spmet_casadi code specifically for generating the
%%% spmet states for use in the Simulink CCCV Charge Controller.

% For Use as an Interpreted Matlab Function (1 Input; 1 Output); can get
% around w/ mux block

% By: Zach Gima 2019-7-11

% INPUTS --- 
% (1) data: struct with fields cur (current), time, V0

% OUTPUTS --- 
% (2) x: 83 (standard discretization params) SPMeT States

function x_sim = spmet_state_sim(x0_nom,Cur,dt)
    import casadi.*
      
    run param/params_NCA_simulink
    
    theta_sx = [];
    theta_0 = [];
    
    %% Initialize ode system setup
    % Define Model Variables
    x = SX.sym('x',size(x0_nom));
    
    % Input
    u = SX.sym('u');
    
    % Build the SPMeT ode System for the intial states; 
    % ode_spmet_casadi holds all of the odes that define the model
    [x_dot, ~, L, ~] = ode_spmet_casadi(x,u,p); % ZTG add 2019-6-17
    
    %% Setup Sensitivity Equations
    % Initial condition
    x0_call = Function('x0_call',{theta_sx},{x0_nom},{'params_sx'},{'result'});
    x0 = full(x0_call(theta_0)); % State I.C
    
    %% Integrator
    % Build ode System
    ode = struct('x', x, 'p', [theta_sx;u], 'ode', x_dot, 'quad', L);
    opts1 = struct('tf',dt, 'abstol',1e-6,'reltol',1e-6,'disable_internal_warnings',1);
    F = integrator('F','cvodes', ode, opts1);

    % Integrate States for Current Time Step 
    Fk = F('x0',x0,'p',[theta_0;Cur]);
    x_sim = full(Fk.xf);
 
end
    
