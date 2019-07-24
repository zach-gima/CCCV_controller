%% V_solve function for CCCV Charge Controller
% Repurposing of ode_spmet for implicitly solving for voltage using
% the fsolve function based on a desired input current
% By: Zach Gima 2019-7-12

function V_out = V_solve(V, x, V_ref, Kp, capacity, C_rate,Int,dt)
    run param/params_NCA_simulink
    
    % Compute control law current, but regulate if above max charge rate
%     I_chrgmax = -capacity*C_rate/1000; %A
%     I_ctrl = Kp*(V-V_ref);
%     
%     if I_ctrl < I_chrgmax % charging current is negative, so higher charging rate is a more negative value
%         scale_Kp = I_chrgmax/I_ctrl; % bound the max current; here use max instead of min because charge current is (-)
%         % recompute Kp if         
%         cur = scale_Kp*Kp*(V-V_ref);
%     else % case where charging current from controller doesn't exceed max rate
%         cur = I_ctrl;
%     end

    [cur,~] = current_ctrl_fcn(V,V_ref,Kp,C_rate,capacity,dt,Int);

    % Parse states
    c_s_n = x(1:(p.Nr-1));
    c_s_p = x(p.Nr : 2*(p.Nr-1));
    c_e = x(2*p.Nr-1 : 2*p.Nr-1+p.Nx-4);
    T1 = x(end-2);
    T2 = x(end-1);
    delta_sei = x(end);

    %% Pre-calculations with current states

    %%% MOLAR FLUXES
    % Compute total molar flux
    jn_tot = cur/(p.Faraday*p.a_s_n*p.Area*p.L_n);
    jp_tot = -cur/(p.Faraday*p.a_s_p*p.Area*p.L_p);

    %%% SOLID PHASE DYNAMICS
    % Solid phase diffusivity temperature dependence
    p.D_s_n = p.D_s_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/T1));
    p.D_s_p = p.D_s_n0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/T1));

    % Construct (A,B) matrices for solid-phase Li diffusion
    [A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p] = spm_plant_obs_mats_nocasadi(p);

    % Compute surface concentrations
    c_ss_n = C_n*c_s_n + D_n*jn_tot;
    c_ss_p = C_p*c_s_p + D_p*jp_tot;
    % Remark: I am cheating slightly here. jn_tot should be jn, but doing so
    % imposes an algebraic equation. This forms a DAE. I am going to
    % approximate jn by jn_tot, which should be ok, since jn and jn_tot have a
    % error on the order of 0.001%

    %%% ELECTROLYTE PHASE DYNAMICS
    % Compute electrolyte Boundary Conditions
    c_e_bcs = p.ce.C * c_e;

    ce0n = c_e_bcs(1);
    cens = c_e_bcs(2);
    cesp = c_e_bcs(3);
    ce0p = c_e_bcs(4);

    % Separate and aggregate electrolyte concentration
    c_en = c_e(1:(p.Nxn-1));
    c_es = c_e((p.Nxn-1)+1:(p.Nxn-1)+(p.Nxs-1));
    c_ep = c_e((p.Nxn-1)+p.Nxs : end);
    c_ex = [ce0n; c_en; cens; c_es; cesp; c_ep; ce0p];

    %% Voltage output

    % Average electrolyte concentrations
    cen_bar = mean(c_ex(1:p.Nxn+1,:));
    ces_bar = mean(c_ex((p.Nxn+1):(p.Nxn+p.Nxs+1),:));
    cep_bar = mean(c_ex((p.Nxn+p.Nxs+1):(p.Nxn+p.Nxs+p.Nxp+1),:));

    %%% ZTG Note: Valøen-reimers / CasADi change
    %%% ElecFactorK just nominally set to 1; allows for
    %%% sensitivity to be computed
    %%% No need for Arrhenius dependence since electrolyteCond is
    %%% temp-dependent

    kap_n = p.ElecFactorK*electrolyteCond(T1,cen_bar);
    kap_s = p.ElecFactorK*electrolyteCond(T1,ces_bar);
    kap_p = p.ElecFactorK*electrolyteCond(T1,cep_bar);
    
%     % Overpotentials due to electrolyte subsystem
%     kap_n_ref = electrolyteCond(T1,cen_bar);
%     kap_s_ref = electrolyteCond(T1,ces_bar);
%     kap_p_ref = electrolyteCond(T1,cep_bar);
% 
%     % Adjustment for Arrhenius temperature dependence
%     kap_n = kap_n_ref * exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/T1));
%     kap_s = kap_s_ref * exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/T1));
%     kap_p = kap_p_ref * exp(p.E.kappa_e/p.R*(1/p.T_ref - 1/T1));

    % Bruggeman relationships
    kap_n_eff = kap_n * p.epsilon_e_n.^(p.brug);
    kap_s_eff = kap_s * p.epsilon_e_s.^(p.brug);
    kap_p_eff = kap_p * p.epsilon_e_p.^(p.brug);

    % Activity coefficient
    dfca_n = electrolyteAct(cen_bar,T1,p);
    dfca_s = electrolyteAct(ces_bar,T1,p);
    dfca_p = electrolyteAct(cep_bar,T1,p);

    % Kinetic reaction rate, adjusted for Arrhenius temperature dependence
    p.k_n = p.k_n0 * exp(p.E.kn/p.R*(1/p.T_ref - 1/T1));
    p.k_p = p.k_p0 * exp(p.E.kp/p.R*(1/p.T_ref - 1/T1));

    % Stochiometric Concentration Ratio
    theta_n = c_ss_n / p.c_s_n_max;
    theta_p = c_ss_p / p.c_s_p_max;

    % Equilibrium Potential
    Unref = refPotentialAnode(theta_n);
    Upref = refPotentialCathode(theta_p);

    % Exchange current density
    c_e_bar = [cen_bar; ces_bar; cep_bar];
    [i_0n,i_0p] = exch_cur_dens_nocasadi(p,c_ss_n,c_ss_p,c_e_bar);

    % Overpotentials
    RTaF=(p.R*T1)/(p.alph*p.Faraday);
    eta_n = RTaF * asinh(cur / (2*p.a_s_n*p.Area*p.L_n*i_0n(1)));
    eta_p = RTaF * asinh(-cur / (2*p.a_s_p*p.Area*p.L_p*i_0p(end)));

    % Total resistance (film + growing SEI layer)
    R_tot_n = p.R_f_n + delta_sei/p.kappa_P;
    R_tot_p = p.R_f_p + 0;

    % SPM Voltage (i.e. w/o electrolyte concentration terms)
    V_noVCE = eta_p - eta_n + Upref - Unref ...
        - (R_tot_n/(p.a_s_n*p.L_n*p.Area) + R_tot_p/(p.a_s_p*p.L_p*p.Area))*cur;

    % Overpotential due to electrolyte conductivity
    V_electrolyteCond = (p.L_n/(2*kap_n_eff) + 2*p.L_s/(2*kap_s_eff) + p.L_p/(2*kap_p_eff))*cur;

    % Overpotential due to electrolyte polarization
    V_electrolytePolar = (2*p.R*T1)/(p.Faraday) * (1-p.t_plus)* ...
            ( (1+dfca_n) * (log(cens) - log(ce0n)) ...
             +(1+dfca_s) * (log(cesp) - log(cens)) ...
             +(1+dfca_p) * (log(ce0p) - log(cesp)));

    % Add 'em up!
    V_out = V_noVCE + V_electrolyteCond + V_electrolytePolar - V;
end