%% Function to implicitly solve for Voltage based on a desired input current
% By: Zach Gima 2019-6-26

function x_dot = ode_spmet_cccv(x, cur)
    run param/params_NCA_simulink

    %% Parse Input Data

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

    %%% ZTG Note: Val�en-reimers / CasADi change
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
    V = V_noVCE + V_electrolyteCond + V_electrolytePolar;

    
    %% Aging Dynamics

    %   SEI Layer Growth model
    %   Eqns Adopted from Ramadass et al (2004) [Univ of South Carolina]
    %   "Development of First Principles Capacity Fade Model for Li-Ion Cells"
    %   DOI: 10.1149/1.1634273
    %   NOTE1: This model has NOT been validated experimentally by eCAL
    %   NOTE2: We assume this submodel only applies to anode

    % Difference btw solid and electrolyte overpotential [V]
    phi_se = eta_n + Unref + p.Faraday*R_tot_n*jn_tot;

    % Side exn overpotential [V]
    eta_s = phi_se - p.Us - p.Faraday*R_tot_n * jn_tot;

    % Molar flux of side rxn [mol/s-m^2]
    j_s = -p.i0s/p.Faraday * exp((-p.alph*p.Faraday)/(p.R*T1)*eta_s);

    % SEI layer growth model [m/s]
    delta_sei_dot = -p.M_P/(p.rho_P) * j_s;

    % Molar flux of intercalation
    jn = (abs(jn_tot) - abs(j_s)) * sign(jn_tot);
    jp = jp_tot;

    %% Solid Phase Dynamics

    % ODE for c_s
    c_s_n_dot = A_n*c_s_n + B_n*jn;
    c_s_p_dot = A_p*c_s_p + B_p*jp;


    %% Electrolyte Dynamics

    % Compute Electrolyte Diffusion Coefficient and Derivative
    [D_en0,dD_en0] = electrolyteDe(T1,c_en);
    [D_es0,dD_es0] = electrolyteDe(T1,c_es);
    [D_ep0,dD_ep0] = electrolyteDe(T1,c_ep);

    % Adjustment for Arrhenius temperature dependence
%     Arrh_De = exp(p.E.De/p.R*(1/p.T_ref - 1/T1));
%     D_en = D_en0 * Arrh_De;
%     D_es = D_es0 * Arrh_De;
%     D_ep = D_ep0 * Arrh_De;
%     dD_en = dD_en0 * Arrh_De;
%     dD_es = dD_es0 * Arrh_De;
%     dD_ep = dD_ep0 * Arrh_De;

    %%% Val�en-Reimers
    D_en = p.ElecFactorD * D_en0;
    D_es = p.ElecFactorD * D_es0;
    D_ep = p.ElecFactorD * D_ep0;
    dD_en = p.ElecFactorD * dD_en0;
    dD_es = p.ElecFactorD * dD_es0;
    dD_ep = p.ElecFactorD * dD_ep0;
    
    % Apply BRUGGEMAN RELATION
    D_en_eff = D_en .* p.epsilon_e_n.^(p.brug-1);
    dD_en_eff = dD_en .* p.epsilon_e_n.^(p.brug-1);

    D_es_eff = D_es .* p.epsilon_e_s.^(p.brug-1);
    dD_es_eff = dD_es .* p.epsilon_e_s.^(p.brug-1);

    D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug-1);
    dD_ep_eff = dD_ep .* p.epsilon_e_p.^(p.brug-1);

    % System Matrices have all been precomputed & stored in param struct "p"

    % Compute derivative
    c_en_dot = dD_en_eff.*(p.ce.M1n*c_en + p.ce.M2n*c_e_bcs(1:2)).^2 ...
        + D_en_eff.*(p.ce.M3n*c_en + p.ce.M4n*c_e_bcs(1:2)) + diag(p.ce.M5n)*jn;

    c_es_dot = dD_es_eff.*(p.ce.M1s*c_es + p.ce.M2s*c_e_bcs(2:3)).^2 ...
        + D_es_eff.*(p.ce.M3s*c_es + p.ce.M4s*c_e_bcs(2:3));

    c_ep_dot = dD_ep_eff.*(p.ce.M1p*c_ep + p.ce.M2p*c_e_bcs(3:4)).^2 ...
        + D_ep_eff.*(p.ce.M3p*c_ep + p.ce.M4p*c_e_bcs(3:4)) + diag(p.ce.M5p)*jp;

    % Assemble c_e_dot
    c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];

    %% Thermal Dynamics

    % State-of-Charge (Bulk)
    r_vec = (0:p.delta_r_n:1)';
    c_n = [c_s_n(1); c_s_n; c_ss_n];
    c_p = [c_s_p(1); c_s_p; c_ss_p];
    SOC_n = 3/p.c_s_n_max * trapz(r_vec,r_vec.^2.*c_n);
    SOC_p = 3/p.c_s_p_max * trapz(r_vec,r_vec.^2.*c_p);

    % Equilibrium potentials
    [Unb,~,dUnbdT] = refPotentialAnode(SOC_n);
    [Upb,~,dUpbdT] = refPotentialCathode(SOC_p);

    % Heat generation
    % disp(cur)
    % disp(V)
    % disp((Upb - Unb))
    % pause;
    Qdot = -cur*(V - (Upb - Unb) + T1*(dUpbdT - dUnbdT));

    % Differential equations
    T1_dot = (p.h12 * (T2-T1) + Qdot) / p.C1;
    T2_dot = (p.h12 * (T1-T2) + p.h2a*(p.T_amb - T2)) / p.C2;


    %% Concatenate time derivatives
    x_dot = [c_s_n_dot; c_s_p_dot; c_e_dot; T1_dot; T2_dot; delta_sei_dot];
end