%% Params for Electrochemical Model
%  Created Apr 26, 2016 by Saehong Park
%
%   Most parameters are from NCA-18650
%   D_e is from Capiglia et al, for c_e = 1000 mol/m^3

% Ref #1
% Uddin et al. (Batteries 2016) (charge degrad.)								
% https://www.mendeley.com/share/document/invite/638008b478/?utm_medium=email&utm_source=transactional&utm_campaign=share%2Finvitation-document								

% Ref #2
% T.R. Ashwin (Power Sources 2016) (mod. El. Par.)
% https://www.mendeley.com/share/document/invite/a3d3f3d0b7/?utm_medium=email&utm_source=transactional&utm_campaign=share%2Finvitation-document

% Ref #3
% Bernadi et al (JPS 2011)
% http://www.sciencedirect.com/science/article/pii/S0378775310011468

%% Geometric Params
% Thickness of each layer
p.L_n = 79e-6;     % Thickness of negative electrode [m]
p.L_s = 80e-6;     % Thickness of separator [m]
p.L_p = 61.5e-6;     % Thickness of positive electrode [m]

% L_ccn = 25e-6;    % Thickness of negative current collector [m]
% L_ccp = 25e-6;    % Thickness of negative current collector [m]

% Particle Radii
p.R_s_n = 10.9e-6;%(True)   % Radius of solid particles in negative electrode [m]
p.R_s_p = 10.9e-6;   % Radius of solid particles in positive electrode [m]

% Volume fractions
p.epsilon_s_n = 0.543889597565723;%0.6987;      % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.666364981170368;%0.5291;      % Volume fraction in solid for pos. electrode

p.epsilon_e_n = 0.3;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 0.5;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.3;   % Volume fraction in electrolyte for pos. electrode

% make element to caclulate phi_{s} by Saehong Park 
p.epsilon_f_n = 1 - p.epsilon_s_n - p.epsilon_e_n;  % Volume fraction of filler in neg. electrode
p.epsilon_f_p = 1 - p.epsilon_s_p - p.epsilon_e_p;  % Volume fraction of filler in pos. electrode

% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]
    
% % % % % % % % % % Mass densities
% % % % % % % % % rho_sn = 1800;    % Solid phase in negative electrode [kg/m^3]
% % % % % % % % % rho_sp = 5010;    % Solid phase in positive electrode [kg/m^3]
% % % % % % % % % rho_e =  1324;    % Electrolyte [kg/m^3]
% % % % % % % % % rho_f = 1800;     % Filler [kg/m^3]
% % % % % % % % % rho_ccn = 8954;   % Current collector in negative electrode
% % % % % % % % % rho_ccp = 2707;   % Current collector in positive electrode
% % % % % % % % % 
% % % % % % % % % % Compute cell mass [kg/m^2]
% % % % % % % % % m_n = p.L_n * (rho_e*p.epsilon_e_n + rho_sn*p.epsilon_s_n + rho_f*epsilon_f_n);
% % % % % % % % % m_s = p.L_s * (rho_e*p.epsilon_e_n);
% % % % % % % % % m_p = p.L_p * (rho_e*p.epsilon_e_p + rho_sp*p.epsilon_s_p + rho_f*epsilon_f_p);
% % % % % % % % % m_cc = rho_ccn*L_ccn + rho_ccp*L_ccp;

% Lumped density [kg/m^2]
% % % % % % % % % p.rho_avg = m_n + m_s + m_p + m_cc;

%% Transport Params
% Diffusion coefficient in solid
p.D_s_n0 = 3.9e-14;%3.9e-14(True)   % Diffusion coeff for solid in neg. electrode, [m^2/s]
p.D_s_p0 = 1e-13; % Diffusion coeff for solid in pos. electrode, [m^2/s]

% Diffusion coefficient for CVM 
p.D_s_n_Factor = 1;
p.D_s_p_Factor = 1;

% Electrolyte Perturbation Factor
p.ElecFactorK = 1; % used for ElectrolyteConductivity
p.ElecFactorD = 1; % used for Electrolyte Diffusivity
p.ElecFactorDA = 1;% used for dactivity

% Diffusional conductivity in electrolyte
p.dactivity = 0;

p.brug = 1.5;       % Bruggeman porosity

% Conductivity of solid
p.sig_n = 100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 100;%(True)    % Conductivity of solid in pos. electrode, [1/Ohms*m]

% Miscellaneous
p.t_plus = 0.38; % 0.36;      % Transference number; use 0.38 for Valoen-Reimers
p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]
p.Area = 0.0744;           % Electrode current collector area [m^2]

%% Kinetic Params
p.R = 8.314472;       % Gas constant, [J/mol-K]

p.alph = 0.5;         % Charge transfer coefficients

p.R_f_n = 5e-4;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 1e-4; % ZTG Change from 1e-3;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_c = 2.5e-03;%5.1874e-05/p.Area; % Contact Resistance/Current Collector Resistance, [Ohms-m^2]
p.R_c = 0; % ZTG Change 2019-4-24

% Nominal Reaction rates
p.k_n0 = 7.5e-4;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
p.k_p0 = 2.3e-3; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]

%% Thermodynamic Params

% Thermal dynamics
p.C_p = 2000;   % Heat capacity, [J/kg-K]
p.h = 0.36;   % Heat transfer coefficient, [W/K-m^2] 0

% Ambient Temperature
p.T_amb = 298.15; % [K]

% Activation Energies
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
% All units are [J/mol]
p.E.kn = 37.48e3;
p.E.kp = 39.57e3;
p.E.Dsn = 42.77e3;
p.E.Dsp = 18.55e3;
p.E.De = 37.04e3;
p.E.kappa_e = 34.70e3;

% Activation Energies
% Energies are from Literature and included Zhang values. Values were
% determined by throwing the highest and lowest values out, then taking the
% mean. The bounds were then selected as the highest/lowest values
% remaining after throwing out the values.

% % For Ea-kn and Ea-kp, only 3 values were chosen from so the middle
% % value was used and the high/low values were used as bounds.
% % All units are [J/mol]
p.E.Dsn = 36.63e3;
p.E.Dsp = 47.98e3;
p.E.kn = 53.4e3;
p.E.kp = 39.57e3;

% Reference temperature
p.T_ref = 298.15; %[K]

% % % Heat transfer parameters
% % % Taken from Zhang et al (2014) [Harbin]
% % % http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
% % p.C1 = 62.7;    % [J/K]
% % p.C2 = 4.5;     % [J/K]
% % p.h12 = 10; % [W/K]
% % p.h2a = 21.45;  % [W/K]

% Heat transfer parameters
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
p.C1 = 62.7/p.Area;    % C_c [J/K]
p.C2 = 4.5/p.Area;     % C_s [J/K]
p.h12 = (1/10)/p.Area; % R_c [W/K]
p.h2a = (1/21.45)/p.Area;  % R_u [W/K]

% % % Heat transfer parameters
% % % Taken from Hector's disseratation
% % p.C1 = 62.7/p.Area;    % [J/(m^2-K)]
% % p.C2 = 4.5/p.Area;     % [J/(m^2-K)]
% % p.h12 = (1/1.94)/p.Area; % [W/(m^2-K)]
% % p.h2a = 1.5*(1/3.08)/p.Area;  % [W/(m^2-K)]

%%% Parameters for side reaction -- SEI layer growth [ZTG Change]
%%% Note: the values may be specific to an LCO chemistry (pulled from spmet
%%% code)
p.kappa_P = 1;      % [S/m] conductivity of side rxn product
p.M_P = 7.3e1;      % [kg/mol] molecular weight of side rxn product
p.rho_P = 2.1e3;    % [kg/m^3] mass density of side rxn product
p.i0s = 0; %1.5e-6;     % [A/m^2] exchange current density of side rxn
p.Us = 0.4;         % [V] reference potential of side rxn

%% Concentrations
% Maxima based on DUALFOIL 
% line 588 in DUALFOIL Fortran code

p.c_s_n_max = 3.71e+04;   % Max concentration in anode, [mol/m^3]

%p.c_s_n_max = 3.6e3 * 372 * 1800 / p.Faraday;   % Max concentration in anode, [mol/m^3]


p.c_s_p_max = 5.10e+04;    % Max concentration in cathode, [mol/m^3]

%p.c_s_p_max = 3.6e3 * 274 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]

p.n_Li_s = 0.1406; %2.5; %2.781;        % Total moles of lithium in solid phase [mol]
p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]

%% Cutoff voltages
p.volt_max = 4.2; %4.1113; %4.7;
p.volt_min = 2.5; %2.6;

%% SPMeT Specific Finite Difference Scheme
p.Nr = 30; % 100 Make this very large so it closely approximates the true model
p.delta_r_n = 1/p.Nr;
p.delta_r_p = 1/p.Nr;
p.r_vec = (0:p.delta_r_n:1)';
%     r_vecx = r_vec(2:end-1);

% Finite difference points along x-coordinate
p.Nxn = 10; %70; 
p.Nxs = 5; %35; 
p.Nxp = 10; %70;
p.Nx = p.Nxn+p.Nxs+p.Nxp;
%     Nx = p.Nx - 3;
%     x_vec_spme = linspace(0,1,Nx+4);

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;p.delta_x_p = 1 / p.Nxp;

%% CE Mats stuff

%%% Lumped Coefficients
Del_xn = p.L_n * p.delta_x_n;
Del_xs = p.L_s * p.delta_x_s;
Del_xp = p.L_p * p.delta_x_p;

%%% Matrices in nonlinear dynamics
M1n = sparse((diag(ones(p.Nxn-2,1),+1) - diag(ones(p.Nxn-2,1),-1))/(2*Del_xn));
M1s = sparse((diag(ones(p.Nxs-2,1),+1) - diag(ones(p.Nxs-2,1),-1))/(2*Del_xs));
M1p = sparse((diag(ones(p.Nxp-2,1),+1) - diag(ones(p.Nxp-2,1),-1))/(2*Del_xp));


M2n = zeros(p.Nxn-1,2);
M2n(1,1) = -1/(2*Del_xn);
M2n(end,end) = 1/(2*Del_xn);
M2n = sparse(M2n);

M2s = zeros(p.Nxs-1,2);
M2s(1,1) = -1/(2*Del_xs);
M2s(end,end) = 1/(2*Del_xs);
M2s = sparse(M2s);

M2p = zeros(p.Nxp-1,2);
M2p(1,1) = -1/(2*Del_xp);
M2p(end,end) = 1/(2*Del_xp);
M2p = sparse(M2p);


M3n = sparse((-2*diag(ones(p.Nxn-1,1),0) + diag(ones(p.Nxn-2,1),+1) + diag(ones(p.Nxn-2,1),-1))/(Del_xn^2));
M3s = sparse((-2*diag(ones(p.Nxs-1,1),0) + diag(ones(p.Nxs-2,1),+1) + diag(ones(p.Nxs-2,1),-1))/(Del_xs^2));
M3p = sparse((-2*diag(ones(p.Nxp-1,1),0) + diag(ones(p.Nxp-2,1),+1) + diag(ones(p.Nxp-2,1),-1))/(Del_xp^2));


M4n = zeros(p.Nxn-1,2);
M4n(1,1) = 1/(Del_xn^2);
M4n(end,end) = 1/(Del_xn^2);
M4n = sparse(M4n);

M4s = zeros(p.Nxs-1,2);
M4s(1,1) = 1/(Del_xs^2);
M4s(end,end) = 1/(Del_xs^2);
M4s = sparse(M4s);

M4p = zeros(p.Nxp-1,2);
M4p(1,1) = 1/(Del_xp^2);
M4p(end,end) = 1/(Del_xp^2);
M4p = sparse(M4p);

M5n = (1-p.t_plus)*p.a_s_n/p.epsilon_e_n * speye(p.Nxn-1);
M5p = (1-p.t_plus)*p.a_s_p/p.epsilon_e_p * speye(p.Nxp-1);

%%% Boundary Conditions
N1 = zeros(4,p.Nx-3); 
N2 = zeros(4);

% BC1
N1(1,1) = +4;
N1(1,2) = -1;
N2(1,1) = -3;

% BC2
N1(2,p.Nxn-2) = (p.epsilon_e_n^p.brug)/(2*Del_xn);
N1(2,p.Nxn-1) = (-4*p.epsilon_e_n^p.brug)/(2*Del_xn);
N2(2,2) = (3*p.epsilon_e_n^p.brug)/(2*Del_xn) + (3*p.epsilon_e_s^p.brug)/(2*Del_xs);
N1(2,p.Nxn) = (-4*p.epsilon_e_s^p.brug)/(2*Del_xs);
N1(2,p.Nxn+1) = (p.epsilon_e_s^p.brug)/(2*Del_xs);

% BC3
N1(3,p.Nxn+p.Nxs-3) = (p.epsilon_e_s^p.brug)/(2*Del_xs);
N1(3,p.Nxn+p.Nxs-2) = (-4*p.epsilon_e_s^p.brug)/(2*Del_xs);
N2(3,3) = (3*p.epsilon_e_s^p.brug)/(2*Del_xs) + (3*p.epsilon_e_p^p.brug)/(2*Del_xp);
N1(3,p.Nxn+p.Nxs-1) = (-4*p.epsilon_e_p^p.brug)/(2*Del_xp);
N1(3,p.Nxn+p.Nxs) = (p.epsilon_e_p^p.brug)/(2*Del_xp);


% BC4
N1(4,end-1) = +1;
N1(4,end) = -4;
N2(4,4) = +3;

C_ce = sparse(-N2\N1); 

p.ce.M1n = M1n;
p.ce.M2n = M2n;
p.ce.M3n = M3n;
p.ce.M4n = M4n;
p.ce.M5n = M5n;

p.ce.M1s = M1s;
p.ce.M2s = M2s;
p.ce.M3s = M3s;
p.ce.M4s = M4s;

p.ce.M1p = M1p;
p.ce.M2p = M2p;
p.ce.M3p = M3p;
p.ce.M4p = M4p;
p.ce.M5p = M5p;

p.ce.C = C_ce;

clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce;
