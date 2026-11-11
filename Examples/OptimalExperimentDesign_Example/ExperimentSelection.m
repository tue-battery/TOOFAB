% Select optimal experiments
% 
% Inputs:
%    - equil is a struct that contains the information on the equilbrium 
%      potential of the cell (or the EMF). When the EMF is used, there should be 
%      2 fields in equil: EMF and Cbat. EMF should be function of SOC and Cbat
%      should be given in C. 
%    - options is a struct where you can specify additional options for the experiment
%      selection procedure. These options contain; the amount of experiment
%      that are generated (expamount), the maximum time of the optimal
%      experiments (exptime) in hours, set this to inf is no limit is
%      wanted, and the maximum cell voltage (cellvoltage) in volts.
%      The settings of the DFN model that is used in the simulations can be
%      changed in options.dfn (e.g., options.dfn.T_amb). The defaults are
%      the same as in parameter_determination.m. The parameterization of the
%      model (parameter ranges and scaling) is also the same as in
%      parameter_determination.m, so keep the two files consistent.
%
% Outputs:
%    - Currents and Currents_cut are all the generated experiments and the
%      experiments cut off at the battery's capacity, respectively.
%    - Sensitivitymatrices is a cell vector containing all the sensitivity
%      analyses of the cut off generated experiments
%    - OptimalExperimentA, OptimalExperimentD, and OptimalExperimentE are
%      the optimal experiments found using hte A-optimality, D-optimality,
%      and E-optimality as a metric of optimality.

function optimalexperiments = ExperimentSelection(equil,options_in)
options.experiment_amount = 10;
options.experiment_time_max = 8.75;
options.cell_voltage = 4.18;
options.safety_limit = [-8, 8];
options.number_of_pulses = 50;
options.debug = 0;

% Settings of the DFN model in the simulations. These are the same defaults
% as in parameter_determination.m, such that the experiments are selected for
% the same model as used in the parameter estimation. Specific fields can be
% changed with options.dfn (e.g., options.dfn.T_amb or options.dfn.set_simp).
options.dfn.thermal_dynamics = 0;
options.dfn.set_simp = [2 2 2 2 1 0];
options.dfn.T_amb = 273.15+26;
options.dfn.verbose = 0;
options.dfn.set_grid = [10 10 10 10 10];

if nargin>1
    if isfield(options_in,'dfn')
        f = fieldnames(options_in.dfn);
        for i = 1:length(f)
            options.dfn.(f{i}) = options_in.dfn.(f{i});
        end
        options_in = rmfield(options_in,'dfn');
    end
    f = fieldnames(options_in);
    for i = 1:length(f)
        options.(f{i}) = options_in.(f{i});
    end
end

N = options.experiment_amount; % Number of experiments


opt.no_of_pulses = options.number_of_pulses;
opt.upper_safety_limit = options.safety_limit(2);
opt.lower_safety_limit = options.safety_limit(1);
L = [];
disp("Generating Experiments")
for i = 1:N
    CC = 0; t_c = 0; desired_cbat = equil.Cbat;
    while CC(end)<desired_cbat || t_c>options.experiment_time_max*3600
        % The next 8 lines are the dials to change the shape of the
        % experiments, change the numbers in front of rand(n) 
        opt.std_pulse_period = 2*abs(randn);    % standard deviation of pulse period
        opt.std_rest_period = 1*abs(randn);     % standard deviation of rest period
        opt.std_random_amplitude = 1/20*abs(randn); % stochastic part of the signal
        opt.mean_pulse_period = 5e2*rand;   % mean of pulse period
        opt.mean_rest_period = 5e2*rand;    % mean of rest period
        opt.amplitude_plus = 4*rand;       % amplitude deviation to positive side
        opt.amplitude_minus = 4*rand;      % amplitude deviation to negative side
        opt.dc_offset = -2;    % DC offset value
        Current{i} = generate_current_profile(opt.mean_pulse_period, opt.std_pulse_period, opt.mean_rest_period, opt.std_rest_period, opt.amplitude_plus, opt.amplitude_minus, opt.std_random_amplitude, opt.dc_offset, opt.no_of_pulses, opt.upper_safety_limit, opt.lower_safety_limit);
        CC = -cumtrapz(Current{i});
        [~, t_c] = min(abs(desired_cbat-CC));
    end
    Current_cut{i} = Current{i}(1:t_c);
    % L_ = length(Current_cut{1});
    L = [L, length(Current_cut{i})];
    disp([i+"/"+N])
end
%% Sensitivity matrices
disp("Starting Sensitivity Analysis")
perturbation = 0.1;
beta_0 = 0.5;
par_names = parameter_names();
N_par = length(par_names);
% The sensitivities are computed for the parameters that influence the
% voltage in an isothermal simulation. The temperature dependencies
% (activation energies), the thermal parameters, the activity coefficient
% derivative dlnfdce (zero range) and the current collector resistance R_cc
% are excluded.
excluded_par = {'Ds_neg_Ea','Ds_pos_Ea','knEa','kpEa','DEa','dlnfdce','R_cc','hc','Cp','V_cell'};
sens_idx = find(~ismember(par_names,excluded_par));
N_par_sens = length(sens_idx);
for i = 1:N
    disp(['Experiment: ',num2str(i),'/',num2str(N)]);
    data.I = Current_cut{i}; data.t = 1:length(data.I); data.V = options.cell_voltage;
    CC = -cumtrapz(data.I); equil.Cbat = CC(end);
    beta_nom = beta_0*ones(1,N_par); %nominal normalized parameters
    y_nom = dfn_model(equil,data,beta_nom,options.dfn); %Obtain the model output with the nominal parameters
    % disp(['y_nom simulated, exp: ',num2str(i),'/',num2str(N)])
    y_nom(isnan(y_nom)) = 2.5;
        parfor k = 1:N_par_sens
            beta_perm_low = beta_nom;
            beta_perm_low(sens_idx(k)) = (beta_0-perturbation);
            y_perm_low{k} = dfn_model(equil,data,beta_perm_low,options.dfn); %Obtain the model output with the lower permutation of the parameters
            y_perm_low{k}(isnan(y_perm_low{k})) = 2.5;
            dydpar_low(k,:) = -(y_perm_low{k}-y_nom)/perturbation; %Compute the sensitivity of the model output to the parameters based on the difference between the lower permutation of the parameters and the nominal parameters
    
            beta_perm_high = beta_nom;
            beta_perm_high(sens_idx(k)) = (beta_0+perturbation);
            y_perm_high{k} = dfn_model(equil,data,beta_perm_high,options.dfn); %Obtain the model output with the higher permutation of the parameters
            y_perm_high{k}(isnan(y_perm_high{k})) = 2.5;
            dydpar_high(k,:) = (y_perm_high{k}-y_nom)/perturbation; %Compute the sensitivity of the model output to the parameters based on the difference between the higher permutation of the parameters and the nominal parameters
            % disp(['Sens ',num2str(k),' has finished'])
        end
    S{i} = (0.5*(abs(dydpar_low)+abs(dydpar_high)))'; %Determine the sensitivity of the model output to the parameters are the average of the two computed sensitivities
    % disp(['S matrix constructed, exp: ',num2str(i),'/',num2str(N)])
    clear y_perm_low y_perm_high dydpar_low dydpar_high
end
% save Generated_Currents_and_Smatrices.mat S Current Current_cut std mean_period mean_amplitude DC N -v7.3
optimalexperiments.Currents = Current;
optimalexperiments.Currents_cut = Current_cut;
optimalexperiments.SensitivityMatrices = S;
%% Fisher information matrix
disp("Starting Experiment Selection")
for i = 1:N
    Sc = S{i};
    F{i} = Sc'*Sc;
    A_opt_inv(i) = trace(pinv(F{i}));
    D_opt_inv(i) = log(det(pinv(F{i})));
    E_opt(i) = max(eig(pinv(F{i})));
    % disp(['Fisher matrix constructed, exp: ',num2str(i),'/',num2str(N)])
end

%% A optimality
[A_sort, A_indx] = sort(A_opt_inv,'ascend');

%% D optimality
[D_sort, D_indx] = sort(D_opt_inv,'ascend');

%% E optimality
[E_sort, E_indx] = sort(E_opt,'ascend');

%% Experiment selection
for i = 1:N
    L_A(i) = length(Current_cut{A_indx(i)})/3600;
    L_D(i) = length(Current_cut{D_indx(i)})/3600;
    L_E(i) = length(Current_cut{E_indx(i)})/3600;
end
aopt = A_indx(find(L_A<=options.experiment_time_max,1));

dopt = D_indx(find(L_D<=options.experiment_time_max,1));

eopt = E_indx(find(L_E<=options.experiment_time_max,1));

optimalexperiments.OptimalExperimentA = Current_cut{aopt};
optimalexperiments.OptimalExperimentD = Current_cut{dopt};
optimalexperiments.OptimalExperimentE = Current_cut{eopt};
if options.debug
    optimalexperiments.debug = {Current,Current_cut,A_sort,D_sort,E_sort,L_A,L_D,L_E,A_indx,D_indx,E_indx,aopt,dopt,eopt};
end

disp("Finished Experiment Selection")

end

%% Functions
function pulse_train = generate_rectangular_pulse_train(pulse_periods)
    % Generate a rectangular pulse train.
    
    % Time instants (with sampling period = 1 sec)
    time = 0:sum(pulse_periods)-1;
    
    % Starting time instants of pulses
    pulse_start_times = cumsum(pulse_periods(1:end-1));
    pulse_start_times = [0, pulse_start_times];
    
    % Generate rectangular shape by setting 1s' between
    % starting times and half of the corresponding period.
    % Later, it is shifted downward by 0.5, and then multiplied
    % with 2 to scale the amplitudes to unity. The minus ensures
    % the discharge pulse phase at the beginning.
    rectangular_shape = (time >= pulse_start_times' & time < (pulse_start_times' + 0.5 * pulse_periods'));
    pulse_train = -2 * (sum(rectangular_shape, 1) - 0.5);
end

function identification_signal = generate_current_profile(mean_pulse_period, ...
    std_pulse_period, mean_rest_period, std_rest_period, ...
    amplitude_plus, amplitude_minus, std_random_amplitude, dc_offset, ...
    no_of_pulses, upper_safety_limit, lower_safety_limit)

    % Generate periods for each pulse
    rng('shuffle'); % Initialize random generator
    pulse_periods = round(abs(randn(1, no_of_pulses) * std_pulse_period + mean_pulse_period));
    
    % Generate pulse train
    pulse_train = generate_rectangular_pulse_train(pulse_periods);

    % Split pulse train into list of phases
    phase_start_times = ceil(cumsum(repmat(pulse_periods, 1, 2) / 2));
    phases = mat2cell(pulse_train, 1, diff([0 phase_start_times]));
    
    % Generate rest periods for each phase
    rest_periods = round(abs(randn(1, numel(phases)) * std_rest_period + mean_rest_period));
    
    % Final amplitude for each phase
    final_amplitudes = repmat(rand(1, ceil(numel(phases)/2)) * (amplitude_plus - amplitude_minus) + amplitude_minus, 1, 2);
    
    identification_signal = [];
    
    for i = 1:numel(phases)
        phase = phases{i};
        
        % Add rest signal (zeros)
        rest_signal = zeros(1, rest_periods(i));
        
        % Add high-frequency noise
        zero_mean_random_signal = randn(1, numel(phase)) * std_random_amplitude;
        
        % Finalize the phase signal
        finalized_phase = final_amplitudes(i) * (phase + dc_offset + zero_mean_random_signal);
        
        % Apply safety limits
        if ~isempty(upper_safety_limit)
            finalized_phase(finalized_phase > upper_safety_limit) = upper_safety_limit;
        end
        if ~isempty(lower_safety_limit)
            finalized_phase(finalized_phase < lower_safety_limit) = lower_safety_limit;
        end
        
        % Concatenate signal with rest period
        identification_signal = [identification_signal, finalized_phase, rest_signal];
    end
end

function y = dfn_model(equil,data,beta,dfn_options) %NB
% Simulate the terminal voltage of the DFN model for the normalized
% parameters beta. The mapping from beta to the model parameters is the same
% as in parameter_determination.m, such that the sensitivities are computed
% for the same parameterization as used in the parameter estimation.
options.equil_mode = 1;
options.U_neg = @(theta_n) 0.7222 + 0.1387*theta_n + 0.029*theta_n.^0.5 - 0.0172./theta_n + 0.0019./...
        theta_n.^1.5 + 0.2808*exp(0.9-15*theta_n)-0.7984*exp(0.4465*theta_n - 0.4108);

par_names = parameter_names();
ranges = determine_parameter_ranges(equil.Cbat);
methods = parameter_methods(par_names);
ph = beta2par(beta,par_names,ranges,methods);

p = phat2p(ph,equil.Cbat);
equil = determine_equilibrium_potentials(equil,ph,options);

p.U_pos = equil.U_pos;
p.U_neg = equil.U_neg;
p.Cbat = equil.Cbat;
p.s0_neg = ph.s0_neg;
p.s0_pos = ph.s0_pos;
p.s100_neg = ph.s100_neg;
p.s100_pos = ph.s100_pos;

% DFN simulation settings (thermal model, ambient temperature, grid, ...)
f = fieldnames(dfn_options);
for i = 1:length(f)
    p.(f{i}) = dfn_options.(f{i});
end

p.ageing = 0;
if isfield(data,'T')
    out = DFN([data.t' data.I' data.T'],data.t(end),data.V(1),p);
else
    out = DFN([data.t' data.I'],data.t(end),data.V(1),p);
end

if length(out.V)>1
    V = interp1(out.t,out.V,data.t);
    % If the model time step is larger than the sample time of the data, the
    % first samples of the data lie before the first model output. These are
    % set to the first model output instead of NaN. A simulation that stops
    % early is still returned as NaN.
    V(data.t < out.t(1)) = out.V(1);
else
    V = out.V;
end
y = V;
end

function par_names = parameter_names()
% Names of the (grouped) parameters, in the same order as in
% parameter_determination.m (for the equilibrium potential mode with the EMF)
par_names = {'s0_neg','s0_pos','s100_neg','s100_pos','Ds_neg', 'Ds_pos', 'De', 'p_neg', 'p_sep', 'p_pos', 't_plus',...
        'sigma_neg', 'sigma_pos','kappa', 'Rf_neg', 'Ds_neg_Ea','Ds_pos_Ea', 'knEa','kpEa','DEa'...
        'alpha', 'k0_neg', 'k0_pos', 'epse_neg', 'epse_sep', 'epse_pos',...
        'dlnfdce','R_cc','hc', 'Cp', 'V_cell'};
end

function methods = parameter_methods(par_names)
% Method that defines how the parameters are determined from beta, as in
% parameter_determination.m (0: linear, 1: logarithmic). All parameters use
% linear scaling between their minimum and maximum.
methods = zeros(1,length(par_names));
end

function p = phat2p(ph,Cbat)
p.F = 96487;
% Free choice parameters
rho = 20; %rho is the capacity per m^2. Taken as a typical value based on Arora et al, 2000 (rho = 23), Ecker et al, 2015 (rho = 15), Sturm et al, 2019 (rho = 26)
p.A_surf = (Cbat/3600)/rho;
p.ce0 = 1000;
p.delta_neg =  8.5e-5;
p.delta_pos = 9.2e-5;
p.delta_sep = 2e-5;
p.L = p.delta_sep+p.delta_neg+p.delta_pos;
p.epss_neg = 0.53;
p.epss_pos = 0.51;
p.R_neg = 8.1e-6;
p.R_pos = 5.8e-6;


p.Ds_neg = @(T) (ph.Ds_neg*exp((ph.Ds_neg_Ea/(8.3145))*((1/298)-(1/T))))*p.R_neg^2;
p.Ds_pos = @(T) (ph.Ds_pos*exp((ph.Ds_pos_Ea/(8.3145))*((1/298)-(1/T))))*p.R_pos^2;

p.t_plus = ph.t_plus;
p.dlnfdce = ph.dlnfdce;


p.De = @(T) (ph.De*exp((ph.DEa/(8.3145))*((1/298)-(1/T))))*(1-p.t_plus)/(p.F*p.A_surf*p.ce0);

p.epse_neg = ph.epse_neg*(1-p.t_plus)/(p.F*p.A_surf*p.delta_neg*p.ce0);
p.epse_sep = ph.epse_sep*(1-p.t_plus)/(p.F*p.A_surf*p.delta_sep*p.ce0);
p.epse_pos = ph.epse_pos*(1-p.t_plus)/(p.F*p.A_surf*p.delta_pos*p.ce0);

p.p_neg = log(ph.p_neg*p.delta_neg)/log(p.epse_neg);
p.p_sep = log(ph.p_sep*p.delta_sep)/log(p.epse_sep);
p.p_pos = log(ph.p_pos*p.delta_pos)/log(p.epse_pos);

p.sigma_neg = ph.sigma_neg*p.delta_neg/(p.epss_neg*p.A_surf);
p.sigma_pos = ph.sigma_pos*p.delta_pos/(p.epss_pos*p.A_surf);

p.kappa = ph.kappa/p.A_surf;

p.Rf_neg = ph.Rf_neg*3*p.A_surf*p.delta_neg*p.epss_neg/p.R_neg;
p.Rf_pos = 0;
p.R_cc = ph.R_cc*p.A_surf;

p.alpha_a = ph.alpha;
p.alpha_c = 1-p.alpha_a;

p.k0_neg = @(T) (ph.k0_neg*exp((ph.knEa/(8.3145))*((1/298)-(1/T))))*p.R_neg*p.F/p.ce0^p.alpha_a;
p.k0_pos = @(T) (ph.k0_pos*exp((ph.kpEa/(8.3145))*((1/298)-(1/T))))*p.R_pos*p.F/p.ce0^p.alpha_a;

p.V_cell = ph.V_cell;
% p.m = 1;
p.hc = ph.hc;
p.Cp = ph.Cp;
end

function par = beta2par(beta,par_names,ranges,method)
for i = 1:length(beta)
    if method(i)==0
        par.(par_names{i}) = (1-beta(i))*ranges.([par_names{i},'_min'])+beta(i)*ranges.([par_names{i},'_max']);
    else
        par.(par_names{i}) = 10.^((1-beta(i))*log10(ranges.([par_names{i},'_min']))+beta(i)*log10(ranges.([par_names{i},'_max'])));
    end
end
end

function ph = determine_parameter_ranges(Cbat)
%% Parameter ranges
%We first define the parmeters in their ungrouped representation - to allow
%users to easily compare with the literature. We then convert the parameter
%ranges to the grouped parameter
% Original parameter ranges
p.dlnfdce_min = 0;
p.dlnfdce_max = 0;

p.s0_neg_min = 0.002;
p.s0_neg_max = 0.13;
p.s0_pos_min = 0.75;
p.s0_pos_max = 0.97;

p.s100_neg_min = 0.5;
p.s100_neg_max = 0.89;

p.s100_pos_min = 0.1;
p.s100_pos_max = 0.5;

p.Ds_neg_min = 2e-18; 
p.Ds_neg_max = 1e-6;

p.Ds_pos_min = 4.123e-16;
p.Ds_pos_max = 1e-9;


p.Ds_neg_Ea_min=1000;
p.Ds_neg_Ea_max=100000;
p.Ds_pos_Ea_min=1000;
p.Ds_pos_Ea_max=100000;

p.DEa_min=1000;
p.DEa_max=100000;

p.epse_neg_min = 0.0043;
p.epse_neg_max = 0.42;

p.epse_sep_min = 0.0740;
p.epse_sep_max = 0.99;

p.epse_pos_min = 0.0053;
p.epse_pos_max = 0.48;

p.t_plus_min = 0.26;
p.t_plus_max = 0.4;


p.sigma_neg_min = 3;
p.sigma_neg_max = 1e5;

p.sigma_pos_min = 0.00352;
p.sigma_pos_max = 72.5;

%The solid volume fractions are fixed (same values as in phat2p), such that
%the grouped ranges of sigma match the conversion in phat2p
p.epss_neg = 0.53;
p.epss_pos = 0.51;

p.kappa_min = 7e-6*Cbat;%
p.kappa_max = 2e-5*Cbat;%

p.Rf_neg_min = 20/Cbat;%
p.Rf_neg_max = 330/Cbat;%

p.R_cc_min = 2e-4;
p.R_cc_max = 2e-2;

p.alpa_min = 0.3;
p.alpa_max = 0.7;

p.knEa_min=1000;
p.knEa_max=100000;
p.kpEa_min=1000;
p.kpEa_max=100000;
p.k0_neg_min = 1.38e-7;
p.k0_neg_max = 1.38e-3;
p.k0_pos_min = 0.64e-7;
p.k0_pos_max = 0.64e-4;

p.hc_min = 1e-3;
p.hc_max = 100;
p.Cp_min = 1;
p.Cp_max = 5000;
p.V_cell_min = 1e-3;
p.V_cell_max = 1;

p.F = 96487;
% Free choice parameters
rho = 20; %rho is the capacity per m^2. Taken as a typical value based on Arora et al, 2000 (rho = 23), Ecker et al, 2015 (rho = 15), Sturm et al, 2019 (rho = 26)
p.A_surf = (Cbat/3600)/rho;
p.ce0 = 1000;
p.delta_neg =  8.5e-5;
p.delta_pos = 9.2e-5;
p.delta_sep = 2e-5;
p.L = p.delta_sep+p.delta_neg+p.delta_pos;
p.R_neg = 8.1e-6;
p.R_pos = 5.8e-6;

% Grouped parameter ranges
ph.dlnfdce_min = p.dlnfdce_min;
ph.dlnfdce_max = p.dlnfdce_max;

ph.s0_neg_min = p.s0_neg_min;
ph.s0_neg_max = p.s0_neg_max;
ph.s0_pos_min = p.s0_pos_min;
ph.s0_pos_max = p.s0_pos_max;

ph.s100_neg_min = p.s100_neg_min;
ph.s100_neg_max = p.s100_neg_max;
ph.s100_pos_min = p.s100_pos_min;
ph.s100_pos_max = p.s100_pos_max;

ph.Ds_neg_min = p.Ds_neg_min/p.R_neg^2;
ph.Ds_neg_max = p.Ds_neg_max/p.R_neg^2;
ph.Ds_pos_min = p.Ds_pos_min/p.R_pos^2;
ph.Ds_pos_max = p.Ds_pos_max/p.R_pos^2;
ph.Ds_neg_Ea_min = p.Ds_neg_Ea_min;
ph.Ds_neg_Ea_max = p.Ds_neg_Ea_max;
ph.Ds_pos_Ea_min = p.Ds_pos_Ea_min;
ph.Ds_pos_Ea_max = p.Ds_pos_Ea_max;

ph.DEa_min = p.DEa_min;
ph.DEa_max = p.DEa_max;

ph.De_min = 2.9e-7*Cbat;
ph.De_max = 7.9e-6*Cbat;

ph.p_neg_min = p.epse_neg_min/p.delta_neg;
ph.p_neg_max = p.epse_neg_max/p.delta_neg;
ph.p_sep_min = p.epse_sep_min/p.delta_sep;
ph.p_sep_max = p.epse_sep_max/p.delta_sep;
ph.p_pos_min = p.epse_pos_min/p.delta_pos;
ph.p_pos_max = p.epse_pos_max/p.delta_pos;

ph.t_plus_min = p.t_plus_min;
ph.t_plus_max = p.t_plus_max;

ph.sigma_neg_min = p.sigma_neg_min*p.epss_neg*p.A_surf/p.delta_neg;
ph.sigma_neg_max = p.sigma_neg_max*p.epss_neg*p.A_surf/p.delta_neg;
ph.sigma_pos_min = p.sigma_pos_min*p.epss_pos*p.A_surf/p.delta_pos;
ph.sigma_pos_max = p.sigma_pos_max*p.epss_pos*p.A_surf/p.delta_pos;

ph.kappa_min = 7e-6*Cbat;
ph.kappa_max = 2e-5*Cbat;

ph.Rf_neg_min = 20/Cbat;
ph.Rf_neg_max = 330/Cbat;

ph.R_cc_min = p.R_cc_min/p.A_surf;
ph.R_cc_max = p.R_cc_max/p.A_surf;

ph.alpha_min = p.alpa_min;
ph.alpha_max = p.alpa_max;

ph.knEa_min = p.knEa_min;
ph.knEa_max = p.knEa_max;
ph.kpEa_min = p.kpEa_min;
ph.kpEa_max = p.kpEa_max;

% phat2p divides by ce0^alpha, and alpha is estimated, so the lower (upper)
% bound uses the lower (upper) bound of alpha
ph.k0_neg_min = p.k0_neg_min*(p.ce0^p.alpa_min)/(p.R_neg*p.F);
ph.k0_neg_max = p.k0_neg_max*(p.ce0^p.alpa_max)/(p.R_neg*p.F);
ph.k0_pos_min = p.k0_pos_min*(p.ce0^p.alpa_min)/(p.R_pos*p.F);
ph.k0_pos_max = p.k0_pos_max*(p.ce0^p.alpa_max)/(p.R_pos*p.F);

ph.epse_neg_min = p.epse_neg_min*p.F*p.A_surf*p.delta_neg*p.ce0/(1-p.t_plus_min);
ph.epse_neg_max = p.epse_neg_max*p.F*p.A_surf*p.delta_neg*p.ce0/(1-p.t_plus_max);
ph.epse_sep_min = p.epse_sep_min*p.F*p.A_surf*p.delta_sep*p.ce0/(1-p.t_plus_min);
% The upper bound of epse_sep uses t_plus_min (instead of t_plus_max) to make
% sure that epse_sep does not exceed 1, which would happen for t_plus<0.4
% with t_plus_max
ph.epse_sep_max = p.epse_sep_max*p.F*p.A_surf*p.delta_sep*p.ce0/(1-p.t_plus_min);
ph.epse_pos_min = p.epse_pos_min*p.F*p.A_surf*p.delta_pos*p.ce0/(1-p.t_plus_min);
ph.epse_pos_max = p.epse_pos_max*p.F*p.A_surf*p.delta_pos*p.ce0/(1-p.t_plus_max);

ph.hc_min = p.hc_min;
ph.hc_max = p.hc_max;
ph.Cp_min = p.Cp_min;
ph.Cp_max = p.Cp_max;
ph.V_cell_min = p.V_cell_min;
ph.V_cell_max = p.V_cell_max;
end

function equil = determine_equilibrium_potentials(equil_data,par,options)
p = par;

if not(options.equil_mode)
    %The electrode potentials and stoichiometries are provided in equil_data
    %(they are not estimated)
    p.U_neg = equil_data.U_neg; p.U_pos = equil_data.U_pos;
    p.s0_neg = equil_data.s0_neg; p.s100_neg = equil_data.s100_neg;
    p.s0_pos = equil_data.s0_pos; p.s100_pos = equil_data.s100_pos;
end

npoints = 1000;
x = linspace(0,1,npoints);
w = linspace(p.s0_neg, p.s100_neg,npoints);                             % negative electrode stoichiometry
z = linspace(p.s0_pos, p.s100_pos,npoints);                             % positive electrode stoichiometry

if options.equil_mode
    p.U_neg = options.U_neg;
    p.EMF = equil_data.EMF;
    EMF = p.EMF(x);
    U_neg = p.U_neg(w);
    U_pos = EMF+U_neg;
    p.EMF = griddedInterpolant(x,EMF);
    x_2080 = x(npoints*0.2+1:npoints*0.8);
    w_2080 = w(npoints*0.2+1:npoints*0.8);
    
    temp = p.EMF(x_2080);
    dtemp  = (temp(2:end)-temp(1:end-1))/(x(2)-x(1));
    EMF_mean = mean(dtemp);
    
    temp = p.U_neg(w_2080);
    dtemp  = (temp(2:end)-temp(1:end-1))/(w(2)-w(1));
    U_neg_mean = mean(dtemp);
    
    U_pos_mean = -(EMF_mean+U_neg_mean);
    
    a = U_pos_mean;
    b = U_pos(npoints*0.2+1)-a*z(npoints*0.2+1);
    
    U_pos_end = @(z) a*z+b;
    
    U_pos_out1 = [U_pos_end(z(1:npoints*0.2+1)) U_pos(npoints*0.2+2:end)];
    
    U_neg_out1 = U_pos_out1-EMF;
    
    U_pos_out = griddedInterpolant(flip(z),flip(U_pos_out1));
    U_neg_out = griddedInterpolant([0 w],[2 U_neg_out1]);
    
    U_pos_full = U_pos_out(x);
    U_neg_full = U_neg_out(x);
    
    p.U_pos_out = griddedInterpolant(x,U_pos_full,'linear','linear');
    p.U_neg_out= griddedInterpolant(x,U_neg_full,'linear','linear');
    
    p.U_pos = @(z) qinterp1(x',U_pos_full',z);
    p.U_neg= @(w) qinterp1(x',U_neg_full',w);
else
    U_neg = p.U_neg(w);
    U_pos = p.U_pos(z);
    EMF = U_pos-U_neg;
    p.EMF = griddedInterpolant(x,EMF);
    p.U_neg_out = p.U_neg;
    p.U_pos_out = p.U_pos;
end

equil.U_neg = p.U_neg_out;
equil.U_pos = p.U_pos_out;
equil.Cbat = equil_data.Cbat;
if not(options.equil_mode)
    equil.s0_neg = equil_data.s0_neg; equil.s100_neg = equil_data.s100_neg;
    equil.s0_pos = equil_data.s0_pos; equil.s100_pos = equil_data.s100_pos;
end

end
