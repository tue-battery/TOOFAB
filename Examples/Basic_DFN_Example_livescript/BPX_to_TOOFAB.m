%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   V 0.3
%
%
%   BPX parameter set extraction to TOOFAB parameter set
%
%   This function takes as the input the path of the BPX parameter set JSON
%   file. Only physical battery parameters are defined (simulation
%   parameters, such as spatial and time discretization, are not defined
%   here)
%
%   This file is a result of collaboration between TNO and the Control 
%   systems group, Electrical Engineering, Eindhoven University of Technology
%   and is included in the TOOlbox for Fast Battery Simulation
%   (TOOFAB)
%   Github: https://github.com/tue-battery/TOOFAB
%
%   Authors: Andres M. Diaz Aguilar, Alessio A. Lodge and Feye S.J. Hoekstra ⁎
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Transport efficiencies from BPX parameter sets are ignored as the
% porosity is already given, and they are equivalent through the relation
% porosity = tranport efficiency ^ (1/1.5)

function p = BPX_to_TOOFAB(jsonFilePath)

    %-------------------------------------------------------------------------%
    %- Non-BPX default parameters --------------------------------------------%
    %-------------------------------------------------------------------------%

    % Default parameters that are not necessarily made explicit in BPX JSON files
    p.F = 96485.33212;   % Faraday's constant [C/mol]
    p.R = 8.314462618;   % Ideal gas constant [J/mol/K]

    %-------------------------------------------------------------------------%
    %- BPX JSON parameters needed to comply with TOOFAB parameter sets -------%
    %-------------------------------------------------------------------------%

    % Load the JSON file
    jsonData = jsondecode(fileread(jsonFilePath));
    cellParams = jsonData.Parameterisation.Cell;
    electrolyteParams = jsonData.Parameterisation.Electrolyte;
    posElectrodeParams = jsonData.Parameterisation.PositiveElectrode;
    negElectrodeParams = jsonData.Parameterisation.NegativeElectrode;    
    sepParams = jsonData.Parameterisation.Separator;

    % Define mappings from JSON data to TOOFAB parameter names

    % ----------------------- 1. Cell parameters ---------------------------- %

    p.T_amb = cellParams.AmbientTemperature_K_;
    p.Vmin = cellParams.LowerVoltageCut_off_V_;
    p.Vmax = cellParams.UpperVoltageCut_off_V_;
    p.A_surf = cellParams.ElectrodeArea_m2_;
    p.T_ref = cellParams.ReferenceTemperature_K_;

    % ---------------------- 2. Electrolyte parameters ---------------------- %

    p.ce0 = electrolyteParams.InitialConcentration_mol_m_3_;    % Electrolyte initial concentration
    p.t_plus = electrolyteParams.CationTransferenceNumber;      % Li+ transference number [-]
    p.electrolyte_conductivity_EA = electrolyteParams.ConductivityActivationEnergy_J_mol_1_;
    p.electrolyte_diffusivity_EA = electrolyteParams.DiffusivityActivationEnergy_J_mol_1_;    

    % Electrolyte conductivity
    % Check if Electrolyte conductivity data is a string, constant or interpolation data
    if ischar(electrolyteParams.Conductivity_S_m_1_)
        % Case 1: Electrolyte conductivity data is given as a string expression
        kappa_ref_str = electrolyteParams.Conductivity_S_m_1_;
        kappa_ref_str = strrep(kappa_ref_str, '**', '.^');  % Replace '**' with '^'
        p.kappa_ref = str2func(['@(x)', kappa_ref_str]);    
        p.kappa = @(c_e, T) exp(p.electrolyte_conductivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.kappa_ref(c_e);
    elseif isnumeric(electrolyteParams.Conductivity_S_m_1_) && isscalar(electrolyteParams.Conductivity_S_m_1_)
        % Case 2: Electrolyte conductivity data is a single number (constant value)
        p.kappa_ref = electrolyteParams.Conductivity_S_m_1_;
        p.kappa = @(c_e, T) exp(p.electrolyte_conductivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.kappa_ref;  
    elseif isstruct(electrolyteParams.Conductivity_S_m_1_) && ...
       isfield(electrolyteParams.Conductivity_S_m_1_, 'x') && ...
       isfield(electrolyteParams.Conductivity_S_m_1_, 'y')
        % Case 3: Electrolyte conductivity data is given as interpolation data in a struct with fields 'x' and 'y'
        x_data_kappa = electrolyteParams.Conductivity_S_m_1_.x;
        y_data_kappa = electrolyteParams.Conductivity_S_m_1_.y;
        p.kappa_ref = @(input_x) interp1(x_data_kappa, y_data_kappa, input_x, 'linear', 'extrap');
        p.kappa = @(c_e, T) exp(p.electrolyte_conductivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.kappa_ref(c_e);
    else 
        disp('Error: The electrolyte conductivity data is neither a constant, a function, or a dataset for interpolation. As such, it cannot be calculated and is not part of the parameter set.');
    end

    % Electrolyte diffusivity
    % Check if electrolyte diffusivity data is a string, constant or interpolation data
    if ischar(electrolyteParams.Diffusivity_m2_s_1_)
        % Case 1: Electrolyte diffusivity data is given as a string expression
        De_ref_str = electrolyteParams.Diffusivity_m2_s_1_;
        De_ref_str = strrep(De_ref_str, '**', '.^');  % Replace '**' with '^'
        p.De_ref = str2func(['@(x)', De_ref_str]);    
        p.De = @(c_e, T) exp(p.electrolyte_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.De_ref(c_e);
    elseif isnumeric(electrolyteParams.Diffusivity_m2_s_1_) && isscalar(electrolyteParams.Diffusivity_m2_s_1_)
        % Case 2: Electrolyte diffusivity data is a single number (constant value)
        p.De_ref = electrolyteParams.Diffusivity_m2_s_1_;
        p.De = @(c_e, T) exp(p.electrolyte_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.De_ref;  
    elseif isstruct(electrolyteParams.Diffusivity_m2_s_1_) && ...
       isfield(electrolyteParams.Diffusivity_m2_s_1_, 'x') && ...
       isfield(electrolyteParams.Diffusivity_m2_s_1_, 'y')
        % Case 3: Electrolyte diffusivity data is given as interpolation data in a struct with fields 'x' and 'y'
        x_data_De = electrolyteParams.Diffusivity_m2_s_1_.x;
        y_data_De = electrolyteParams.Diffusivity_m2_s_1_.y;
        p.De_ref = @(input_x) interp1(x_data_De, y_data_De, input_x, 'linear', 'extrap');
        p.De = @(c_e, T) exp(p.electrolyte_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.De_ref(c_e);
    else 
        disp('Error: The electrolyte diffusivity data is neither a constant, a function, or a dataset for interpolation. As such, it cannot be calculated and is not part of the parameter set.');
    end

    % ------------------ 3. Negative electrode parameters ------------------- %

    p.delta_neg = negElectrodeParams.Thickness_m_;
    p.R_neg = negElectrodeParams.ParticleRadius_m_;

    % Check if OCP data is a string expression or interpolation data
    if ischar(negElectrodeParams.OCP_V_)
        % Case 1: OCP is given as a string expression
        U_neg_str = negElectrodeParams.OCP_V_;
        U_neg_str = strrep(U_neg_str, '**', '.^');  % Replace '**' with '^'
        p.U_neg = str2func(['@(x)', U_neg_str]);
    elseif isstruct(negElectrodeParams.OCP_V_) && ...
       isfield(negElectrodeParams.OCP_V_, 'x') && ...
       isfield(negElectrodeParams.OCP_V_, 'y')
        % Case 2: OCP is given as interpolation data in a struct with fields 'x' and 'y'
        x_data_OCP_neg = negElectrodeParams.OCP_V_.x;
        y_data_OCP_neg = negElectrodeParams.OCP_V_.y;
        p.U_neg = @(x) interp1(x_data_OCP_neg, y_data_OCP_neg, x, 'linear', 'extrap');
    else 
        disp('Error: The Negative OCP data is neither a function nor a dataset for interpolation. As such, it cannot be imported into the parameter set.');
    end

    p.sigma_neg = negElectrodeParams.Conductivity_S_m_1_;
    p.s100_neg = negElectrodeParams.MaximumStoichiometry; % Stoichiometry at 100% SoC at the neg. electrode [-]
    p.s0_neg = negElectrodeParams.MinimumStoichiometry;
    p.cs_max_neg = negElectrodeParams.MaximumConcentration_mol_m_3_;

    % Volume Fractions
    p.epse_neg = negElectrodeParams.Porosity;    %Electrolyte volume fraction at the neg. electrode [-] 

    p.neg_surface_area_per_unit_volume = negElectrodeParams.SurfaceAreaPerUnitVolume_m_1_; % unused as TOOFAB uses 3*p.epss_neg/p.R_neg (3*active_volume_frac/particle_radius) instead
    p.epss_neg = (p.neg_surface_area_per_unit_volume * p.R_neg)/3; %Active material volume fraction at the neg. electrode [-]
    
    %%%%%

    % Negative electrode diffusivity
    p.neg_electrode_diffusivity_EA = negElectrodeParams.DiffusivityActivationEnergy_J_mol_1_;    
    % Check if electrolyte diffusivity data is a string, constant or interpolation data
    if ischar(negElectrodeParams.Diffusivity_m2_s_1_)
        % Case 1: Negative electrode diffusivity data is given as a string expression
        Ds_neg_ref_str = negElectrodeParams.Diffusivity_m2_s_1_;
        Ds_neg_ref_str = strrep(Ds_neg_ref_str, '**', '^');  % Replace '**' with '^'
        p.Ds_neg_ref = str2func(['@(x)', Ds_neg_ref_str]);    
        p.Ds_neg = @(stoich, T) exp(p.neg_electrode_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.Ds_neg_ref(stoich); 
    elseif isnumeric(negElectrodeParams.Diffusivity_m2_s_1_) && isscalar(negElectrodeParams.Diffusivity_m2_s_1_)
        % Case 2: Negative electrode diffusivity data is a single number (constant value)
        p.Ds_neg_ref = negElectrodeParams.Diffusivity_m2_s_1_;

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        p.Ds_neg = @(stoich, T) exp(p.neg_electrode_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.Ds_neg_ref; % Original version with stoich dependency  
        % p.Ds_neg = @(T) exp(p.neg_electrode_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.Ds_neg_ref;  
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

    elseif isstruct(negElectrodeParams.Diffusivity_m2_s_1_) && ...
       isfield(negElectrodeParams.Diffusivity_m2_s_1_, 'x') && ...
       isfield(negElectrodeParams.Diffusivity_m2_s_1_, 'y')
        % Case 3: Negative electrode diffusivity data is given as interpolation data in a struct with fields 'x' and 'y'
        x_data_Ds_neg = negElectrodeParams.Diffusivity_m2_s_1_.x;
        y_data_Ds_neg = negElectrodeParams.Diffusivity_m2_s_1_.y;
        p.Ds_neg_ref = @(input_x) interp1(x_data_Ds_neg, y_data_Ds_neg, input_x, 'linear', 'extrap');
        p.Ds_neg = @(stoich, T) exp(p.neg_electrode_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.Ds_neg_ref(stoich);
    else 
        disp('Error: The negative electrode diffusivity data is neither a constant, a function, or a dataset for interpolation. As such, it cannot be calculated and is not part of the parameter set.');
    end

    % Reaction rate constant --> different units on TOOFAB than on BPX %
    p.k0_neg_EA = negElectrodeParams.ReactionRateConstantActivationEnergy_J_mol_1_; % Reaction Rate Constant Activation Energy
    p.k0_neg_ref = negElectrodeParams.ReactionRateConstant_mol_m_2_s_1_ * p.F / (p.cs_max_neg * p.ce0^0.5); % Reference reaction rate constant to covert to TOOFAB's form
    p.k0_neg = @(T) exp(p.k0_neg_EA / p.R * (1 / p.T_ref - 1 / T)) * p.k0_neg_ref; % Kinetic constant in the neg. electrode [mol^(-3/2)*m^(-1/2)*s^(-1)]
    
    % Entropic Change Coefficient (function of stoich)
    % Check if Entropic change data is a string expression or interpolation data
    if ischar(negElectrodeParams.EntropicChangeCoefficient_V_K_1_)
        % Case 1: Entropic change data is givepn as a string expression
        dUdt_neg_str = negElectrodeParams. EntropicChangeCoefficient_V_K_1_;
        dUdt_neg_str = strrep(dUdt_neg_str, '**', '.^');  % Replace '**' with '^'
        p.dU_dT_neg = str2func(['@(x)', dUdt_neg_str]);
    elseif isnumeric(negElectrodeParams.EntropicChangeCoefficient_V_K_1_) && isscalar(negElectrodeParams.EntropicChangeCoefficient_V_K_1_)
        % Case 2: Entropic change data is a single number (constant value)
        constant_value = negElectrodeParams.EntropicChangeCoefficient_V_K_1_;
        p.dU_dT_neg = @(x) constant_value * ones(size(x));  
    elseif isstruct(negElectrodeParams.EntropicChangeCoefficient_V_K_1_) && ...
       isfield(negElectrodeParams.EntropicChangeCoefficient_V_K_1_, 'x') && ...
       isfield(negElectrodeParams.EntropicChangeCoefficient_V_K_1_, 'y')
        % Case 3: Entropic change data is given as interpolation data in a struct with fields 'x' and 'y'
        x_data_neg = negElectrodeParams.EntropicChangeCoefficient_V_K_1_.x;
        y_data_neg = negElectrodeParams.EntropicChangeCoefficient_V_K_1_.y;
        p.dU_dT_neg = @(x) interp1(x_data_neg, y_data_neg, x, 'linear', 'extrap');
    else 
        disp('Error: The Negative Electrode Entropic Change Coefficient Data is neither a constant, a function, or a dataset for interpolation. As such, it cannot be calculated and is not part of the parameter set.');
    end

    % ------------------ 4. Positive electrode parameters ------------------- %

    p.delta_pos = posElectrodeParams.Thickness_m_;
    p.R_pos = posElectrodeParams.ParticleRadius_m_;

    % Check if OCP data is a string expression or interpolation data
    if ischar(posElectrodeParams.OCP_V_)
        % Case 1: OCP is given as a string expression
        U_pos_str = posElectrodeParams.OCP_V_;
        U_pos_str = strrep(U_pos_str, '**', '.^');  % Replace '**' with '^'
        p.U_pos = str2func(['@(x)', U_pos_str]);    
    elseif isstruct(posElectrodeParams.OCP_V_) && ...
       isfield(posElectrodeParams.OCP_V_, 'x') && ...
       isfield(posElectrodeParams.OCP_V_, 'y')
        % Case 2: OCP is given as interpolation data in a struct with fields 'x' and 'y'
        x_data_OCP_pos = posElectrodeParams.OCP_V_.x;
        y_data_OCP_pos = posElectrodeParams.OCP_V_.y;
        p.U_pos = @(x) interp1(x_data_OCP_pos, y_data_OCP_pos, x, 'linear', 'extrap');
    else 
        disp('Error: The Negative OCP data is neither a function nor a dataset for interpolation. As such, it cannot be imported into the parameter set.');
    end

    p.sigma_pos = posElectrodeParams.Conductivity_S_m_1_;
    p.s100_pos = posElectrodeParams.MinimumStoichiometry;           % Stoichiometry at 100% SoC at the pos. electrode [-]
    p.s0_pos = posElectrodeParams.MaximumStoichiometry;
    p.cs_max_pos = posElectrodeParams.MaximumConcentration_mol_m_3_;

    % Volume fractions
    p.epse_pos = posElectrodeParams.Porosity;       % Electrolyte volume fraction at the neg. electrode [-]

    p.pos_surface_area_per_unit_volume = posElectrodeParams.SurfaceAreaPerUnitVolume_m_1_; 
    % unused as TOOFAB uses 3*p.epss_pos/p.R_pos (3*active_volume_frac/particle_radius) to calculate the surface area
    % per unit volume in their DFN code

    p.epss_pos = (p.pos_surface_area_per_unit_volume * p.R_pos)/3; %Active material volume fraction at the pos. electrode [-] 
    
    %%%%%

    % Positive electrode diffusivity
    p.pos_electrode_diffusivity_EA = posElectrodeParams.DiffusivityActivationEnergy_J_mol_1_;    
    % Check if positive electrode diffusivity data is a string, constant or interpolation data
    if ischar(posElectrodeParams.Diffusivity_m2_s_1_)
        % Case 1: Positive electrode diffusivity data is given as a string expression
        Ds_pos_ref_str = posElectrodeParams.Diffusivity_m2_s_1_;
        Ds_pos_ref_str = strrep(Ds_pos_ref_str, '**', '.^');  % Replace '**' with '^'
        p.Ds_pos_ref = str2func(['@(x)', Ds_pos_ref_str]);    
        p.Ds_pos = @(stoich, T) exp(p.pos_electrode_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.Ds_pos_ref(stoich);
    elseif isnumeric(posElectrodeParams.Diffusivity_m2_s_1_) && isscalar(posElectrodeParams.Diffusivity_m2_s_1_)
        % Case 2: Positive electrode diffusivity data is a single number (constant value)
        p.Ds_pos_ref = posElectrodeParams.Diffusivity_m2_s_1_;
        p.Ds_pos = @(stoich, T) exp(p.pos_electrode_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.Ds_pos_ref;  
    elseif isstruct(posElectrodeParams.Diffusivity_m2_s_1_) && ...
       isfield(posElectrodeParams.Diffusivity_m2_s_1_, 'x') && ...
       isfield(posElectrodeParams.Diffusivity_m2_s_1_, 'y')
        % Case 3: Positive electrode diffusivity data is given as interpolation data in a struct with fields 'x' and 'y'
        x_data_Ds_pos = posElectrodeParams.Diffusivity_m2_s_1_.x;
        y_data_Ds_pos = posElectrodeParams.Diffusivity_m2_s_1_.y;
        p.Ds_pos_ref = @(input_x) interp1(x_data_Ds_pos, y_data_Ds_pos, input_x, 'linear', 'extrap');
        p.Ds_pos = @(stoich, T) exp(p.pos_electrode_diffusivity_EA / p.R * (1 / p.T_ref - 1 / T)) * p.Ds_pos_ref(stoich);
    else 
        disp('Error: The positive electrode diffusivity data is neither a constant, a function, or a dataset for interpolation. As such, it cannot be calculated and is not part of the parameter set.');
    end

    % Reaction rate constant --> different units on TOOFAB than on BPX %
    p.k0_pos_EA = posElectrodeParams.ReactionRateConstantActivationEnergy_J_mol_1_; %Reaction Rate Constant Activation Energy
    p.k0_pos_ref = posElectrodeParams.ReactionRateConstant_mol_m_2_s_1_ * p.F / (p.cs_max_pos * p.ce0^0.5);
    p.k0_pos = @(T) exp(p.k0_pos_EA / p.R * (1 / p.T_ref - 1 / T)) * p.k0_pos_ref; %Kinetic constant in the neg. electrode [mol^(-3/2)*m^(-1/2)*s^(-1)]

    % Entropic Change Coefficient (function of stoich)
    % Check if Entropic change data is a string expression or interpolation data
    if ischar(posElectrodeParams.EntropicChangeCoefficient_V_K_1_)
        % Case 1: Entropic change data is given as a string expression
        dUdt_pos_str = posElectrodeParams. EntropicChangeCoefficient_V_K_1_;
        dUdt_pos_str = strrep(dUdt_pos_str, '**', '.^');  % Replace '**' with '^'
        p.dU_dT_pos = str2func(['@(x)', dUdt_pos_str]);
    elseif isnumeric(posElectrodeParams.EntropicChangeCoefficient_V_K_1_) && isscalar(posElectrodeParams.EntropicChangeCoefficient_V_K_1_)
        % Case 2: Entropic change data is a single number (constant value)
        constant_value = posElectrodeParams.EntropicChangeCoefficient_V_K_1_;
        p.dU_dT_pos = @(x) constant_value * ones(size(x));  
    elseif isstruct(posElectrodeParams.EntropicChangeCoefficient_V_K_1_) && ...
       isfield(posElectrodeParams.EntropicChangeCoefficient_V_K_1_, 'x') && ...
       isfield(posElectrodeParams.EntropicChangeCoefficient_V_K_1_, 'y')
        % Case 3: Entropic change data is given as interpolation data in a struct with fields 'x' and 'y'
        x_data_pos = posElectrodeParams.EntropicChangeCoefficient_V_K_1_.x;
        y_data_pos = posElectrodeParams.EntropicChangeCoefficient_V_K_1_.y;
        p.dU_dT_pos = @(x) interp1(x_data_pos, y_data_pos, x, 'linear', 'extrap');
    else 
        disp('The Positive Electrode Entropic Change Coefficient Data is neither a constant, a function, or a dataset for interpolation. As such, it cannot be calculated and is not part of the parameter set.');
    end

    % ---------------------- 5. Separator parameters ------------------------ %

    p.delta_sep = sepParams.Thickness_m_;
    p.epse_sep = sepParams.Porosity;  

    % ----------------------- 6. Thermal parameters ------------------------- %

    % p.hc = 0.8915; %Convective heat transfer coefficient [W/K]

    p.volume = cellParams.Volume_m3_;   
    p.A_external = cellParams.ExternalSurfaceArea_m2_; % Cell cooling surface area
    p.hc = cellParams.ThermalConductivity_W_m_1_K_1_*p.A_external; %Convective heat transfer coefficient [W/K]
    p.Cp = cellParams.SpecificHeatCapacity_J_K_1_kg_1_; % Heat capacity [J/kg/K]

    %-------------------------------------------------------------------------%
    %- Extra BPX parameters that are included in the JSON file, not needed to %
    %- comply with TOOFAB parameter sets (here for completeness, and thus) ---%
    %-------------- will not be used by TOOFAB's solver) ---------------------%

    % ----------------------- 1. Cell parameters ---------------------------- %

    p.T_init = cellParams.InitialTemperature_K_;
    p.density = cellParams.Density_kg_m_3_;
    p.num_electrode_pairs_parallel = cellParams.NumberOfElectrodePairsConnectedInParallelToMakeACell;
    p.m = p.density * p.volume;     % Is used by TOOFAB!

    % ------------------- 3/4. Electrode Parameters ------------------------ %
    
    p.t_eff_pos = posElectrodeParams.TransportEfficiency; % Transport efficiency in positive electrode (unused by TOOFAB)
    p.t_eff_neg = negElectrodeParams.TransportEfficiency; % Transport efficiency in negative electrode (unused by TOOFAB)
    p.t_eff_sep = sepParams.TransportEfficiency; % Transport efficiency in positive electrode (unused by TOOFAB)

% with density and volume you can find the required mass

    % ---- ASSUMMED PARAMETERS THAT ARE NOT ON THE JSON FILE NEEDED FOR TOOFAB ---- %
   
    p.R_cc = 0;  % Overall current collector resistance
    
    % Transfer coefficient of surface reaction
    p.alpha_a = 0.5;                                                        % Charge transfer coefficent (anodic) [-]
    p.alpha_c = (1-p.alpha_a);                                              % Charge transfer coefficient (cathodic) [-]

    % Bruggeman constant (numbers from PyBaMM translation of BPX parameter set)
    p.p_neg = 1.5;  %Bruggeman constant at the neg. electrode [-]
    p.p_pos = 1.5;  %Bruggeman constant at the pos. electrode [-]
    p.p_sep = 1.5;  %Bruggeman constant at the seperator [-] 

    p.Rf0_neg = 0; 
    p.Rf0_pos = 0;  
    p.Rf_neg = p.Rf0_neg; %SEI film resistance in the positive electrode
    p.Rf_pos = p.Rf0_pos; %SEI film resistance in the negative electrode

    p.Cbat = -(p.s100_pos-p.s0_pos)*p.epss_pos*p.delta_pos*p.A_surf*p.F*p.cs_max_pos; % A.s
    p.nominal_Cbat = cellParams.NominalCellCapacity_A_h_; % Ah

    % Reversible capacity of the battery [C]

    % ---- ASSUMMED PARAMETERS THAT ARE NOT ON THE JSON FILE (AS DONE ON PYBAMM) ---- %
    
    p.dlnfdce = 0; % Often assumed to be 0 in literature
    
    % Missing: electrode width [m], electrode height [m]
    
    % Negative electrode density, positive electrode density, separator density, 
    % negative current collector density, positive current collector density
    % (all assumed to be the same as the regular p.density on PyBaMM)

    % Same story with thermal conductivity (same for all sections of cell)
    % Same story with specific heat capacity

    % CC thickness = 0
    % CC conductivity = 40000000.0

    % Total heat transfer coefficient [W.m-2.K-1]: 0

    % Initial concentration in negative electrode [mol.m-3]: 22496.096400000002
    % Initial concentration in positive electrode [mol.m-3]: 19599.888
    % --> check how these are calculated!

    % ---------- Display loaded parameters -----------%

    disp('BPX parameters loaded');
    % disp(p);
    
end