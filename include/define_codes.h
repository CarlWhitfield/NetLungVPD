#ifndef DEFINE_CODES_H
#define DEFINE_CODES_H
#include <vector>
#include <string>

//option codes for tree structure
#define TREE_KEY "Tree"
#define TREE_OPTION_COUNT 2
#define MODEL_M 'm'
#define MODEL_FROM_FILE 'f'
const char Tree_option_list[] = {MODEL_M, MODEL_FROM_FILE};
const std::string Tree_option_name_list[] = {"Model_M", "Model_from_file"};

//option codes for acinus
#define ACINUS_TYPE_KEY "Acinus_type"
#define ACINUS_TYPE_OPTION_COUNT 3
#define WELL_MIXED_BAG 'b'
#define TREE_ACINUS 't'
#define DETAILED_ACINUS 'd'
const char Acinus_type_option_list[] = {WELL_MIXED_BAG, TREE_ACINUS, DETAILED_ACINUS};
const std::string Acinus_type_option_name_list[] = {"Well_mixed_bag", "Tree_acinus", "Detailed_acinus"};

//option codes for acinus
#define ACINUS_MODEL_KEY "Acinus_model"
#define ACINUS_MODEL_OPTION_COUNT 3
#define WEIBEL_ACIN 'w'
#define DUTRIEUE_ACIN 'd'
#define HENRY_ACIN 'h'
const char Acinus_model_option_list[] = {WEIBEL_ACIN, DUTRIEUE_ACIN, HENRY_ACIN};
const std::string Acinus_model_option_name_list[] = {"Weibel_symm_acinus", "Dutrieue_acinus", "Henry_acinus"};

//option codes for flow type
#define FLOW_TYPE_KEY "Flow_Type"
#define FLOW_TYPE_OPTION_COUNT 2
#define POISEUILLE 'p'
#define PEDLEY 'n'
const char FlowType_option_list[] = {POISEUILLE, PEDLEY};
const std::string FlowType_option_name_list[] = {"Poiseuille_flow", "Pedley_flow"};

//option codes for flow constraint
#define FLOW_INPUT_KEY "Flow_Input"
#define FLOW_INPUT_OPTION_COUNT 3
#define SIN_FLOW 's'
#define STEP_FLOW 'd'
#define FLOW_FROM_FILE 'f'
const char FlowInput_option_list[] = {SIN_FLOW, STEP_FLOW, FLOW_FROM_FILE};
const std::string FlowInput_option_name_list[] = {"Sinusoidal_breathing", "Step_function_flow", "Breathing_pattern_from_file"};

//option codes for elasticity model
#define ELASTICITY_MODEL_KEY "Elasticity_Model"
#define ELASTICITY_MODEL_OPTION_COUNT 2
#define LINEAR_ELASTICITY 'l'
#define TAWHAI_ELASTICITY 't'
const char ElasticityModel_option_list[] = {LINEAR_ELASTICITY, TAWHAI_ELASTICITY};
const std::string  ElasticityModel_option_name_list[] = {"Linear_elasticity", "Nonlinear_elasticity_tawhai"};

//option codes for gravity effect
#define GRAVITY_MODEL_KEY "Gravity_Model"
#define GRAVITY_MODEL_OPTION_COUNT 2
#define NO_GRAVITY 'n'
#define LINEAR_GRAVITY_APPROX 'l'
const char GravityModel_option_list[] = {NO_GRAVITY, LINEAR_GRAVITY_APPROX};
const std::string  GravityModel_option_name_list[] = {"No_gravity", "Linear_gravity_approximation"};

//option codes for gravity direction
#define GRAVITY_DIRECTION_KEY "Gravity_Direction"
#define GRAVITY_DIRECTION_OPTION_COUNT 3
#define GRAVITY_X 'x'
#define GRAVITY_Y 'y'
#define GRAVITY_Z 'z'
const char GravityDirection_option_list[] = {GRAVITY_X, GRAVITY_Y, GRAVITY_Z};
const std::string  GravityDirection_option_name_list[] = {"X", "Y", "Z"};

//option codes for diffusivity type
#define DIFFUSIVITY_TYPE_KEY "Diffusivity_Type"
#define DIFFUSIVITY_TYPE_OPTION_COUNT 4
#define MOLECULAR_DIFFUSIVITY 'm'
#define TAYLOR_DISPERSION 't'
#define SCHERER_DISPERSION 's'
#define TAULBEE_DISPERSION 'b'
const char DiffusivityType_option_list[] = {MOLECULAR_DIFFUSIVITY, TAYLOR_DISPERSION, SCHERER_DISPERSION, TAULBEE_DISPERSION};
const std::string DiffusivityType_option_name_list[] = {"Molecular_diffusivity", "Taylor_dispersion", "Scherer_dispersion","Taulbee_dispersion"};

//option codes for diffusion deposition term
#define DEP_FORMULA_DIFF_KEY "Deposition_formula_diffusion"
#define DEP_DIFF_OPTION_COUNT 4
#define SIMPLE_DIFF_DEP 's'
#define INGHAM_DIFF_DEP 'i'
#define YU_DIFF_DEP 'y'
#define INGHAM91_DIFF_DEP 'j'
const char DiffFormula_option_list[] = {SIMPLE_DIFF_DEP, INGHAM_DIFF_DEP, YU_DIFF_DEP, INGHAM91_DIFF_DEP};
const std::string DiffFormula_option_name_list[] = {"Simple_diffusion_deposition",
                                                    "Ingham_diffusion_deposition", "Yu_diffusion_deposition",
                                                    "Ingham91_diffusion_deposition"};

//option codes for sedimentation deposition term
#define DEP_FORMULA_SED_KEY "Deposition_formula_sedimentation"
#define DEP_SED_OPTION_COUNT 3
#define SIMPLE_SED_DEP 's'
#define PICH_SED_DEP 'p'
#define WANG_SED_DEP 'w'
const char SedFormula_option_list[] = {SIMPLE_SED_DEP, PICH_SED_DEP, WANG_SED_DEP};
const std::string SedFormula_option_name_list[] = {"Simple_sedimentation_deposition",
                                                    "Pich_sedimentation_deposition",
                                                    "Wang_sedimentation_deposition"};

//option codes for diffusion deposition term
#define DEP_FORMULA_IMP_KEY "Deposition_formula_impaction"
#define DEP_IMP_OPTION_COUNT 2
#define YEH_IMP_DEP 'y'
#define ZHANG_IMP_DEP 'z'
const char ImpFormula_option_list[] = {YEH_IMP_DEP, ZHANG_IMP_DEP};
const std::string ImpFormula_option_name_list[] = {"Yeh_impaction_deposition",
                                                    "Zhang_impaction_deposition"};

//option codes for inspiration mode
#define INHALATION_MODE_KEY "Inhalation_mode"
#define INHALATION_MODE_OPTION_COUNT 2
#define CONTINUOUS_CONC 'c'
#define BOLUS_CONC 'b'
const char InhalationMode_option_list[] = { CONTINUOUS_CONC, BOLUS_CONC };
const std::string InhalationMode_option_name_list[] = { "Continuous conc inhalation",
														"Bolus conc inhalation" };

#define DISCRETISATION_KEY "Discretisation_mode"
#define DISCRETISATION_OPTION_COUNT 2
#define MAX_MIN_DISC 'm'
#define MAX_FV_PECLET_DISC 'p'
const char DiscretisationMode_option_list[] = { MAX_MIN_DISC, MAX_FV_PECLET_DISC };
const std::string DiscretisationMode_option_name_list[] = { "Max size min no. discretisation",
															"Max Peclet discretisation" };

//conversions
#define L_TO_M3 0.001
#define MM_TO_M 0.001
#define CMH20_TO_PA 98.0665

//conversion keys
#define VOL_KEY "Volume"
#define TIME_KEY "Time"
#define PRESSURE_KEY "Pressure"
#define LENGTH_KEY "Length"
#define RECIPROCAL_LENGTH_KEY "Reciprocal_length"
#define ELASTANCE_KEY "Elastance"
#define RESISTANCE_KEY "Resistance"
#define INERTANCE_KEY "Inertance"
#define VISCOSITY_KEY "Viscosity"
#define DENSITY_KEY "Density"
#define DIFFUSIVITY_KEY "Diffusivity"
#define FORCE_DENSITY_KEY "ForceDensity"
#define ACCELERATION_KEY "Acceleration"
#define DENSITY_KEY "Density"
#define VOL_FLUX_KEY "Volume_flux"

#define OPTIONS_FILE_EXT "options"
#define PARAMS_FILE_EXT "params"
#define NODE_FILE_EXT "nodes"
#define BRANCH_FILE_EXT "branches"
#define TERM_NODE_FILE_EXT "termnodes"
#define FLOW_INPUT_FILE_EXT "flow"
#define TNODE_MAP_FILE_EXT "map"

//parameter keys
#define MIN_FV_PER_GEN_KEY "Min_elements_per_airway"
#define AIRWAY_DEAD_SPACE_KEY "Airway_dead_space"
#define FRC_VOLUME_KEY "FRC"
#define TIDAL_VOLUME_KEY "Tidal_volume"
#define PLEURAL_PRESSURE_CHANGE_KEY "Pleural_DP"
#define MOUTH_VOLUME_KEY "Mouth_volume"
#define BOLUS_START_INSP_KEY "Bolus_start_volume"
#define BOLUS_VOL_DURATION_KEY "Bolus_volume"
#define TIDAL_TIME_KEY "Tidal_breath_time"
#define LUNG_ELASTANCE_KEY "Lung_elastance"
#define BAG_RESISTANCE_KEY "Bag_resistance"
#define BAG_INERTANCE_KEY "Bag_inertance"
#define MOUTH_RESISTANCE_KEY "Mouth_resistance"
#define AIR_VISCOSITY_KEY "Air_viscosity"
#define AIR_DENSITY_KEY "Air_density"
#define TIME_STEP_KEY "Time_step"
#define MAX_FV_LENGTH_KEY "Max_element_length"
#define MAX_ELEMENT_PE_KEY "Max_element_peclet"
#define GAS_DIFFUSIVITY_KEY "Gas_diffusivity"
#define MUCUS_AIR_PARTITION_COEFFICIENT_KEY "Mucus_air_partition_coefficient"
#define PARTICLE_SIZE_KEY "Particle_size"
#define PARTICLE_DENSITY_KEY "Particle_density"
#define GRAV_ACCEL_KEY "Gravity_acceleration"
#define BOLTZMANN_TEMP_VISCOSITY_KEY "Boltzmann_temp_viscosity"
#define MEAN_FREE_PATH_KEY "Mean_free_path"

#define WASHIN_END_KEY "Washin_end_time"
#define ACIN_ASYMM_FACTOR_KEY "Acin_asymm_factor"
#define ACIN_SCALE_FACTOR_KEY "Acin_scale_factor"
#define ACIN_FV_PER_GEN_KEY "Acin_elements_per_airway"
#define ACIN_LR_RATIO_KEY "Acin_length_rad_ratio"
#define GRAVITY_PRESSURE_GRADIENT_KEY "Gravity_press_gradient"

#define PRINT_TIME_KEY "Print_timestep"
#define BREATH_COUNT_KEY "Total_breaths"

//bool option keys
#define ACINUS_FROM_FILE "Acinus_files"
#define SIMULATE_WASHIN_KEY "Simulate_washin"
#define GRAVITY_NEGATIVE_KEY "Gravity_negative"
#define GAS_EXCHANGE_KEY "Gas_exchange"
#define RAD_DEPENDENT_ACIN_VOL_KEY "Acinus_volumes_vary_with_term_radius"
#define PRINT_ACINUS_VTK_KEY "Print_acini_vtk"
#define PRINT_ACINUS_CSV_KEY "Print_acini_csv"
#define FLOW_BC_OPTION_KEY "Flow_BC"
#define SIMULATE_DEPOSITION_KEY "Sim_deposition"
#define SIMULATE_ALVEOLAR_DEPOSITION_KEY "Alv_deposition"
#define DISCRETISE_BY_SCHERER_PE_KEY "Peclet_discretisation"
#define PRINT_FLOW_VTK_KEY "Print_flow_vtk"
#define PRINT_TRANSPORT_VTK_KEY "Print_transport_vtk"
#define PRINT_ACINAR_AIRWAYS_KEY "Print_acinar_airways"
#define IMPACTION_BOTH_DIRS_KEY "Impaction_both_directions"


#endif
