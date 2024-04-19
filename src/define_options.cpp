#include<iostream>
#include "option_defaults.h"
#include "define_options.h"
#include "define_codes.h"

using namespace inlist;

//sets up default options (all arguments defined in define_codes.h)
SimOptionList::SimOptionList():OptionList()
{
	//mutliple choice options
	this->add(TREE_KEY, std::make_shared<Option<char>>(TREE_DEFAULT, std::string(TREE_KEY), 
		                       Tree_option_list, Tree_option_name_list, TREE_OPTION_COUNT));
	this->add(ACINUS_TYPE_KEY, std::make_shared<Option<char>>(ACINUS_TYPE_DEFAULT, 
		                              std::string(ACINUS_TYPE_KEY), Acinus_type_option_list,
		                              Acinus_type_option_name_list, ACINUS_TYPE_OPTION_COUNT));
	this->add(ACINUS_MODEL_KEY, std::make_shared<Option<char>>(ACINUS_MODEL_DEFAULT, 
		                                std::string(ACINUS_MODEL_KEY), Acinus_model_option_list, 
		                                Acinus_model_option_name_list, ACINUS_MODEL_OPTION_COUNT));
	this->add(FLOW_TYPE_KEY, std::make_shared<Option<char>>(FLOW_TYPE_DEFAULT, std::string(FLOW_TYPE_KEY), 
		                 FlowType_option_list, FlowType_option_name_list, FLOW_TYPE_OPTION_COUNT));
	this->add(FLOW_INPUT_KEY, std::make_shared<Option<char>>(FLOW_INPUT_DEFAULT, std::string(FLOW_INPUT_KEY), 
		                       FlowInput_option_list, FlowInput_option_name_list, FLOW_INPUT_OPTION_COUNT));
	this->add(ELASTICITY_MODEL_KEY, std::make_shared<Option<char>>(ELASTICITY_MODEL_DEFAULT, 
		                           std::string(ELASTICITY_MODEL_KEY), ElasticityModel_option_list,
		                           ElasticityModel_option_name_list, ELASTICITY_MODEL_OPTION_COUNT));
	this->add(GRAVITY_MODEL_KEY, std::make_shared<Option<char>>(GRAVITY_MODEL_DEFAULT, 
		                         std::string(GRAVITY_MODEL_KEY), GravityModel_option_list,
		                         GravityModel_option_name_list, GRAVITY_MODEL_OPTION_COUNT));
	this->add(GRAVITY_DIRECTION_KEY, std::make_shared<Option<char>>(GRAVITY_DIRECTION_DEFAULT, 
		                         std::string(GRAVITY_DIRECTION_KEY), GravityDirection_option_list,
		                         GravityDirection_option_name_list, GRAVITY_DIRECTION_OPTION_COUNT));
	this->add(DIFFUSIVITY_TYPE_KEY, std::make_shared<Option<char>>(DIFFUSIVITY_TYPE_DEFAULT, 
		                             std::string(DIFFUSIVITY_TYPE_KEY), DiffusivityType_option_list, 
		                             DiffusivityType_option_name_list, DIFFUSIVITY_TYPE_OPTION_COUNT));
    this->add(DEP_FORMULA_SED_KEY, std::make_shared<Option<char>>(DEP_FORMULA_SED_DEFAULT,
                                                                   std::string(DEP_FORMULA_SED_KEY), SedFormula_option_list,
                                                                   SedFormula_option_name_list, DEP_SED_OPTION_COUNT));
    this->add(DEP_FORMULA_DIFF_KEY, std::make_shared<Option<char>>(DEP_FORMULA_DIFF_DEFAULT,
                                                                   std::string(DEP_FORMULA_DIFF_KEY), DiffFormula_option_list,
                                                                   DiffFormula_option_name_list, DEP_DIFF_OPTION_COUNT));
    this->add(DEP_FORMULA_IMP_KEY, std::make_shared<Option<char>>(DEP_FORMULA_IMP_DEFAULT,
                                                                   std::string(DEP_FORMULA_IMP_KEY), ImpFormula_option_list,
                                                                   ImpFormula_option_name_list, DEP_IMP_OPTION_COUNT));
	this->add(INHALATION_MODE_KEY, std::make_shared<Option<char>>(INHALATION_MODE_DEFAULT,
																std::string(INHALATION_MODE_KEY), InhalationMode_option_list,
																InhalationMode_option_name_list, INHALATION_MODE_OPTION_COUNT));
	this->add(DISCRETISATION_KEY, std::make_shared<Option<char>>(DISCRETISATION_MODE_DEFAULT,
				std::string(DISCRETISATION_KEY), DiscretisationMode_option_list,
				DiscretisationMode_option_name_list, DISCRETISATION_OPTION_COUNT));
	
	//boolean options
	this->add(ACINUS_FROM_FILE, std::make_shared<Option<bool>>(ACIN_FROM_FILE_DEFAULT, std::string(ACINUS_FROM_FILE)));
	this->add(SIMULATE_WASHIN_KEY, std::make_shared<Option<bool>>(SIMULATE_WASHIN_DEFAULT, std::string(SIMULATE_WASHIN_KEY)));	
	this->add(GRAVITY_NEGATIVE_KEY, std::make_shared<Option<bool>>(GRAVITY_NEGATIVE_DEFAULT,
		                                          std::string(GRAVITY_NEGATIVE_KEY)));
	this->add(RAD_DEPENDENT_ACIN_VOL_KEY, std::make_shared<Option<bool>>(RAD_DEPENDENT_ACIN_VOL_DEFAULT, 
		                                               std::string(RAD_DEPENDENT_ACIN_VOL_KEY)));
	this->add(DISCRETISE_BY_SCHERER_PE_KEY, std::make_shared<Option<bool>>(DISCRETISE_BY_SCHERER_PE_DEFAULT,
												std::string(DISCRETISE_BY_SCHERER_PE_KEY)));
	this->add(PRINT_ACINUS_VTK_KEY, std::make_shared<Option<bool>>(PRINT_ACINUS_VTK_DEFAULT,
		                                        std::string(PRINT_ACINUS_VTK_KEY)));
	this->add(PRINT_FLOW_VTK_KEY, std::make_shared<Option<bool>>(PRINT_FLOW_VTK_DEFAULT,
		std::string(PRINT_FLOW_VTK_KEY)));
	this->add(PRINT_TRANSPORT_VTK_KEY, std::make_shared<Option<bool>>(PRINT_TRANSPORT_VTK_DEFAULT,
		std::string(PRINT_TRANSPORT_VTK_KEY)));
	this->add(PRINT_ACINAR_AIRWAYS_KEY, std::make_shared<Option<bool>>(PRINT_ACINAR_AIRWAYS_DEFAULT,
		std::string(PRINT_ACINAR_AIRWAYS_KEY)));
	this->add(PRINT_ACINUS_CSV_KEY, std::make_shared<Option<bool>>(PRINT_ACINUS_CSV_DEFAULT, 
		                                      std::string(PRINT_ACINUS_CSV_KEY)));
	this->add(FLOW_BC_OPTION_KEY, std::make_shared<Option<bool>>(FLOW_BC_OPTION_DEFAULT, 
		                                      std::string(FLOW_BC_OPTION_KEY)));
    this->add(SIMULATE_DEPOSITION_KEY, std::make_shared<Option<bool>>(SIMULATE_DEPOSITION_DEFAULT,
            std::string(SIMULATE_DEPOSITION_KEY)));
	this->add(SIMULATE_ALVEOLAR_DEPOSITION_KEY, std::make_shared<Option<bool>>(SIMULATE_ALVEOLAR_DEPOSITION_DEFAULT,
		std::string(SIMULATE_ALVEOLAR_DEPOSITION_KEY)));
    this->add(IMPACTION_BOTH_DIRS_KEY, std::make_shared<Option<bool>>(IMPACTION_BOTH_DIRS_DEFAULT,
                                                                      std::string(IMPACTION_BOTH_DIRS_KEY)));
}

network::Position SimOptionList::get_grav_direction()
{
	double sign = 1.0;
	if(this->get_option<bool>(GRAVITY_NEGATIVE_KEY)->get_value()) sign = -1.0;
	switch(this->get_option<char>(GRAVITY_DIRECTION_KEY)->get_value())
	{
	case GRAVITY_X:
		{
			return network::Position(sign,0,0);
		} break;
	
	case GRAVITY_Y:
		{
			return network::Position(0,sign,0);
		} break;

	case GRAVITY_Z:
		{
			return network::Position(0,0,sign);
		} break;

	default:
		{
			std::cout << "Gravity direction not recognised.\n";
			return network::Position(0,0,0);
		} break;
	}
}

