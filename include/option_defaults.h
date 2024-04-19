#ifndef OPTION_DEFAULTS_H
#define OPTION_DEFAULTS_H

#include "define_codes.h"

#define TREE_DEFAULT MODEL_FROM_FILE
#define ACINUS_TYPE_DEFAULT TREE_ACINUS
#define ACINUS_MODEL_DEFAULT WEIBEL_ACIN
#define FLOW_TYPE_DEFAULT POISEUILLE
#define FLOW_INPUT_DEFAULT SIN_FLOW
#define ELASTICITY_MODEL_DEFAULT LINEAR_ELASTICITY
#define GRAVITY_MODEL_DEFAULT NO_GRAVITY
#define GRAVITY_DIRECTION_DEFAULT GRAVITY_Z
#define DIFFUSIVITY_TYPE_DEFAULT TAYLOR_DISPERSION
#define DEP_FORMULA_DIFF_DEFAULT YU_DIFF_DEP
#define DEP_FORMULA_SED_DEFAULT SIMPLE_SED_DEP
#define DEP_FORMULA_IMP_DEFAULT ZHANG_IMP_DEP
#define INHALATION_MODE_DEFAULT	CONTINUOUS_CONC
#define DISCRETISATION_MODE_DEFAULT MAX_MIN_DISC
#define SIMULATE_DEPOSITION_DEFAULT true
#define SIMULATE_ALVEOLAR_DEPOSITION_DEFAULT false
#define ACIN_FROM_FILE_DEFAULT false
#define GRAVITY_NEGATIVE_DEFAULT true
#define SIMULATE_WASHIN_DEFAULT false
#define RAD_DEPENDENT_ACIN_VOL_DEFAULT false
#define DISCRETISE_BY_SCHERER_PE_DEFAULT false
#define PRINT_ACINUS_VTK_DEFAULT false
#define PRINT_FLOW_VTK_DEFAULT true
#define PRINT_TRANSPORT_VTK_DEFAULT true
#define PRINT_ACINAR_AIRWAYS_DEFAULT false
#define PRINT_ACINUS_CSV_DEFAULT false
#define FLOW_BC_OPTION_DEFAULT true
#define IMPACTION_BOTH_DIRS_DEFAULT true
	//RespFunc;       
	//BreathOp;     
	//ShapeOp;        
	//SolverOp;
	//InitOp;        
	//InputOp;         
	//TaylorDisp;    
	//OutputOp;   
	//PrintOp;    
	//output_perts;  
	//output_lung_volumes;
	//include_acin_perts;

#endif