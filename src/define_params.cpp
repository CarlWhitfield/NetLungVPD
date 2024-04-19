#include "define_codes.h"
#include "define_params.h"
#include "param_defaults.h"

void PhysParam::update_value( const double & val ) //this is to be used to update the sim value parameter directly
{
	Parameter<double>::update_value(val);
	if(this->phys_to_sim_conversion == NULL)
	{
		std::cout << "Warning, failed to update Phys and SI values when updating " << name << "\n";
	}
	else
	{
		PhysValue = val / (*(this->phys_to_sim_conversion));
		this->calc_SI_from_Phys_value();
	}
}

std::string PhysParam::phys_value_string()
{
	std::stringstream ss;
	ss << PhysValue << PhysUnits;
	return ss.str();
}

SimParameterList::SimParameterList():ParameterList()
{
	using namespace inlist;
	//integer parameters
	this->add(MIN_FV_PER_GEN_KEY, std::make_shared<Parameter<int>>(MIN_FV_PER_GEN_DEFAULT, 
		                                                      std::string(MIN_FV_PER_GEN_KEY)));
	this->add(ACIN_FV_PER_GEN_KEY, std::make_shared<Parameter<int>>(ACIN_FV_PER_GEN_DEFAULT, 
		                                                      std::string(ACIN_FV_PER_GEN_KEY)));
	//double parameters need to be converted to common units
	//all set to 1 initially
	this->set_conversion(VOL_KEY, 1.0);
	this->set_conversion(LENGTH_KEY, 1.0);
	this->set_conversion(RECIPROCAL_LENGTH_KEY, 1.0);
	this->set_conversion(TIME_KEY, 1.0);
	this->set_conversion(PRESSURE_KEY, 1.0);
	this->set_conversion(ELASTANCE_KEY, 1.0);
	this->set_conversion(RESISTANCE_KEY, 1.0);
	this->set_conversion(INERTANCE_KEY, 1.0);
	this->set_conversion(VISCOSITY_KEY, 1.0);
	this->set_conversion(DENSITY_KEY, 1.0);
	this->set_conversion(DIFFUSIVITY_KEY, 1.0);
	this->set_conversion(FORCE_DENSITY_KEY, 1.0);
    this->set_conversion(ACCELERATION_KEY, 1.0);
    this->set_conversion(VOL_FLUX_KEY,1.0);

	// volume params
	this->add(FRC_VOLUME_KEY, std::make_shared<VolumeParam>(FRC_VOLUME_DEFAULT_LITRES, 
		                                    std::string(FRC_VOLUME_KEY), this->get_conversion_ptr(VOL_KEY)));
	this->add(AIRWAY_DEAD_SPACE_KEY, std::make_shared<VolumeParam>(AIRWAY_DEAD_SPACE_DEFAULT_LITRES,
		                                           std::string(AIRWAY_DEAD_SPACE_KEY), this->get_conversion_ptr(VOL_KEY)));
	this->add(TIDAL_VOLUME_KEY, std::make_shared<VolumeParam>(TIDAL_VOLUME_DEFAULT_LITRES, 
											  std::string(TIDAL_VOLUME_KEY), this->get_conversion_ptr(VOL_KEY)));
	this->add(PLEURAL_PRESSURE_CHANGE_KEY, std::make_shared<PressureParam>(PLEURAL_PRESSURE_CHANGE_DEFAULT_CMH20, 
											  std::string(PLEURAL_PRESSURE_CHANGE_KEY), this->get_conversion_ptr(PRESSURE_KEY)));											  
	this->add(MOUTH_VOLUME_KEY, std::make_shared<VolumeParam>(MOUTH_VOLUME_DEFAULT_LITRES, 
											  std::string(MOUTH_VOLUME_KEY), this->get_conversion_ptr(VOL_KEY)));
	this->add(BOLUS_START_INSP_KEY, std::make_shared<VolumeParam>(BOLUS_START_VOL_DEFAULT_LITRES,
		std::string(BOLUS_START_INSP_KEY), this->get_conversion_ptr(VOL_KEY)));
	this->add(BOLUS_VOL_DURATION_KEY, std::make_shared<VolumeParam>(BOLUS_DURATION_VOL_DEFAULT_LITRES,
		std::string(BOLUS_VOL_DURATION_KEY), this->get_conversion_ptr(VOL_KEY)));

	//time
	this->add(TIDAL_TIME_KEY, std::make_shared<TimeParam>(TIDAL_TIME_DEFAULT_SECONDS, 
											  std::string(TIDAL_TIME_KEY), this->get_conversion_ptr(TIME_KEY)));
	this->add(TIME_STEP_KEY, std::make_shared<TimeParam>(TIME_STEP_DEFAULT_SECONDS, 
		                                      std::string(TIME_STEP_KEY), this->get_conversion_ptr(TIME_KEY)));
	this->add(WASHIN_END_KEY, std::make_shared<TimeParam>(WASHIN_END_DEFAULT_SECONDS, 
		                                      std::string(WASHIN_END_KEY), this->get_conversion_ptr(TIME_KEY)));
	this->add(PRINT_TIME_KEY, std::make_shared<TimeParam>(PRINT_TIME_DEFAULT_SECONDS, 
		                                      std::string(PRINT_TIME_KEY), this->get_conversion_ptr(TIME_KEY)));
	
	//pressure per length
	this->add(GRAVITY_PRESSURE_GRADIENT_KEY, std::make_shared<PressureParam>(GRAVITY_PRESSURE_GRADIENT_DEFAULT_CMH20_PER_MM,
		                                     std::string(GRAVITY_PRESSURE_GRADIENT_KEY), 
											 this->get_conversion_ptr(FORCE_DENSITY_KEY)));

	//elastance
	this->add(LUNG_ELASTANCE_KEY, std::make_shared<ElastanceParam>(LUNG_ELASTANCE_DEFAULT_CMH20_PER_L, 
											       std::string(LUNG_ELASTANCE_KEY), this->get_conversion_ptr(ELASTANCE_KEY)));
	
	//resistance
	this->add(BAG_RESISTANCE_KEY, std::make_shared<ResistanceParam>(BAG_RESISTANCE_DEFAULT_CMH20S_PER_L, 
											       std::string(BAG_RESISTANCE_KEY), this->get_conversion_ptr(RESISTANCE_KEY)));
	this->add(MOUTH_RESISTANCE_KEY, std::make_shared<ResistanceParam>(MOUTH_RESISTANCE_DEFAULT_CMH20S_PER_L, 
											       std::string(MOUTH_RESISTANCE_KEY), this->get_conversion_ptr(RESISTANCE_KEY)));
	
	//inertance
	this->add(BAG_INERTANCE_KEY, std::make_shared<InertanceParam>(BAG_INERTANCE_DEFAULT_CMH20S2_PER_L,
		                                            std::string(BAG_INERTANCE_KEY), this->get_conversion_ptr(INERTANCE_KEY)));

	//viscosity
	this->add(AIR_VISCOSITY_KEY, std::make_shared<ViscosityParam>(AIR_VISCOSITY_DEFAULT_CMH20S, 
												   std::string(AIR_VISCOSITY_KEY), this->get_conversion_ptr(VISCOSITY_KEY)));
	
	//density
	this->add(AIR_DENSITY_KEY, std::make_shared<DensityParam>(AIR_DENSITY_DEFAULT_KG_PER_M3, 
												   std::string(AIR_DENSITY_KEY), this->get_conversion_ptr(DENSITY_KEY)));

    //density
    this->add(PARTICLE_DENSITY_KEY, std::make_shared<DensityParam>(PARTICLE_DENSITY_DEFAULT_KG_PER_M3,
                                                 std::string(PARTICLE_DENSITY_KEY), this->get_conversion_ptr(DENSITY_KEY)));

	//length
	this->add(MAX_FV_LENGTH_KEY, std::make_shared<LengthParam>(MAX_FV_LENGTH_DEFAULT_MM, std::string(MAX_FV_LENGTH_KEY), 
		                                       this->get_conversion_ptr(LENGTH_KEY)));

    this->add(PARTICLE_SIZE_KEY, std::make_shared<LengthParam>(PARTICLE_SIZE_DEFAULT_MM, std::string(PARTICLE_SIZE_KEY),
                                                               this->get_conversion_ptr(LENGTH_KEY)));
	//diffusivity
	this->add(GAS_DIFFUSIVITY_KEY, std::make_shared<DiffusivityParam>(GAS_DIFFUSIVITY_DEFAULT_CM2_PER_S, std::string(GAS_DIFFUSIVITY_KEY), 
													  this->get_conversion_ptr(DIFFUSIVITY_KEY)));

    //grav accel
    this->add(GRAV_ACCEL_KEY, std::make_shared<AccelerationParam>(GRAV_ACCEL_MS2, std::string(GRAV_ACCEL_KEY),
                                                                      this->get_conversion_ptr(ACCELERATION_KEY)));

    //Volume flux for deposition diffusivity: Boltzmann constant * Air temperature / Air viscosity
    this->add(BOLTZMANN_TEMP_VISCOSITY_KEY, std::make_shared<VolumeFluxParam>(BOLTZMANN_TEMP_VISCOSITY_DEFAULT_M3,std::string(BOLTZMANN_TEMP_VISCOSITY_KEY),
                                                                          this->get_conversion_ptr(VOL_FLUX_KEY)));
    //Mean free path
    this->add(MEAN_FREE_PATH_KEY,std::make_shared<LengthParam>(MEAN_FREE_PATH_DEFAULT_MM,std::string(MEAN_FREE_PATH_KEY),
                                                               this->get_conversion_ptr(LENGTH_KEY)));

	//ratios (between 0 and 1)
	this->add(ACIN_ASYMM_FACTOR_KEY, std::make_shared<RatioParam>(ACIN_ASYMM_FACTOR_DEFAULT, std::string(ACIN_ASYMM_FACTOR_KEY)));
	this->add(ACIN_SCALE_FACTOR_KEY, std::make_shared<RatioParam>(ACIN_SCALE_FACTOR_DEFAULT, std::string(ACIN_SCALE_FACTOR_KEY)));
	
	//non-dimensional
	this->add(ACIN_LR_RATIO_KEY, std::make_shared<Parameter<double>>(ACIN_LR_RATIO_DEFAULT, std::string(ACIN_LR_RATIO_KEY)));
	this->add(BREATH_COUNT_KEY, std::make_shared<Parameter<double>>(BREATH_COUNT_DEFAULT, std::string(BREATH_COUNT_KEY)));
	this->add(MAX_ELEMENT_PE_KEY, std::make_shared<Parameter<double>>(MAX_ELEMENT_PE_DEFAULT, std::string(MAX_ELEMENT_PE_KEY)));
}

void SimParameterList::check_and_convert(SimOptionList *o)
{
	this->check_OK();

	//check for contradictory parameters -- HERE
	if(!o->get_option<bool>(SIMULATE_WASHIN_KEY)->get_value())
	{
		dict2[WASHIN_END_KEY]->update_value(-1.0);
	}
	
	if (dict2[AIRWAY_DEAD_SPACE_KEY]->get_phys_value() > dict2[FRC_VOLUME_KEY]->get_phys_value()) {
        std::cerr << "Error, airway dead-space large than FRC.\n";
        dict2[AIRWAY_DEAD_SPACE_KEY]->set_to_default();
        if (dict2[AIRWAY_DEAD_SPACE_KEY]->get_phys_value() > dict2[FRC_VOLUME_KEY]->get_phys_value()) {
            dict2[AIRWAY_DEAD_SPACE_KEY]->set_phys_value(0.1 * dict2[FRC_VOLUME_KEY]->get_phys_value());
            std::cerr << "Setting deadspace to 10% of FRC value: " << dict2[AIRWAY_DEAD_SPACE_KEY]->phys_value_string()
                      << "L.\n";
        } else {
            std::cerr << "Setting deadspace to default value: " << dict2[AIRWAY_DEAD_SPACE_KEY]->phys_value_string()
                      << ".\n";
        }
    }

    //check diffusivity option -- calculate in cm2/s and change the phys value to particle rather than gas diff
    if (o->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value())
    {
        //do calculation in SI units (m^2/s)
        double particle_diam_m = dict2[PARTICLE_SIZE_KEY]->get_SI_value();
        std::cout << dict2[PARTICLE_SIZE_KEY]->get_SI_value() << std::endl;
		//std::cout << "Particle diam m: " << particle_diam_m << std::endl;
        double mean_free_path_m = this->get_param<double>(MEAN_FREE_PATH_KEY)->get_SI_value();        //needs to be calculated
        //std::cout << "Mean free path m: " << mean_free_path_m << std::endl;
        double cslip_factor = 1 + (2*mean_free_path_m/particle_diam_m)*(1.257 + 0.4*exp(-0.55*particle_diam_m/mean_free_path_m));
		//std::cout << "C_slip factor: " << cslip_factor << std::endl;
        double particle_diff_ms2 = cslip_factor * BOLTZMANN_CONSTANT * BODY_TEMP_K / (3 * M_PI * particle_diam_m * dict2[AIR_VISCOSITY_KEY]->get_SI_value());
		//std::cout << "Particle diff m2s: " << particle_diff_ms2 << std::endl;
        //set value in Phys units (cm^2/s)
        this->dict2[GAS_DIFFUSIVITY_KEY]->set_phys_value(particle_diff_ms2*10000);

		if (o->get_option<char>(DIFFUSIVITY_TYPE_KEY)->get_value() == TAYLOR_DISPERSION) {
			std::cerr << "Error, Taylor dispersion leads to unphysical behaviour for aerosols, " <<
				"changing to modified Scherer diffusivity." << std::endl;
			o->get_option<char>(DIFFUSIVITY_TYPE_KEY)->update_value(SCHERER_DISPERSION);
		}
    }

	//check that bolus is within tidal volume range (if tidal volume is being used)
	if ((o->get_option<char>(INHALATION_MODE_KEY)->get_value() == BOLUS_CONC) && o->get_option<bool>(FLOW_BC_OPTION_KEY)->get_value() && o->get_option<char>(FLOW_INPUT_KEY)->get_value() != FLOW_FROM_FILE)
	{
		double VT = dict2[TIDAL_VOLUME_KEY]->get_phys_value();
		double VB0 = dict2[BOLUS_START_INSP_KEY]->get_phys_value();
		double DVB = dict2[BOLUS_VOL_DURATION_KEY]->get_phys_value();
		if (VB0 + DVB > VT)
		{
			std::cout << "Warning end of bolus volume larger than tidal volume, changing to maximum value:" << std::endl;
			double DVBnew = std::min(DVB,VT);
			double VB0new = VT - DVBnew;
			std::cout << dict2[BOLUS_START_INSP_KEY]->get_name() << " = " << VB0new << " and " 
					  << dict2[BOLUS_VOL_DURATION_KEY]->get_name() << " = " << DVBnew << "." << std::endl;
			dict2[BOLUS_START_INSP_KEY]->set_phys_value(VB0new);
			dict2[BOLUS_VOL_DURATION_KEY]->set_phys_value(DVBnew);
		}
	}

	//normalising conversions
	this->set_conversion(VOL_KEY, 1.0/this->dict2[FRC_VOLUME_KEY]->get_phys_value());   //lung volume set to 1
	this->set_conversion(TIME_KEY, 1.0/this->dict2[TIDAL_TIME_KEY]->get_phys_value());  //breath time set to 1
	this->set_conversion(PRESSURE_KEY, 1.0/(this->dict2[LUNG_ELASTANCE_KEY]->get_phys_value()*this->dict2[FRC_VOLUME_KEY]->get_phys_value()));   //lung volume set to 1
	//simulation values are defined so that FRC volume is 1.0, breath cycle time is 1.0, and pleural pressure is 1.0
    VolumeParam V1(1.0, std::string("scale"), this->get_conversion_ptr(VOL_KEY));
	TimeParam T1(1.0, std::string("scale"), this->get_conversion_ptr(TIME_KEY));
	PressureParam P1(1.0, std::string("scale"), this->get_conversion_ptr(PRESSURE_KEY));
    //convert everything else accordingly by working in SI
	double Vscale = V1.get_value() / V1.get_SI_value();
	double Tscale = T1.get_value() / T1.get_SI_value();
	double Pscale = P1.get_value() / P1.get_SI_value();

	//derived conversions -- calc from si units, as these match up
	//length = volume^(1/3)
	LengthParam L1(1.0, std::string("scale"));
	this->set_conversion(LENGTH_KEY,L1.get_SI_value() * pow( Vscale , 1.0/3.0));
	L1.set_conversion(this->get_conversion_ptr(LENGTH_KEY));  //use this in future definitions too
	double Lscale = L1.get_value() / L1.get_SI_value();

	//pressure per length
	PressurePerLengthParam PL1(1.0, std::string("scale"));
	this->set_conversion(FORCE_DENSITY_KEY, PL1.get_SI_value() * Pscale / Lscale);

	//reciprocal length
	ReciprocalLengthParam RL1(1.0, std::string("scale"));
	this->set_conversion(RECIPROCAL_LENGTH_KEY, RL1.get_SI_value() / Lscale);

	//density = pressure * time^2 / length^2
	DensityParam D1(1.0, std::string("scale"));
	this->set_conversion(DENSITY_KEY,D1.get_SI_value() * Pscale * Tscale * Tscale   / (Lscale * Lscale));

	//elastance = pressure / volume
	ElastanceParam E1(1.0, std::string("scale"));
	this->set_conversion(ELASTANCE_KEY,E1.get_SI_value() * Pscale / Vscale);

	//resistance = pressure * time / volume
	ResistanceParam R1(1.0, std::string("scale"));
	this->set_conversion(RESISTANCE_KEY,R1.get_SI_value() * Pscale * Tscale / Vscale);

	//inertance = pressure * time^2 / volume
	InertanceParam I1(1.0, std::string("scale"));
	this->set_conversion(INERTANCE_KEY, I1.get_SI_value() * Pscale * Tscale * Tscale / Vscale);

	//viscosity = pressure * time
	ViscosityParam Visc1(1.0, std::string("scale"));
	this->set_conversion(VISCOSITY_KEY,Visc1.get_SI_value() * Pscale * Tscale);

	//diffusivity = lenght^2 / time
	DiffusivityParam Diff1(1.0, std::string("scale"));
	this->set_conversion(DIFFUSIVITY_KEY,Diff1.get_SI_value() * Lscale * Lscale / Tscale);

    //acceleration
    AccelerationParam Acc1(1.0, std::string("scale"));
    this->set_conversion(ACCELERATION_KEY,Acc1.get_SI_value() * Lscale / (Tscale * Tscale));

    //volume flux = length^3/time
    VolumeFluxParam Flux1(1.0,std::string("scale"));
    this->set_conversion(VOL_FLUX_KEY,Flux1.get_SI_value() * Vscale / Tscale);

	//convert and get simulation values
	for (auto it = dict2.begin(); it!=dict2.end(); ++it)
	{
		//update sim values due to change in conversions
        it->second->calc_sim_from_phys_value();
	}
}

void VolumeParam::define_units()
{
	this->PhysUnits = "L";
	this->SIunits = "m^3";
}

void VolumeParam::calc_SI_from_Phys_value()
{ 
	this->SIValue = this->get_phys_value() * L_TO_M3; 
}

void TimeParam::define_units()
{
	this->SIunits = "s";
	this->PhysUnits = "s";
}

void TimeParam::calc_SI_from_Phys_value()
{ 
	this->SIValue = this->get_phys_value();
}

void LengthParam::define_units()
{
	this->SIunits = "m";
	this->PhysUnits = "mm";
}

void LengthParam::calc_SI_from_Phys_value()
{
	this->SIValue = MM_TO_M * this->get_phys_value();
}

void ReciprocalLengthParam::define_units()
{
	this->SIunits = "m^(-1)";
	this->PhysUnits = "mm^(-1)";
}

void ReciprocalLengthParam::calc_SI_from_Phys_value()
{ 
	this->SIValue = this->get_phys_value() / MM_TO_M; 
}

void PressureParam::define_units()
{
	this->SIunits = "Pa";
	this->PhysUnits = "cmH20";
}

void PressureParam::calc_SI_from_Phys_value()
{
	this->SIValue = CMH20_TO_PA * this->get_phys_value(); 
}

void PressurePerLengthParam::define_units()
{
	this->SIunits = "Pa m^(-1)";
	this->PhysUnits = "cmH20 mm^(-1)";
}

void PressurePerLengthParam::calc_SI_from_Phys_value()
{
	this->SIValue = 1000 * CMH20_TO_PA * this->get_phys_value(); 
}

void ElastanceParam::define_units()
{
	this->SIunits = "Pa/m^3";
	this->PhysUnits = "cmH20/L";
}

void ElastanceParam::calc_SI_from_Phys_value()
{ 
	this->SIValue = ( CMH20_TO_PA / L_TO_M3 ) * this->get_phys_value(); 
}

void ResistanceParam::define_units()
{
	this->SIunits = "Pas/m^3";
	this->PhysUnits = "cmH20s/L";
}

void ResistanceParam::calc_SI_from_Phys_value()
{
	this->SIValue = ( CMH20_TO_PA / L_TO_M3 ) * this->get_phys_value(); 
}

void InertanceParam::define_units()
{
	this->SIunits = "Pas^2/m^3";
	this->PhysUnits = "cmH20s^2/L";
}

void InertanceParam::calc_SI_from_Phys_value()
{
	this->SIValue = ( CMH20_TO_PA / L_TO_M3 ) * this->get_phys_value(); 
}

void ViscosityParam::define_units()
{
	this->SIunits = "Pas";
	this->PhysUnits = "cmH20s";
}

void ViscosityParam::calc_SI_from_Phys_value()
{
	this->SIValue = ( CMH20_TO_PA ) * this->get_phys_value(); 
}

void DiffusivityParam::define_units()  
{
	this->SIunits = "m2/s";
	this->PhysUnits = "cm2/s";
}

void DiffusivityParam::calc_SI_from_Phys_value()
{
	this->SIValue = 0.0001 * this->get_phys_value(); 
}

void AccelerationParam::define_units()
{
    this->SIunits = "m/s^2";
    this->PhysUnits = "m/s^2";
}

void AccelerationParam::calc_SI_from_Phys_value()
{
    this->SIValue = this->get_phys_value();
}

void DensityParam::define_units()
{
    this->SIunits = "kg/m^3";
    this->PhysUnits = "kg/m^3";
}

void DensityParam::calc_SI_from_Phys_value()
{
    this->SIValue = this->get_phys_value();
}

void VolumeFluxParam::define_units()
{
    this->SIunits = "m^3/s";
    this->PhysUnits = "m^3/s";
}

void VolumeFluxParam::calc_SI_from_Phys_value()
{
    this->SIValue = this->get_phys_value();
}

