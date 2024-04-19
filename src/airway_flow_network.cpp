#include "airway_flow_network.h"
#include "acinus_models.h"
#include "define_codes.h"
#include "define_options.h"
#include "define_params.h"
#include "globals.h"
#include<iostream>
#include<chrono>

AirwayFlowNetwork::AirwayFlowNetwork(std::map<long int, std::shared_ptr<AirwayNetworkNode>> & node, 
									 std::map<long int, std::shared_ptr<AirwayNetworkEdge>> & edge,
									 SimOptionList *o, SimParameterList *p):AirwayNetwork(node, edge)
{
	//assign parameters to network
	this->setup_afn(o,p);
}

void AirwayFlowNetwork::setup_afn(SimOptionList *o, SimParameterList *p)
{
	this->init_var_vectors();   //initialise blank variable lists
	this->assign_params_to_flow_network(o,p);
	std::cout << "Creating pv solver...\n";
	this->pv_solver = create_pv_solver(o, p);
	std::cout << "PV solver set up. \n";
}

std::shared_ptr<PressureVolumeSolver> AirwayFlowNetwork::create_pv_solver(SimOptionList * o, SimParameterList *p)
{
	using namespace inlist;

	bool is_linear;
	std::shared_ptr<PressureVolumeSolver> rv;
	// Decide whether model is linear and no feedback
	if(o->get_option<char>(FLOW_TYPE_KEY)->get_value() == POISEUILLE && o->get_option<char>(ELASTICITY_MODEL_KEY)->get_value() == LINEAR_ELASTICITY)
	{
		is_linear = true;
	}
	else is_linear = false;
	// Is flow also sinusoidal?
	switch(o->get_option<char>(FLOW_INPUT_KEY)->get_value())
	{
	case SIN_FLOW:
		{
			if(is_linear)
			{
				double bag_vol = p->get_param<double>(FRC_VOLUME_KEY)->get_value() - p->get_param<double>(AIRWAY_DEAD_SPACE_KEY)->get_value();
				rv = std::make_shared<LinearSinusoidalPressureVolumeSolver>
					                     (this, p->get_param<double>(TIDAL_VOLUME_KEY)->get_value(),
										  p->get_param<double>(PLEURAL_PRESSURE_CHANGE_KEY)->get_value(),
			                              p->get_param<double>(TIDAL_TIME_KEY)->get_value(),
										  bag_vol, p->get_param<double>(TIME_STEP_KEY)->get_value());
			}
			else
			{
				rv = std::make_shared<GenericPressureVolumeSolver>
					                      (this, p->get_param<double>(TIDAL_VOLUME_KEY)->get_value(), 
					                       p->get_param<double>(TIDAL_TIME_KEY)->get_value(),
										   p->get_param<double>(TIME_STEP_KEY)->get_value(),
										   SIN_FLOW);
			}
		} break;

	case FLOW_FROM_FILE:
		{
			
			std::vector<double> flux_in = read_flow_file(o->get_filename(FLOW_INPUT_FILE_EXT), p->get_param<double>(TIME_STEP_KEY)->get_SI_value(),
														 p->get_conversion(TIME_KEY), p->get_conversion(VOL_KEY),
														 p->get_param<double>(BREATH_COUNT_KEY));
			if(is_linear)
			{
				rv = std::make_shared<LinearGenericPressureVolumeSolver>(this, 
					                  p->get_param<double>(TIME_STEP_KEY)->get_value(), flux_in);
			}
			else
			{
				rv = std::make_shared<GenericPressureVolumeSolver>(this, 
					                  p->get_param<double>(TIME_STEP_KEY)->get_value(), flux_in);
			}
		} break;

	default:
		{
			if(is_linear)
			{
				rv = std::make_shared<LinearGenericPressureVolumeSolver>
					             (this, p->get_param<double>(TIDAL_VOLUME_KEY)->get_value(), 
					              p->get_param<double>(TIDAL_TIME_KEY)->get_value(),
								  p->get_param<double>(TIME_STEP_KEY)->get_value(),
								  o->get_option<char>(FLOW_INPUT_KEY)->get_value());
			}
			else
			{
				rv = std::make_shared<GenericPressureVolumeSolver>
					            (this, p->get_param<double>(TIDAL_VOLUME_KEY)->get_value(), 
					             p->get_param<double>(TIDAL_TIME_KEY)->get_value(),
								 p->get_param<double>(TIME_STEP_KEY)->get_value(),
								 o->get_option<char>(FLOW_INPUT_KEY)->get_value());
			}
		} break;
	}
	auto start = std::chrono::system_clock::now();
	rv->setup_flow_calc(o->get_option<bool>(FLOW_BC_OPTION_KEY)->get_value());
	auto end = std::chrono::system_clock::now();
	std::cout << "Flow initialisation took: " << (std::chrono::duration<double>(end-start)).count() << "s.\n"; 
	return rv;
}

void AirwayFlowNetwork::assign_params_to_flow_network(SimOptionList *o, SimParameterList *p)
{
	using namespace inlist;

	double Nt = 0;

	for (size_t k = 0; k < this->count_nodes(); k++)
	{
		this->node_vals[FN_NODE_PRESSURE][k] = 0;    //initialise
		if(k >= this->term_start)    //term nodes
		{
			Nt +=  this->get_node(k)->get_original_npts();
		}
	}

	//calculate airway dead space
	double AirwayDeadSpace = 0;
	for (size_t j = 0; j < this->count_edges(); j++)
	{
		this->edge_vals[FN_EDGE_FLUX][j] = 0;    //initialise
		AirwayDeadSpace += this->get_edge(j)->get_inner_volume();
	}
	//update parameter
	p->get_param<double>(AIRWAY_DEAD_SPACE_KEY)->update_value(AirwayDeadSpace);
	//std::cout << "Dead space L: " << p->get_param<double>(AIRWAY_DEAD_SPACE_KEY)->get_phys_value() << '\n';

	//mouth volume does not count as part of FRC..
	double bag_vol = p->get_param<double>(FRC_VOLUME_KEY)->get_value() - p->get_param<double>(AIRWAY_DEAD_SPACE_KEY)->get_value();
	//std::cout << "Bag vol L: " << bag_vol/p->get_conversion(VOL_KEY) << '\n';
	this->e_calc = parse_elast_option(o->get_option<char>(ELASTICITY_MODEL_KEY)->get_value());

	//gravity gradients
	network::Position grav_dir;
	std::vector<double> grav_dist;
	this->setup_gravity_dist(o->get_option<char>(GRAVITY_MODEL_KEY)->get_value(),
		                     o->get_grav_direction(), 
							 p->get_param<double>(GRAVITY_PRESSURE_GRADIENT_KEY)->get_value(), Nt,
							 grav_dist);

	std::vector<double> intrinsic_acinus_vols;
	this->setup_acinus_vols_dist(o->get_option<bool>(RAD_DEPENDENT_ACIN_VOL_KEY)->get_value(),
		                         bag_vol, Nt, intrinsic_acinus_vols);


	if (o->get_option<char>(ELASTICITY_MODEL_KEY)->get_value() == LINEAR_ELASTICITY){
		//if linear, elastance won't depend on extension, and so it is more accurate to assume a compliance gradient
		std::shared_ptr<ElastanceObject> tawhai_elast = parse_elast_option(TAWHAI_ELASTICITY);
		double K0 = p->get_param<double>(LUNG_ELASTANCE_KEY)->get_value() * this->count_term_nodes();
		for (size_t k = this->term_start; k < this->NodeVec.size(); k++)   //loop over terminal nodes
		{
			size_t kt = k - this->term_start;   //term index
			double Knew = K0;
			double V0 = intrinsic_acinus_vols[kt];
			double vol_old = V0;
			double vol_new = V0 - grav_dist[kt] / K0;
			while (abs((vol_new - vol_old) / V0) > 1E-03)
			{
				vol_old = vol_new;
				Knew = tawhai_elast->calc_elastance(pow(2.0 * vol_new / V0, 1.0 / 3.0), K0);
				vol_new = V0 - grav_dist[kt] / Knew;
			}
			//std::cout << "z = " << this->get_node(k)->get_pos().x[2] << " Pgrad = " << grav_dist[kt] << " DV = " << vol_new/V0 << " DK = " << Knew/K0 << std::endl;
			this->term_node_vals[FN_TNODE_ELASTANCE][kt] = Knew;  //use this as initial elastance
		}
	}
	else {
		double K0 = p->get_param<double>(LUNG_ELASTANCE_KEY)->get_value() * this->count_term_nodes();
		for (size_t k = this->term_start; k < this->NodeVec.size(); k++)   //loop over terminal nodes
		{
			size_t kt = k - this->term_start;   //term index
			this->term_node_vals[FN_TNODE_ELASTANCE][kt] = K0;
		}
	}

	for (size_t k = this->term_start; k < this->NodeVec.size(); k++)   //loop over terminal nodes
	{
		size_t kt = k-this->term_start;   //term index

		double Nth;
		Nth =  this->get_node(k)->get_original_npts();

		//elastance already defined in gravity function
		this->term_node_vals[FN_TNODE_ELASTANCE][kt] = (bag_vol / (this->count_term_nodes()*intrinsic_acinus_vols[kt])) * 
														this->term_node_vals[FN_TNODE_ELASTANCE][kt];
		this->term_node_vals[FN_TNODE_RESISTANCE][kt] = (bag_vol / intrinsic_acinus_vols[kt]) * 
			                                          p->get_param<double>(BAG_RESISTANCE_KEY)->get_value();
		this->term_node_vals[FN_TNODE_INERTANCE][kt] = (bag_vol / intrinsic_acinus_vols[kt]) * 
			                                          p->get_param<double>(BAG_INERTANCE_KEY)->get_value();
		this->term_node_vals[FN_TNODE_PRESSGRAD][kt] = grav_dist[kt];   //positive if in direction of gravity relative to centre of mass
		this->term_node_vals[FN_TNODE_INTRINSIC_VOL][kt] = intrinsic_acinus_vols[kt]; 

		this->term_node_vals[FN_TNODE_VOLUME][kt] =  intrinsic_acinus_vols[kt]
		     - this->term_node_vals[FN_TNODE_PRESSGRAD][kt]/this->term_node_vals[FN_TNODE_ELASTANCE][kt];
	
		//std::cout << "Termnode: " << kt << '\n';
		//std::cout << "Elastance = " << this->term_node_vals[FN_TNODE_ELASTANCE][kt] << '\n';
		//std::cout << "Intrinsic Vol per endpt = " << this->term_node_vals[FN_TNODE_INTRINSIC_VOL][kt]/this->get_node(k)->get_original_npts() << '\n';
		//std::cout << "Vol per endpt = " << this->term_node_vals[FN_TNODE_VOLUME][kt]/this->get_node(k)->get_original_npts() << "\n\n";

	}
	this->Ppl = -this->e_calc->calc_elastance(pow(2.0,1.0/3.0), p->get_param<double>(LUNG_ELASTANCE_KEY)->get_value());
	//std::cout<< "Total intrinsic vol: " << this->term_node_vals[FN_TNODE_INTRINSIC_VOL].sum() / p->get_conversion(VOL_KEY) << '\n';
	//std::cout<< "Total vol: " << this->term_node_vals[FN_TNODE_VOLUME].sum() / p->get_conversion(VOL_KEY) << '\n';
	//mouth resistance
	this->mouth_resistance = p->get_param<double>(MOUTH_RESISTANCE_KEY)->get_value();
	//same resistance object for all
	this->res_calc = parse_res_option(o->get_option<char>(FLOW_TYPE_KEY)->get_value(), p->get_param<double>(AIR_VISCOSITY_KEY)->get_value(), 
		                              p->get_param<double>(AIR_DENSITY_KEY)->get_value());

	this->update_all_edge_resistances_by_rule();
	std::cout << "Resistance calcs applied.\n";
}

std::vector<double> read_flow_file(const std::string & filename, const double & dt, 
								   const double & sec_to_simtime, const double & L_to_simvol,
								   inlist::Parameter<double> *stop_time)
{
	//need to read file and interpolate to correct time-step
	using namespace std;
	vector<double> flux_interp;

	if(check_infile(filename))
	{
		vector<vector<double>> data_in = parse_csv_file<double>(filename);
		vector<double> time, flux;   //flux in litres per second
		time.resize(data_in.size());
		flux.resize(data_in.size());
		for(size_t irow=0; irow < data_in.size(); irow++)
		{
			if(irow == 0) time[irow] = 0;
			else time[irow] = sec_to_simtime * (data_in[irow][0] - data_in[0][0]);
			flux[irow] = L_to_simvol * data_in[irow][1] / sec_to_simtime;
		}

		
		flux_interp.resize(size_t( time[time.size()-1] / dt ));
		size_t i_time = 0;
		for(size_t n = 0; n < flux_interp.size(); n++)
		{
			while(i_time < time.size()-1 && n*dt > time[i_time+1])
			{
				i_time++;
			}
			flux_interp[n] = (flux[i_time]*(time[i_time+1] - n*dt) + flux[i_time+1]*(n*dt - time[i_time])) / (time[i_time+1] - time[i_time]);
		}
		stop_time->update_value(flux_interp.size() * dt);
	}
	else
	{
		cout << "Problem reading flux input file, aborting...\n";
		abort_on_failure();
	}

	return flux_interp;
}

int AirwayFlowNetwork::print_flow_vtk(const std::string & filename, SimParameterList *p)
{
	using namespace std;
	std::unordered_map<string, vector<double>> extra_vals;

	for(size_t n = 0; n < FN_NODE_COUNT; n++)
	{
		extra_vals[fn_node_strings[n]] = vector<double>();
		extra_vals[fn_node_strings[n]].resize(this->count_nodes());
	}
	for(size_t k = 0; k < this->count_nodes(); k++) 
	{
		extra_vals[fn_node_strings[FN_NODE_PRESSURE]][k] = this->node_vals[FN_NODE_PRESSURE][k] / p->get_conversion(PRESSURE_KEY);
	}
	for(size_t n = 0; n < FN_EDGE_COUNT; n++)
	{
		extra_vals[fn_edge_strings[n]] = vector<double>();
		extra_vals[fn_edge_strings[n]].resize(this->count_edges());
	}
	for(size_t j = 0; j < this->count_edges(); j++) 
	{
		extra_vals[fn_edge_strings[FN_EDGE_FLUX]][j] = this->edge_vals[FN_EDGE_FLUX][j] * p->get_conversion(TIME_KEY) / p->get_conversion(VOL_KEY);
		extra_vals[fn_edge_strings[FN_EDGE_RESISTANCE]][j] = this->edge_vals[FN_EDGE_RESISTANCE][j] / p->get_conversion(RESISTANCE_KEY);
	}

	return this->print_vtk(filename, 1.0 / p->get_conversion(LENGTH_KEY), extra_vals);
}

void AirwayFlowNetwork::setup_gravity_dist(const char & option, const network::Position & grav_dir,
										   const double & press_gradient,
										   const double & Ntot, std::vector<double> & grav_dist)
{
	double press_gradient_per_length = 0;
	grav_dist.resize(this->NodeVec.size() - this->term_start);
	std::fill(grav_dist.begin(), grav_dist.end(), 0.0);
	double mean_dist = 0;
	if(option == LINEAR_GRAVITY_APPROX)
	{
		//gravity direction
		grav_dist.resize(this->NodeVec.size() - this->term_start);
		for (size_t k = this->term_start; k < this->NodeVec.size(); k++)   //loop over terminal nodes
		{
			size_t kt = k-this->term_start;   //term index
			grav_dist[kt] = this->get_node(k)->get_original_pos().dot(grav_dir);
			mean_dist += this->get_node(k)->get_original_npts()*grav_dist[kt];
		}	
		mean_dist = mean_dist/Ntot;
		for (size_t k = this->term_start; k < this->NodeVec.size(); k++)   //loop over terminal nodes
		{
			size_t kt = k-this->term_start;   //term index
			grav_dist[kt] = press_gradient*(grav_dist[kt] - mean_dist);
		}
	}
}

void AirwayFlowNetwork::setup_acinus_vols_dist(const bool & option, const double & bag_vol, 
											   const double & Ntot,
											   std::vector<double> & acinus_vols)
{
	acinus_vols.resize(this->count_term_nodes());

	if(option)
	{
		//acinus volume is dependent on terminal radius -- use weibel estimates for
		//gradient and normalise
		double mean_rad =  0;
		for (size_t k = this->term_start; k < this->NodeVec.size(); k++)   //loop over terminal nodes
		{
			if(this->get_node(k)->has_seed_rad_changed())
			{
				mean_rad += this->get_node(k)->get_original_npts()*this->get_node(k)->get_seed_rad();
			}
			else
			{
				mean_rad += this->get_edge(this->get_edge_in_index(k,0))->get_geom()->get_inner_radius()*
					        this->get_node(k)->get_original_npts();
			}
			
		}
		mean_rad /= Ntot;
		//assume max is 3 * mean and min is (1/3.0) * mean
		double adjusted_mean = 0;
		double Weibel_corr = 0.54;  //rad proportional to this * ac vol (Haefeli Bleur 89)
		std::vector<double> sr;
		sr.resize(this->count_term_nodes());
		for (size_t k = this->term_start; k < this->NodeVec.size(); k++)   //loop over terminal nodes
		{
			size_t kt = k - this->term_start;
			double seed_rad;
			if(this->get_node(k)->has_seed_rad_changed())
			{

				seed_rad = this->get_node(k)->get_seed_rad();
			}
			else
			{
				seed_rad = this->get_edge(this->get_edge_in_index(k,0))->get_geom()->get_inner_radius();
			}
			if(seed_rad > 3 * Weibel_corr * mean_rad)
			{
				seed_rad = 3 * Weibel_corr * mean_rad;
			}
			if(seed_rad < (1.0/3.0) * Weibel_corr * mean_rad)
			{
				seed_rad = (1.0/3.0) * Weibel_corr * mean_rad;
			}
			adjusted_mean += seed_rad * this->get_node(k)->get_original_npts();
			sr[kt] = seed_rad;
		}
		adjusted_mean /= Ntot;

		for (size_t k = this->term_start; k < this->NodeVec.size(); k++)   //loop over terminal nodes
		{
			size_t kt = k - this->term_start;
			acinus_vols[kt] = (this->get_node(k)->get_original_npts()*bag_vol/Ntot)
				                           *(1 + (sr[kt] - adjusted_mean)/(Weibel_corr));
		}
	}
	else
	{
		for (size_t k = this->term_start; k < this->NodeVec.size(); k++)   //loop over terminal nodes
		{
			size_t kt = k - this->term_start;
			acinus_vols[kt] = (this->get_node(k)->get_original_npts()*bag_vol/Ntot);
		}
	}
}
