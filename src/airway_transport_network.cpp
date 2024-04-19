#include "airway_transport_network.h"
#include "list_template.h"
#include "define_codes.h"
#include "globals.h"
#include <Eigen/Sparse>
//#include <boost/math/quadrature/trapezoidal.hpp>

void AirwayTransportNetwork::init_var_vectors()
{
	this->node_vals.resize(TN_NODE_COUNT);
	for (size_t n = 0; n < TN_NODE_COUNT; n++){
		this->node_vals[n] = Eigen::VectorXd::Zero(this->count_nodes());
	}

	this->edge_vals.resize(TN_EDGE_COUNT);
	for (size_t n = 0; n < TN_EDGE_COUNT; n++){
		this->edge_vals[n] = Eigen::VectorXd::Zero(this->count_edges());
	}

	for (size_t j = 0; j < this->count_edges(); j++) {
		this->edge_vals[TN_EDGE_HORSFIELD_ORDER][j] = this->edge_horsfield_order[j];
	}
}

AirwayTransportNetwork::AirwayTransportNetwork(AirwayNetwork * flow_tree, const Eigen::VectorXd & intrinsic_vols,
                                               SimOptionList *o, SimParameterList *p, AirwayNetwork *acin_template) {
    double init_conc;

	std::cout << "Md at start of transport network build: " << p->get_param<double>(GAS_DIFFUSIVITY_KEY)->get_value() << std::endl;

    //set initial concentration from option
    if (o->get_option<bool>(SIMULATE_WASHIN_KEY)->get_value()) {
        init_conc = 0;    //0 if simulating washin
    } else {
        init_conc = 1;    //1 if only simulating washout
    }
	this->print_acinar_airways = o->get_option<bool>(PRINT_ACINAR_AIRWAYS_KEY)->get_value();
    this->washin_mode = true;
    if (p->get_param<double>(MOUTH_VOLUME_KEY)->get_value() > 0) {
        this->Mouth = std::make_shared<DSVolume>(std::make_shared<FlexibleVolumeElement>
                                                         (p->get_param<double>(MOUTH_VOLUME_KEY)->get_value(),
                                                          init_conc, init_conc));
    } else Mouth = NULL;
	std::shared_ptr<DiscretisationOperator> cond_disc, acin_disc;
	if (o->get_option<char>(DISCRETISATION_KEY)->get_value() == MAX_MIN_DISC){
		double dx = p->get_param<double>(MAX_FV_LENGTH_KEY)->get_value();
		int minfv = p->get_param<int>(MIN_FV_PER_GEN_KEY)->get_value();
		int acinminfv = p->get_param<int>(ACIN_FV_PER_GEN_KEY)->get_value();
		cond_disc = std::make_shared<MaxMinDiscretisationOperator>(minfv,dx);
		acin_disc = std::make_shared<MaxMinDiscretisationOperator>(acinminfv, dx);
	}
	else{
		double maxPe = p->get_param<double>(MAX_ELEMENT_PE_KEY)->get_value();
		cond_disc = std::make_shared<SchererPeDiscretisationOperator>(maxPe);
		acin_disc = std::make_shared<SchererPeDiscretisationOperator>(maxPe);
	}

	std::cout << "Discretising for transport." << std::endl;


	this->discretise_for_transport(*(flow_tree), cond_disc.get());

    this->conducting_nodes.resize(this->NodeVec.size());
    this->conducting_edges.resize(this->EdgeVec.size());


    this->conducting_diffusivity = parse_diff_option(o->get_option<char>(DIFFUSIVITY_TYPE_KEY)->get_value(),
                                                     p->get_param<double>(
                                                             GAS_DIFFUSIVITY_KEY)->get_value(),
													 o->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value());  //same diffusivity object for all conducting branches

    //acinus diffusivity rule
    this->acinar_diffusivity = parse_diff_option(TAULBEE_DISPERSION,
                                                 p->get_param<double>(GAS_DIFFUSIVITY_KEY)->get_value(),
												 o->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value());
	std::cout << "Creating acinus" << std::endl;
    double acin_template_vol = 0;
    if (o->get_option<bool>(ACINUS_FROM_FILE)->get_value()) acin_template_vol = acin_template->get_total_edge_volume();
    this->acinus_vector.resize(flow_tree->count_term_nodes());
    this->acin_node_map.resize(flow_tree->count_term_nodes());
    this->acin_edge_map.resize(flow_tree->count_term_nodes());
    for (size_t k = flow_tree->get_first_term_index(); k < flow_tree->count_nodes(); k++) {
        size_t kt = k - flow_tree->get_first_term_index();
        double seed_rad;
        if (flow_tree->get_node(k)->has_seed_rad_changed())   //if network has been averaged, seed rad for acinus may be different to actual rad
        {
            seed_rad = flow_tree->get_node(k)->get_seed_rad();
        } else    //otherwise use the (unconstricted) radius
        {
            seed_rad = flow_tree->get_edge(flow_tree->get_edge_in_index(k, 0))->get_original_geom()->get_inner_radius();
        }
        if (o->get_option<bool>(ACINUS_FROM_FILE)->get_value())  //import acinus from file
        {
            //scale factor based on intrinsic vol
            double scale_factor = pow(intrinsic_vols[kt] / (flow_tree->get_node(k)->point_count() * acin_template_vol),
                                      1.0 / 3.0);

            //scale factor based on radius
            //double mean_rad = 0;
            //for(size_t jo = 0; jo < acin_template->count_edges_out(0); jo++)
            //{
            //	mean_rad += acin_template->get_edge(acin_template->get_edge_out_index(0,jo))->get_geom()->get_inner_radius();
            //}
            //mean_rad /= acin_template->count_edges_out(0);
            //double scale_factor = seed_rad / mean_rad;
            //double acin_temp_inner_vol = acin_template->get_total_inner_edge_volume();
            //double Npt = flow_tree->get_node(k)->point_count();
            //std::cout << kt << ' ' << seed_rad << ' ' << mean_rad << ' ' << scale_factor*scale_factor*scale_factor*Npt*acin_temp_inner_vol
            //	      << ' ' << intrinsic_vols[kt] << "\n\n";
            //
            //if(scale_factor*scale_factor*scale_factor*Npt*acin_temp_inner_vol > 0.3*intrinsic_vols[kt])  //at most inner volume can be 30% the intrinsic outer volume
            //{
            //	scale_factor =  pow(0.3*intrinsic_vols[kt]/(Npt*acin_temp_inner_vol), 1.0/3.0);
            //}
            //if(scale_factor*scale_factor*scale_factor*Npt*acin_temp_inner_vol < 0.2*intrinsic_vols[kt])  //at least it can be 20%
            //{
            //	scale_factor =  pow(0.2*intrinsic_vols[kt]/(Npt*acin_temp_inner_vol), 1.0/3.0);
            //}

            this->acinus_vector[kt] = import_acin_template(o->get_option<char>(ACINUS_MODEL_KEY)->get_value(),
                                                           scale_factor,
                                                           acin_disc.get(),
                                                           flow_tree->get_node_smart_ptr(k), *(acin_template));
            this->acinus_vector[kt]->update_volume(intrinsic_vols[kt]);
        } else   //build acinus from scratch
        {
			this->acinus_vector[kt] = parse_acinus_option(o->get_option<char>(ACINUS_TYPE_KEY)->get_value(),
                                                          flow_tree->get_node_smart_ptr(k), seed_rad,
                                                          intrinsic_vols[kt],
                                                          o->get_option<char>(ACINUS_MODEL_KEY)->get_value(), 
														  acin_disc.get(), p);
        }

        //the entrance edges already point to end_node
        this->acin_edge_map[kt].resize(this->acinus_vector[kt]->tree->count_edges());
        for (size_t ja = 0; ja < this->acinus_vector[kt]->tree->count_edges(); ja++)   //loop over acinus edges
        {
            this->EdgeVec.push_back(this->acinus_vector[kt]->tree->get_edge_smart_ptr(
                    ja));  //these were set up when acinus was built -- should be transport edges
            this->acin_edge_map[kt][ja] = this->EdgeVec.size() - 1;
        }
        this->acin_node_map[kt].resize(this->acinus_vector[kt]->tree->count_nodes());
        this->acin_node_map[kt][0] = this->get_node_index(flow_tree->get_node(k));
        for (size_t ka = 1; ka < this->acinus_vector[kt]->tree->count_nodes(); ka++) //loop over acinus edges
        {
            this->NodeVec.push_back(this->acinus_vector[kt]->tree->get_node_smart_ptr(ka));
            this->acin_node_map[kt][ka] = this->NodeVec.size() - 1;
        }

    }
    this->update_node_edge_maps();

    std::vector<size_t> nodes_old_to_new, edges_old_to_new;
    if (this->reorder_network(nodes_old_to_new, edges_old_to_new)) abort_on_failure();   //sort out network
    for (size_t k = 0; k < this->conducting_nodes.size(); k++) {
        this->conducting_nodes[k] = nodes_old_to_new[k];
    }
    for (size_t j = 0; j < this->conducting_edges.size(); j++) {
        this->conducting_edges[j] = edges_old_to_new[j];
    }

    bool incidence_needs_update = false;
    for (size_t kt = 0; kt < flow_tree->count_term_nodes(); kt++) {
        if(incidence_needs_update == false && this->acinus_vector[kt]->tree->count_edges() > 0) {
            incidence_needs_update = true;
        }
        for (size_t ka = 0; ka < this->acinus_vector[kt]->tree->count_nodes(); ka++) //loop over acinus edges
        {
            this->acin_node_map[kt][ka] = nodes_old_to_new[this->acin_node_map[kt][ka]];
        }
        for (size_t ja = 0; ja < this->acinus_vector[kt]->tree->count_edges(); ja++) {
            this->acin_edge_map[kt][ja] = edges_old_to_new[this->acin_edge_map[kt][ja]];
        }

    }
    //initialisation
    if(incidence_needs_update) this->fill_incidence_matrix();
    this->init_var_vectors();

    for (size_t k = 0; k < this->count_nodes(); k++) {
        this->node_vals[TN_NODE_CONC][k] = init_conc;
    }

    this->update_conducting_diffusivities();
    this->update_acinar_diffusivities();
    this->update_node_volumes();
	this->update_inner_node_volumes();
    //TODO: if statement here?
    this->update_branching_angles();
    this->transport_solver = create_transport_solver('n', this,
                                                     o->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value(),
													 o->get_option<bool>(SIMULATE_ALVEOLAR_DEPOSITION_KEY)->get_value());
    if (o->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value()) {
        this->transport_solver->set_deposition_params(p->get_param<double>(PARTICLE_SIZE_KEY)->get_value(),
                                                      p->get_param<double>(PARTICLE_DENSITY_KEY)->get_value(),
                                                      p->get_param<double>(GRAV_ACCEL_KEY)->get_value(),
                                                      p->get_param<double>(BOLTZMANN_TEMP_VISCOSITY_KEY)->get_value(),
                                                      p->get_param<double>(AIR_VISCOSITY_KEY)->get_value(),
                                                      p->get_param<double>(MEAN_FREE_PATH_KEY)->get_value(),
                                                      p->get_param<double>(AIR_DENSITY_KEY)->get_value(),
                                                      o->get_option<char>(DEP_FORMULA_DIFF_KEY)->get_value(),
                                                      o->get_option<char>(DEP_FORMULA_SED_KEY)->get_value(),
                                                      o->get_option<char>(DEP_FORMULA_IMP_KEY)->get_value(),
                                                      o->get_option<bool>(IMPACTION_BOTH_DIRS_KEY)->get_value());
        //std::cout << this->get_edge(0)->get_geom()->get_inner_radius() << std::endl;
    }

    //assign to airways
    if (this->is_tree) {
        std::vector<size_t> next_edges;    //stores next set of downstream edges to loop over
        this->edge_airway_no.resize(this->count_edges());    //stores an airway number corresponding to each edge
        next_edges.resize(this->count_edges_out(0));   //start at inlet node
        size_t n_airways = 0;           //counts no. of airways
        for (size_t jo = 0; jo < this->count_edges_out(0); jo++) {
            n_airways++;
            size_t j = this->get_edge_out_index(0, jo);
            next_edges[jo] = j;
            this->edge_airway_no[j] = n_airways - 1;
        }

        while (next_edges.size() > 0) {    //loop through tree
            std::vector<size_t> new_next_edges;
            new_next_edges.reserve(2 * next_edges.size());
            for (size_t ji = 0; ji < next_edges.size(); ++ji) {
                size_t j_in = next_edges[ji];
                size_t i = this->get_node_out_index(j_in);
                //check for branching event (means creation of new airway)
                bool branching_event = false;
                if (this->count_edges_out(i) > 1) branching_event = true;   //either there are >1 edges downstream
                else {
                    if (this->count_edges_out(i) == 1) {
                        size_t j = this->get_edge_out_index(i, 0);  //or the number of branches changes
                        if (this->get_edge(j_in)->branch_count() != this->get_edge(j)->branch_count()) {
                            branching_event = true;
                        }
                    }
                }
                //add next generation of airways
                if (branching_event) { //create new airways
                    for (size_t jo = 0; jo < this->count_edges_out(i); ++jo) {
                        n_airways++;
                        size_t j = this->get_edge_out_index(i, jo);
                        new_next_edges.push_back(j);
                        this->edge_airway_no[j] = n_airways - 1;
                    }
                } else {   //add elements to existing airways
                    for (size_t jo = 0; jo < this->count_edges_out(i); ++jo) {
                        size_t j = this->get_edge_out_index(i, jo);
                        this->edge_airway_no[j] = this->edge_airway_no[j_in];
                        new_next_edges.push_back(j);
                    }
                }
            }
            next_edges = new_next_edges;
        }
        //add these to vector of vector
        this->airway_to_edge_map.resize(n_airways);
        for(size_t j = 0; j < this->count_edges(); ++j){
            this->airway_to_edge_map[this->edge_airway_no[j]].push_back(j);
        }
        airway_count = n_airways;//store total number of airways
        for (int J = 0; J < n_airways; ++J) {
            no_airway_edges.push_back(airway_to_edge_map[J].size());
        }
    }
	std::cout << "Number of elements in transport tree: " << this->count_nodes() << std::endl;
}

double AirwayTransportNetwork::get_tot_IG_volume(const bool & UseInnerVol)
{
	double vol = 0;
	if (UseInnerVol) {
		for (size_t k = 0; k < this->NodeVec.size(); k++)
		{
			vol += this->get_node_conc(k) * this->node_vals[TN_INNER_NODE_VOLUME][k];
		}
	}
	else {
		for (size_t k = 0; k < this->NodeVec.size(); k++)
		{
			vol += this->get_node_conc(k) * this->node_vals[TN_NODE_VOLUME][k];
		}
	}
	
	return vol;
}

double AirwayTransportNetwork::get_acinus_IG_volume(const size_t & kt, const bool & UseInnerVol)  //kt is terminal node index
{
	double v = 0;
	for(size_t ka = 0; ka < this->acinus_vector[kt]->tree->count_nodes(); ka++)
	{
		size_t k = this->acin_node_map[kt][ka];
		if (UseInnerVol) v += this->get_node_conc(k) * this->node_vals[TN_INNER_NODE_VOLUME][k];
		else v += this->get_node_conc(k) * this->node_vals[TN_NODE_VOLUME][k];
	}
	return v;
}

void AirwayTransportNetwork::compute_all_fluxes(const double & dt)
{
	for(size_t ho = 0; ho < this->count_horsfield_orders(); ho++)  //loop over horsfield orders
	{
		for(size_t ji = 0; ji < this->count_edges_in_horsfield_order(ho); ji++)    //loop over edges in each order
		{
			size_t j = this->get_edge_index_from_horsfield_order(ho,ji);     //index of this edge
			double tot_flux_out = 0;
			size_t k_out = this->get_node_out_index(j);
			for (size_t jo = 0; jo < this->count_edges_out(k_out); jo++)
			{
				size_t eo_index = this->get_edge_out_index(this->get_node_out_index(j),jo);
				tot_flux_out += this->edge_vals[TN_EDGE_FLUX][eo_index];
			}
			tot_flux_out += (this->node_vals[TN_NODE_VOLUME][k_out] - this->transport_solver->get_node_vol_old(k_out)) / dt;
			this->edge_vals[TN_EDGE_FLUX][j] = tot_flux_out;
			//if(tot_flux_out != tot_flux_out)
			//{
			//	std::cout << "Error, flux out = " << ho << ' ' << j << ' ' << tot_flux_out << std::endl;
			//	system("pause");
			//}
		}
	}
}

void AirwayTransportNetwork::update_node_volumes()  //outer volume
{
    this->node_vals[TN_NODE_VOLUME] = Eigen::VectorXd::Zero(this->count_nodes());
    for(size_t k = 0; k < this->count_nodes(); k++)
    {
        this->node_vals[TN_NODE_VOLUME][k] = this->get_node(k)->get_intrinsic_vol();
    }
    for(size_t j = 0; j < this->count_edges(); j++)  //add half of each edge volume to connected nodes
    {
        double half_vol = 0.5*this->get_edge(j)->branch_count()*this->get_edge(j)->get_geom()->outer_volume();
        this->node_vals[TN_NODE_VOLUME][this->get_node_in_index(j)] += half_vol;
        this->node_vals[TN_NODE_VOLUME][this->get_node_out_index(j)] += half_vol;
    }
}

void AirwayTransportNetwork::update_inner_node_volumes()  //outer volume
{
    this->node_vals[TN_INNER_NODE_VOLUME] = Eigen::VectorXd::Zero(this->count_nodes());
    for(size_t k = 0; k < this->count_nodes(); k++)
    {
        this->node_vals[TN_INNER_NODE_VOLUME][k] = this->get_node(k)->get_intrinsic_vol();
    }
    for(size_t j = 0; j < this->count_edges(); j++)  //add half of each edge volume to connected nodes
    {
        double half_vol = 0.5*this->get_edge(j)->branch_count()*this->get_edge(j)->get_geom()->inner_volume();
        this->node_vals[TN_INNER_NODE_VOLUME][this->get_node_in_index(j)] += half_vol;
        this->node_vals[TN_INNER_NODE_VOLUME][this->get_node_out_index(j)] += half_vol;
    }
}

void AirwayTransportNetwork::update_acinar_volume(const size_t &kt, const double &v)
{
	this->acinus_vector[kt]->update_volume(v);
    for (size_t kh = 0; kh < this->acin_node_map[kt].size(); kh++) {
        size_t k = this->acin_node_map[kt][kh];
        this->node_vals[TN_NODE_VOLUME][k] = this->get_node(k)->get_intrinsic_vol();
        for(size_t ji = 0; ji < this->count_edges_in(k); ji++)
        {
            size_t j = this->get_edge_in_index(k,ji);
            this->node_vals[TN_NODE_VOLUME][k] += 0.5*this->get_edge(j)->branch_count()*this->get_edge(j)->get_geom()->outer_volume();
        }
        for(size_t jo = 0; jo < this->count_edges_out(k); jo++)
        {
            size_t j = this->get_edge_out_index(k,jo);
            this->node_vals[TN_NODE_VOLUME][k] += 0.5*this->get_edge(j)->branch_count()*this->get_edge(j)->get_geom()->outer_volume();
        }
    }
}




void AirwayTransportNetwork::update_branching_angles()
{
    for(size_t j = 0; j < this->count_edges(); j++)
    {

        pos::Position<double> edge_dir = this->get_edge(j)->get_direction_vector();
        for(size_t jout = 0; jout < this->count_edges_out(this->get_node_out_index(j)); jout++)
        {
            pos::Position<double> edge_out_dir = this->get_edge(this->get_edge_out_index(this->get_node_out_index(j),jout))->get_direction_vector();
            //send to map to store branch angle
            double bangle = edge_out_dir.angle(edge_dir);
	    if(this->get_edge(j)->branch_count() == 1){
            if (bangle < M_PI)
            {
                this->set_branching_angle(j,this->get_edge_out_index(this->get_node_out_index(j),jout),edge_out_dir.angle(edge_dir));
	    }
        }
	    else
	      {
            //Default branching angle for trumpets is M_PI/4, although note that this may not actually be used in the
            //impaction term anyway.
		this->set_branching_angle(j,this->get_edge_out_index(this->get_node_out_index(j),jout),M_PI/4.0);
          }
	}
    }
}

double integrand_acin_grangle(const double & theta_par, const double & theta, const double & phi, const double & theta_min,
        const double & theta_max, const double & phi_min, const double & phi_max)
{
    //integrand for finding expected gravity angle for acinar edges
    return sin(acos(cos(theta) * cos(theta_par) - sin(theta) * cos(phi) * sin(theta_par))) / (M_PI * M_PI);
}

void AirwayTransportNetwork::compute_all_edge_gravity_angles(const pos::Position<double> grav_dir) {
    this->edge_vals[TN_EDGE_GRAVITY_ANGLE] = Eigen::VectorXd::Zero(this->count_edges());
    for(size_t j = 0; j < this->count_edges(); j++) {
        if (this->get_edge(j)->branch_count() == 1) {
            this->edge_vals[TN_EDGE_GRAVITY_ANGLE][j] = this->get_edge(j)->get_direction_vector().angle(grav_dir);
        } else {
            //Default gravity angle for trumpets is
            this->edge_vals[TN_EDGE_GRAVITY_ANGLE][j] = 0.0;
        }
    }
    //Assign gravity angles to weibel tree acinus units
        for (int i = 0; i < this->count_term_nodes(); ++i) {
            std::vector<size_t> node_vec = this->acin_node_map[i];
            size_t node_no = this->acin_node_map[i][0];
            //loop through number of nodes in the tree, excpet the last one
            for (int k = 0; k < node_vec.size()-1; ++k) {
                //for each node, starting with the top node, get the edge in and out
                size_t edge_in = this->get_edge_in_index(node_no,0);
                size_t edge_out = this->get_edge_out_index(node_no,0);
                //if the branch count is the same for edges in and out, set gravity angle to be equal to the gravity
                //angle of the edge in, since there is not branching event
                if (this->get_edge(edge_in)->branch_count() == this->get_edge(edge_out)->branch_count())
                {
                    this->edge_vals[TN_EDGE_GRAVITY_ANGLE][edge_out] = this->get_edge_gravity_angle(edge_in);
                }
                //if the branch count is different determine the new gravity angle based on the previous
                else
                {
                    //formula to determine new gravity angle after branching
                    //Calculates the expected value of sin(g_angle) based on the g_angle of the parent edge
                    /*
                    double g_angle_up = this->get_edge_gravity_angle(edge_in);//gravity angle of edge above

                    double theta_min = 0.0;
                    double theta_max = M_PI / 2.0;
                    double phi_min = 0.0;
                    double phi_max = 2 * M_PI;
                    double theta_p = g_angle_up;
                    auto inner_integral = [theta_p,theta_min,theta_max,phi_min,phi_max](double phi)
                            {
                            auto f = [theta_p, phi, theta_min, theta_max, phi_min, phi_max](double theta)
                                    {
                                        return integrand_acin_grangle(theta_p, theta, phi, theta_min, theta_max, phi_min, phi_max);
                                    };
                                return boost::math::quadrature::trapezoidal(f, theta_min, theta_max, 1e-3);
                            };
                    double result = boost::math::quadrature::trapezoidal(inner_integral, phi_min, phi_max, 1e-3);

                    //std::cout << g_angle_up << ' ' << result << std::endl;
                    this->edge_vals[TN_EDGE_GRAVITY_ANGLE][edge_out] = asin(result);
                    */

                    this->edge_vals[TN_EDGE_GRAVITY_ANGLE][edge_out] = asin(2.0/M_PI);
                    //this->edge_vals[TN_EDGE_GRAVITY_ANGLE][edge_out] = M_PI/3.0;
                }
                //move to next node down in the tree
                node_no = this->get_node_out_index(edge_out);
            }
        }



}

void AirwayTransportNetwork::print_end_transport_headers(const std::string &filename)
{
    std::ofstream output;
    output.open(filename);
    output << "AirwayNumber, ";
    output << "ParentAirway, ";
    output << "TotalDeposition_" << VolumeParam().get_phys_units();
    output << "Radius" << LengthParam().get_phys_units();
    output << "Length" << LengthParam().get_phys_units();
    output << "Nbranches";
    output << "WeibelGeneration";
    output << "HorsfieldOrder";
    output << "MaxFlux" << VolumeFluxParam().get_phys_units();
    output << "ImpactionDeposition_" << VolumeParam().get_phys_units() << '\n';
    output.close();
}

void AirwayTransportNetwork::print_end_transport_csv(const std::string & filename, SimParameterList *p)
{
    std::ofstream output;
    output.open(filename);
	Eigen::VectorXd deps = this->get_all_airway_depositions();
    Eigen::VectorXd deps_imp = this->get_all_airway_imp_depositions();
    Eigen::VectorXd deps_sed = this->get_all_airway_sed_depositions();
    Eigen::VectorXd deps_dif = this->get_all_airway_diff_depositions();
	Eigen::VectorXd max_fluxes = this->get_all_edge_max_fluxes();
	Eigen::VectorXd lengths = this->get_all_airway_lengths();
    for (int J = 0; J < this->count_airways(); ++J) {
        //if (this->get_edge(this->get_edges_in_airway(J,0))->branch_count() == 1)
        //{
            //Print airway number
            output << J << ", ";
            //Print parent airway number (-1 for trachea)
            for (int j = 0; j < this->count_airway_edges(J); ++j) {
                size_t kup = this->get_node_in_index(this->get_edges_in_airway(J,j));
                if (kup == 0)
                {
                    output << -1 << ", ";
                    break;
                }
                else {
					size_t jup = this->get_edge_in_index(kup, 0);
                    size_t Jup = this->get_airway_from_edge(jup);
                    if (Jup != J) {
                        output << Jup << ", ";
                        break;
                    }
                }
            }
            double dep = deps(this->get_edges_in_airway(J,0)) / p->get_conversion(VOL_KEY);
            double dep_imp = deps_imp(this->get_edges_in_airway(J,0)) / p->get_conversion(VOL_KEY);
            //Print deposition in airway
            output << dep << ", ";
            //Print radius
            output << this->get_edge(this->get_edges_in_airway(J,0))->get_geom()->get_inner_radius() / p->get_conversion(LENGTH_KEY) << ", ";
            //Print length
            double airway_length = 0;
                for (int i = 0; i < this->count_airway_edges(J); ++i) {
                    airway_length += this->get_edge(this->get_edges_in_airway(J,i))->get_geom()->get_length();
                }
                //std::cout << lengths[this->get_edges_in_airway(J,0)] / p->get_conversion(LENGTH_KEY) << ", ";
            output << airway_length / p->get_conversion(LENGTH_KEY) << ", ";
            //output << lengths[this->get_edges_in_airway(J,0)] / p->get_conversion(LENGTH_KEY) << ", ";
            //Print Nbranches
            output << this->get_edge(this->get_edges_in_airway(J,0))->branch_count() << ", ";
            //Print Weibel generation
            output << this->get_weibel_order(this->get_edges_in_airway(J,0)) << ", ";
            //Print Horsfield order
            output << this->get_horsfield_order(this->get_edges_in_airway(J,0)) << ", ";
            //Print max flux
            output << max_fluxes(J) / p->get_conversion(VOL_FLUX_KEY) << ", ";
            //Impaction Deposition
            output << dep_imp << '\n';
        //}
    }
    output.close();
}

//airway no, parent airway, deposition info, ...
//do it same as printing csv

int AirwayTransportNetwork::print_transport_vtk(const std::string & filename, SimParameterList *p)
{
	using namespace std;
	std::unordered_map<string, vector<double>> extra_vals;

	for(size_t n = 0; n < TN_NODE_COUNT; n++)
	{
		extra_vals[tn_node_strings[n]] = vector<double>();
		extra_vals[tn_node_strings[n]].resize(this->count_nodes());
	}
    for(size_t k = 0; k < this->count_nodes(); k++)
    {
        extra_vals[tn_node_strings[TN_NODE_CONC]][k] = this->node_vals[TN_NODE_CONC][k];
        extra_vals[tn_node_strings[TN_NODE_DEPOSITION]][k] = this->node_vals[TN_NODE_DEPOSITION][k] * p->get_conversion(TIME_KEY) / p->get_conversion(VOL_KEY);
        extra_vals[tn_node_strings[TN_NODE_VOLUME]][k] = this->node_vals[TN_NODE_VOLUME][k] / p->get_conversion(VOL_KEY);
        extra_vals[tn_node_strings[TN_NODE_CUMULATIVE_DEPOSITION]][k] = this->node_vals[TN_NODE_CUMULATIVE_DEPOSITION][k] / p->get_conversion(VOL_KEY);
        extra_vals[tn_node_strings[TN_NODE_CUMULATIVE_DIFFUSIVE_DEPOSITION]][k] = this->node_vals[TN_NODE_CUMULATIVE_DIFFUSIVE_DEPOSITION][k] / p->get_conversion(VOL_KEY);
        extra_vals[tn_node_strings[TN_NODE_CUMULATIVE_SEDIMENTATION_DEPOSITION]][k] = this->node_vals[TN_NODE_CUMULATIVE_SEDIMENTATION_DEPOSITION][k] / p->get_conversion(VOL_KEY);
        extra_vals[tn_node_strings[TN_NODE_CUMULATIVE_IMPACTION_DEPOSITION]][k] = this->node_vals[TN_NODE_CUMULATIVE_IMPACTION_DEPOSITION][k] / p->get_conversion(VOL_KEY);
		extra_vals[tn_node_strings[TN_NODE_CUMULATIVE_ALVEOLAR_DEPOSITION]][k] = this->node_vals[TN_NODE_CUMULATIVE_ALVEOLAR_DEPOSITION][k] / p->get_conversion(VOL_KEY);
        extra_vals[tn_node_strings[TN_NODE_TOTAL_DEP_OVER_VOL]][k] = this->node_vals[TN_NODE_TOTAL_DEP_OVER_VOL][k];
        extra_vals[tn_node_strings[TN_NODE_BRANGLE]][k] = this->node_vals[TN_NODE_BRANGLE][k];
        //extra_vals[tn_node_strings[TN_NODE_TOTAL_DEPOSITION]][k] = this->node_vals[TN_NODE_TOTAL_DEPOSITION][k] / p->get_conversion(VOL_KEY);
        extra_vals[tn_node_strings[TN_NODE_GENERATIONAL_DEPOSITION]][k] = this->node_vals[TN_NODE_GENERATIONAL_DEPOSITION][k] / p->get_conversion(VOL_KEY);
        extra_vals[tn_node_strings[TN_NODE_PENDELLUFT]][k] = this->node_vals[TN_NODE_PENDELLUFT][k] / p->get_conversion(VOL_KEY);
        extra_vals[tn_node_strings[TN_NODE_SPECT_DEPOSITION]][k] = this->node_vals[TN_NODE_SPECT_DEPOSITION][k] / p->get_conversion(VOL_KEY);
    }
	for(size_t n = 0; n < TN_EDGE_COUNT; n++)
	{
		extra_vals[tn_edge_strings[n]] = vector<double>();
		extra_vals[tn_edge_strings[n]].resize(this->count_edges());
	}
	for(size_t j = 0; j < this->count_edges(); j++) 
	{
		extra_vals[tn_edge_strings[TN_EDGE_FLUX]][j] = this->edge_vals[TN_EDGE_FLUX][j] * p->get_conversion(TIME_KEY) / p->get_conversion(VOL_KEY);
		extra_vals[tn_edge_strings[TN_EDGE_DIFFUSIVITY]][j] = this->edge_vals[TN_EDGE_DIFFUSIVITY][j] * p->get_conversion(TIME_KEY) / p->get_conversion(VOL_KEY);
        extra_vals[tn_edge_strings[TN_EDGE_GRAVITY_ANGLE]][j] = this->edge_vals[TN_EDGE_GRAVITY_ANGLE][j];
        extra_vals[tn_edge_strings[TN_EDGE_AIRWAY_DEPOSITION]][j] = this->edge_vals[TN_EDGE_AIRWAY_DEPOSITION][j] / p->get_conversion(VOL_KEY);
        extra_vals[tn_edge_strings[TN_EDGE_AIRWAY_IMP_DEPOSITION]][j] = this->edge_vals[TN_EDGE_AIRWAY_IMP_DEPOSITION][j] / p->get_conversion(VOL_KEY);
        extra_vals[tn_edge_strings[TN_EDGE_AIRWAY_DIFF_DEPOSITION]][j] = this->edge_vals[TN_EDGE_AIRWAY_DIFF_DEPOSITION][j] / p->get_conversion(VOL_KEY);
        extra_vals[tn_edge_strings[TN_EDGE_AIRWAY_SED_DEPOSITION]][j] = this->edge_vals[TN_EDGE_AIRWAY_SED_DEPOSITION][j] / p->get_conversion(VOL_KEY);
        extra_vals[tn_edge_strings[TN_EDGE_AIRWAY_NO]][j] = this->edge_vals[TN_EDGE_AIRWAY_NO][j];
        extra_vals[tn_edge_strings[TN_EDGE_GENERATION_NO]][j] = this->edge_vals[TN_EDGE_GENERATION_NO][j];
        extra_vals[tn_edge_strings[TN_EDGE_MAX_FLUX]][j] = this->edge_vals[TN_EDGE_MAX_FLUX][j] * p->get_conversion(TIME_KEY) / p->get_conversion(VOL_KEY);
        extra_vals[tn_edge_strings[TN_EDGE_AIRWAY_LENGTH]][j] = this->edge_vals[TN_EDGE_AIRWAY_LENGTH][j] / p->get_conversion(LENGTH_KEY);
        extra_vals[tn_edge_strings[TN_EDGE_HORSFIELD_ORDER]][j] = this->edge_vals[TN_EDGE_HORSFIELD_ORDER][j];
    }
	
	return (this->print_vtk(filename, 1.0 / p->get_conversion(LENGTH_KEY), extra_vals, !print_acinar_airways));
}

void AirwayTransportNetwork::update_conducting_diffusivities()
{
	for(size_t i_e = 0; i_e < this->conducting_edges.size(); i_e++)
	{
		size_t j = this->conducting_edges[i_e];

		this->edge_vals[TN_EDGE_DIFFUSIVITY][j] = this->get_edge(j)->branch_count() * this->get_edge(j)->get_geom()->inner_area()
			                                      * this->conducting_diffusivity->calc_diffusivity(this->get_edge(j)->get_geom(), 
												  this->edge_vals[TN_EDGE_FLUX][j] / this->get_edge(j)->branch_count()) /
												  this->get_edge(j)->get_geom()->get_length();


    }
}

void AirwayTransportNetwork::update_acinar_diffusivities()
{
	for(size_t kt = 0; kt < this->acin_edge_map.size(); kt++)
	{
		for(size_t ja = 0; ja < this->acin_edge_map[kt].size(); ja++)
		{
			size_t j = this->acin_edge_map[kt][ja];
			this->edge_vals[TN_EDGE_DIFFUSIVITY][j] = this->get_edge(j)->branch_count() * this->get_edge(j)->get_geom()->inner_area()
			                                      * this->acinar_diffusivity->calc_diffusivity(this->get_edge(j)->get_geom(), 
												  this->edge_vals[TN_EDGE_FLUX][j] / this->get_edge(j)->branch_count()) /
												  this->get_edge(j)->get_geom()->get_length();
		}
	}
}