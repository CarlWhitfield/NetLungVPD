#include <string>
#include <iostream>
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    #define IS_WINDOWS 1
#else
    #define IS_WINDOWS 0
#endif
#if IS_WINDOWS
    #include <direct.h>
#else
    #include <unistd.h>
#endif
#include "lung_simulation.h"
#include "list_template.h"
#include "define_codes.h"
#include "define_params.h"
//#include<vtkSmartPointer.h>
//#include<vtkSphereSource.h>
//#include<vtkPolyData.h>
//#include<vtkFieldData.h>
//#include<vtkDoubleArray.h>
//#include<vtkMultiBlockDataSet.h>
//#include<vtkXMLMultiBlockDataWriter.h>

using namespace inlist;

LungSimulation::LungSimulation(const int &argc, char** argv)  //constructs lung simulation from command line args
{
	options = std::make_shared<SimOptionList>();         //start by assigning default option and paramter lists
	params = std::make_shared<SimParameterList>();

	//5 possible arguments in total, all filenames, sort by extensions
	std::vector<std::string> extensions;
	extensions.resize(7);
	extensions[0] = OPTIONS_FILE_EXT;
	extensions[1] = PARAMS_FILE_EXT;
	extensions[2] = NODE_FILE_EXT;
	extensions[3] = TERM_NODE_FILE_EXT;
	extensions[4] = BRANCH_FILE_EXT;
	extensions[5] = FLOW_INPUT_FILE_EXT;
	extensions[6] = TNODE_MAP_FILE_EXT;

	this->options->get_filenames_from_args(extensions, argc, argv);
	this->sort_files();
	if( this->options->filename_exists(OPTIONS_FILE_EXT) ) this->options->read_file(options->get_filename(OPTIONS_FILE_EXT));  //if options file given, read it in
	if( this->options->filename_exists(PARAMS_FILE_EXT) ) this->params->read_file(options->get_filename(PARAMS_FILE_EXT));  //if options file given, read it in
	this->params->check_and_convert(options.get());  //run checks and conversions here
	if( this->options->get_option<bool>(ACINUS_FROM_FILE)->get_value() ) this->read_acin_template();
}

int LungSimulation::initialise()   //inititalise networks
{
	std::cout << "Initialising...\n";

	this->create_flow_network();//Loads airway network from file (assuming network is from file)
	//flow_tree = new AirwayFlowNetwork(options, params);
#if IS_WINDOWS
	if(this->options->get_option<bool>(PRINT_ACINUS_VTK_KEY)->get_value()) _mkdir("acini");
	if (this->options->get_option<bool>(PRINT_FLOW_VTK_KEY)->get_value()) _mkdir("flow_tree");
#else
    //Make directories for printing vtks
    if(this->options->get_option<bool>(PRINT_ACINUS_VTK_KEY)->get_value()) mkdir("acini",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(this->options->get_option<bool>(PRINT_FLOW_VTK_KEY)->get_value()) mkdir("flow_tree",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
	std::cout << "Flow Tree Initialised.\n";
    //Build transport network
	if( this->options->get_option<bool>(ACINUS_FROM_FILE)->get_value() )
	{
		this->transport_tree = std::make_shared<AirwayTransportNetwork>(this->flow_tree.get(), 
			             this->flow_tree->get_all_termnode_intrinsic_volumes(), this->options.get(), this->params.get(), 
			             this->acin_template.get());
	}
	else
	{
		this->transport_tree = std::make_shared<AirwayTransportNetwork>(this->flow_tree.get(), 
			this->flow_tree->get_all_termnode_intrinsic_volumes(), this->options.get(), this->params.get());
	}
	std::cout << "Transport tree built.\n";
#if IS_WINDOWS
	if (this->options->get_option<bool>(PRINT_TRANSPORT_VTK_KEY)->get_value()) _mkdir("transport_tree");
#else
    //Make directory for printing transport vtk
	if(this->options->get_option<bool>(PRINT_TRANSPORT_VTK_KEY)->get_value()) mkdir("transport_tree", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
	//this->create_tree_maps();
	this->copy_fluxes_to_transport_network();
	return 0;
}

int LungSimulation::read_acin_template()
{
	//read in acinus 
	this->acin_template = std::make_shared<AirwayNetwork>(this->acin_files[NODE_FILE_EXT], 
		           this->acin_files[BRANCH_FILE_EXT], this->acin_files[TERM_NODE_FILE_EXT], 1.0);

	//check for extra args
	std::vector<char> extra_node_args = this->acin_template->get_extra_node_args();
	std::vector<char> extra_edge_args = this->acin_template->get_extra_edge_args();
	//deal with extra args
	for(size_t i_arg = 0; i_arg < extra_node_args.size(); i_arg++)
	{
		for(size_t k = 0; k < this->acin_template->count_nodes(); k++)
		{
			std::vector<double> input_vals = this->acin_template->get_extra_node_inputs(extra_node_args[i_arg], this->acin_template->get_node(k));
			switch(extra_node_args[i_arg])
			{
			case 'p':   //node position has changed
				{
					//do not really need to do anything
					if(input_vals.size() > 2)
					{
						this->acin_template->get_node(k)->set_original_pos(
							   network::Position(input_vals[0], input_vals[1], input_vals[2]));
					}
				} break;

			case 'r':   //node seed rad
				{
					//do not really need to do anything
					if(input_vals.size() > 0)
					{
						this->acin_template->get_node(k)->set_seed_rad(input_vals[0]);
					}
				} break;

			default:
				{
					std::cout << "Error did not recognise node option: " << extra_node_args[i_arg] << '\n';
				} break;
			}
		}
	}

	for(size_t i_arg = 0; i_arg < extra_edge_args.size(); i_arg++)
	{
		for(size_t j = 0; j < this->acin_template->count_edges(); j++)
		{
			std::vector<double> input_vals = this->acin_template->get_extra_edge_inputs(extra_edge_args[i_arg], this->acin_template->get_edge(j));
			switch(extra_edge_args[i_arg])
			{
			case 'o':    //outer radius
				{
					if(input_vals.size() > 0)
					{
						this->acin_template->get_edge(j)->update_outer_radius(input_vals[0]);
					}
				} break;

			default:
				{
					std::cout << "Error did not recognise edge option: " << extra_edge_args[i_arg] << '\n';
				}

			}
		}
	}


	return 0;
}

int LungSimulation::create_flow_network()
{
	//whether network is to be read from file or built from default
	if (options->get_option<char>(TREE_KEY)->get_value() == MODEL_FROM_FILE)
	{
		//check all files have been given and they exist
		bool filesOK = true;
		std::string types[] = {NODE_FILE_EXT, TERM_NODE_FILE_EXT, BRANCH_FILE_EXT};
		for(int i = 0; i < 3; i++)
		{
			if(!file_exists(this->tree_files[types[i]]))
			{
				//either file was not inputted or the path doesn't exist
				//revert to default
				options->get_option<char>(TREE_KEY)->set_to_default();
				std::cout << "Error, cannot find " << types[i] << " file, cannot build airway network. Using default tree model (" 
						<< options->get_option<char>(TREE_KEY)->print() << ") instead.\n";
				filesOK = false;
				break;   //stop trying to process these
			}
		}

		//continue if no problems exist
		if(filesOK) 
		{
			if(this->import_network())
			{
				options->get_option<char>(TREE_KEY)->set_to_default();
				std::cout << "Error, problem importing network. Using default tree model (" 
					 << options->get_option<char>(TREE_KEY)->print() << ") instead.\n";
			}
		}
		else
		{
			std::cout << "Problem with input files, aborting...\n";
			return 1;
		}

	}
	else
	{
		if (options->get_option<char>(TREE_KEY)->get_value() == MODEL_M)
		{
		//build_model_m(); 
		}
		else
		{
			std::cout << "Invalid tree option, aborting...\n";
			return 1;
		}
	}

	return 0;
}

int LungSimulation::import_network()  //import airway network from file
{
	LengthParam length_scale(1.0, std::string("scale"), params->get_conversion_ptr(LENGTH_KEY));
	std::cout << "Reading...\n";
	this->flow_tree = std::make_shared<AirwayFlowNetwork>(this->tree_files[NODE_FILE_EXT], this->tree_files[BRANCH_FILE_EXT],
							 this->tree_files[TERM_NODE_FILE_EXT], length_scale.get_value());
	this->flow_tree->setup_afn(options.get(), params.get());

	return 0;
}

int LungSimulation::simulate()  //run simulation
{
    //Load time step
	double dt = params->get_param<double>(TIME_STEP_KEY)->get_value();
    //Time step counter
	size_t itime = 0;
	int nf = 0;
	std::string washout_filename("washout.csv");
	std::string acinus_filename("acinus.csv");
	if(this->options->get_option<bool>(PRINT_ACINUS_CSV_KEY)->get_value()) this->print_acinus_headers(acinus_filename);
	this->print_washout_file_headers(washout_filename);
    //Counter for printing
	double next_print = 0.0;
    //Start timer for simulation
	auto start = std::chrono::system_clock::now();
	double Nb = params->get_param<double>(BREATH_COUNT_KEY)->get_value();
	int Nout = (Nb / (100*dt));
	if (options->get_option<char>(INHALATION_MODE_KEY)->get_value() == BOLUS_CONC){
		this->transport_tree->set_source_mode(false);  //switch bolus source off initially
	}
	else{
		this->transport_tree->set_source_mode(true);  //switch source on initially if not bolus
	}

	while (itime*dt < Nb + 1E-09)
	{
		this->time = itime*dt;
		if(this->transport_tree->get_washin_mode() && time >= params->get_param<double>(WASHIN_END_KEY)->get_value())
		{
			this->transport_tree->set_washin_mode(false);  //switch to washout
		}
		//switch bolus source on or off
		if (this->transport_tree->get_washin_mode() &&  options->get_option<char>(INHALATION_MODE_KEY)->get_value() == BOLUS_CONC) //if bolus conc
		{
			double DV = this->flow_tree->get_lung_expansion_volume(itime);
			double VB0 = params->get_param<double>(BOLUS_START_INSP_KEY)->get_value();
			double DVB = params->get_param<double>(BOLUS_VOL_DURATION_KEY)->get_value();
			if ((this->transport_tree->get_source_mode() == false) && (DV >= VB0) && (DV <= VB0 + DVB))
			{
				this->transport_tree->set_source_mode(true);
			}
			if ((this->transport_tree->get_source_mode() == true) && ((DV < VB0) || (DV > VB0 + DVB)))
			{
				this->transport_tree->set_source_mode(false);
			}
		}



        //Print vtk flow file if at the specified time to print
		if (time >= next_print - 1E-09)
		{
			std::cout << "Printing outputs..." << std::endl;
			std::stringstream ss;
			if (this->options->get_option<bool>(PRINT_FLOW_VTK_KEY)->get_value()) {
				ss << "flow_tree/flow_tree_" << nf;
				this->flow_tree->print_flow_vtk(ss.str(), this->params.get());
				ss.clear();
				ss.str("");
			}
			if(options->get_option<bool>(PRINT_ACINUS_VTK_KEY)->get_value())
			{
				ss << "acini/acini_" << nf;
				this->print_acinus_vtk(ss.str(), nf);
				ss.clear();
				ss.str("");
			}
			if (this->options->get_option<bool>(PRINT_TRANSPORT_VTK_KEY)->get_value())
			{
				ss << "transport_tree/transport_tree_" << nf;
				this->transport_tree->print_transport_vtk(ss.str(), this->params.get());
				ss.clear();
				ss.str("");
			}
			if (this->options->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value()) {
				ss << "deposited_" << nf << ".csv";
				std::string deposited_filename(ss.str());
				this->transport_tree->print_end_transport_csv(deposited_filename, this->params.get());
			}
			nf++;
			next_print += this->params->get_param<double>(PRINT_TIME_KEY)->get_value();
			std::cout << "Outputs printed." << std::endl;
		}
		if(this->options->get_option<bool>(PRINT_ACINUS_CSV_KEY)->get_value()) this->append_acinus_vals(acinus_filename, time);
		//auto pstart = std::chrono::system_clock::now();
        this->append_washout_file(washout_filename, time);
		//auto pend = std::chrono::system_clock::now();
		//std::cout << "Print took " << (std::chrono::duration<double>(pend-pstart)).count() << '\n';
		/*std::cout << this->flow_tree->get_termnode_volume(0)/this->flow_tree->get_node(this->flow_tree->get_first_term_index())->point_count() << ", ";
		std::cout << this->flow_tree->get_node_pressure(this->flow_tree->get_first_term_index()) << ", ";
		std::cout << this->flow_tree->get_termnode_volume(6)/this->flow_tree->get_node(this->flow_tree->get_first_term_index()+6)->point_count() << ", ";
		std::cout << this->flow_tree->get_node_pressure(this->flow_tree->get_first_term_index() + 6) << "\n";*/
		if((itime+1)*dt < params->get_param<double>(BREATH_COUNT_KEY)->get_value() + 1E-09)
		{
	//		auto fsstart = chrono::system_clock::now();
            //Time-step the ventilation solution
			this->flow_tree->solve_pressure_volume_update(itime+1);
	//		auto fsend = chrono::system_clock::now();
	//		std::cout << "Flow update took " << (chrono::duration<double>(fsend-fsstart)).count() << '\n';
		
	//		auto cfstart = chrono::system_clock::now();
            //Convert flux solution to transport network
			this->copy_fluxes_to_transport_network();
	//		auto cfend = chrono::system_clock::now();
	//		std::cout << "Flow copy took " << (chrono::duration<double>(cfend-cfstart)).count() << '\n';
		
	//		auto tsstart = chrono::system_clock::now();
            //Time-step the transport solution
			this->transport_tree->solve_concentration_update(dt);
	//		auto tsend = chrono::system_clock::now();
	//		std::cout << "Transport solve took " << (chrono::duration<double>(tsend-tsstart)).count() << '\n';
		}
		if ((itime+1) % Nout == 0) std::cout << "Transport simulation " << 100 * (this->time + dt) / Nb << "% complete." << std::endl;
        //Move to next time step
		itime++;
	}
	auto end = std::chrono::system_clock::now();
	std::cout << "Simulation took " << (std::chrono::duration<double>(end-start)).count() << '\n';

	return 0;
}

//sort filenames into acin and tree files
void LungSimulation::sort_files()
{
	std::string types[] = {NODE_FILE_EXT, TERM_NODE_FILE_EXT, BRANCH_FILE_EXT};
	for(int i = 0; i < 3; i++)
	{
		if(options->filename_exists(types[i])) 
		{
			size_t file_count = this->options->count_files_with_ext(types[i]);
			//if file name contains "acin" it must be acinus file
			for(size_t n=0; n<file_count; n++)
			{
				if(this->options->get_filename(types[i],n).find(std::string("acin")) != std::string::npos)  //if string contains "acin"
				{
					this->acin_files[types[i]] = this->options->get_filename(types[i],n);
					std::cout << "acin " << types[i] << " filename: " << this->options->get_filename(types[i],n) << '\n';
				}
				else 
				{
					this->tree_files[types[i]] = this->options->get_filename(types[i],n);
					std::cout << types[i] << " filename: " << this->options->get_filename(types[i],n) << '\n';
				}
			}
		}		
	}
	//TODO: Simplify by removing acin files and take it as .acin_nodes
}

////the map should be stored in lung simulation object and be directed from flow -> transport
//void LungSimulation::create_tree_maps()
//{
//	//create map from flow network to transport network
//	this->transport_edge_map.clear();
//	this->transport_edge_map.resize(this->flow_tree->count_edges());
//	for(size_t j = 0; j < this->flow_tree->count_edges(); j++)    //loop over flow tree edges
//	{
//		size_t kfin = this->flow_tree->get_node_in_index(j);        //node into flow edge
//		size_t kfout = this->flow_tree->get_node_out_index(j);       //node out of flow edge
//		size_t ktin = this->transport_tree->get_node_index(this->flow_tree->get_node(kfin));     //index of kfin on transport tree
//		size_t ktout = this->transport_tree->get_node_index(this->flow_tree->get_node(kfout));   //index of kfout on transport tree
//
//		//work way back up from ktout to ktin
//		this->transport_edge_map[j].push_back(this->transport_tree->get_edge_in_index(ktout, 0));   //edge feeding ktout maps to j
//		size_t kh = this->transport_tree->get_node_in_index(this->transport_tree->get_edge_in_index(ktout, 0));    //index of next node
//		std::vector<size_t> intermediate_nodes;
//		while(kh != ktin)  //only works so long as each node has max 1 edge in
//		{
//			intermediate_nodes.push_back(kh);
//			this->transport_edge_map[j].push_back(this->transport_tree->get_edge_in_index(kh, 0));     //intermediate edges map to j
//			kh = this->transport_tree->get_node_in_index(this->transport_tree->get_edge_in_index(kh, 0));  //index of next node
//		}
//	}
//}

//and then fix these
void LungSimulation::copy_fluxes_to_transport_network() //this should only happen for conducting network -> not for acinus
{
	for(size_t kt = 0; kt < this->flow_tree->count_term_nodes(); kt++)
	{
		this->transport_tree->update_acinar_volume(kt, this->flow_tree->get_termnode_volume(kt));
	}
	this->transport_tree->compute_all_fluxes(this->params->get_param<double>(TIME_STEP_KEY)->get_value());

	//diffusivity depends on updated flow problem if diffusivity type is not molecular everywhere
	if(this->options->get_option<char>(DIFFUSIVITY_TYPE_KEY)->get_value() != MOLECULAR_DIFFUSIVITY)
	{
		this->transport_tree->update_conducting_diffusivities();
        this->transport_tree->update_acinar_diffusivities();
	}
}

void LungSimulation::print_acinus_headers(const std::string &filename)
{
	std::ofstream output;
	output.open(filename);
	output << "Time_s";
	for (size_t kt = 0; kt < this->transport_tree->count_acini(); kt++)
	{
		output << ", " << "Volume_L_" << kt;
	}
	for (size_t kt = 0; kt < this->transport_tree->count_acini(); kt++)
	{
		output << ", " << "IG_Volume_L_" << kt;
	}
	output << '\n';
	output.close();

	output.open("acinus_positions.csv");
	output << "Acinus number, x(mm), y(mm), z(mm)";
	double L_scale = this->params->get_conversion(LENGTH_KEY);
	for (size_t kt = 0; kt < this->transport_tree->count_acini(); kt++)
	{
		network::Position pos;
		this->transport_tree->get_acinus_pos(kt, pos);
		output << kt << ", " << pos.x[0]*L_scale <<  ", " << pos.x[1]*L_scale 
			   <<  ", " << pos.x[2]*L_scale << '\n';
	}
	output.close();
}

//outputs should be done differently 
void LungSimulation::append_acinus_vals(const std::string &filename, const double & t)
{
	std::ofstream output;
	output.open(filename, std::ios_base::app);
	double vol_conv = this->params->get_conversion(VOL_KEY);
	double time_conv = this->params->get_conversion(TIME_KEY);
	output << std::fixed << std::setprecision(OUTPUT_PRECISION_DEFAULT);
	output << t / time_conv;
	for (size_t kt = 0; kt < this->transport_tree->count_acini(); kt++)
	{
		output << ", " << this->transport_tree->get_acinus_volume(kt) / vol_conv;
	}
	for (size_t kt = 0; kt < this->transport_tree->count_acini(); kt++)
	{
		output << ", " << this->transport_tree->get_acinus_IG_volume(kt, 
								this->options->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value()) / vol_conv;
	}
	output << '\n';
	output.close();
}

void LungSimulation::print_washout_file_headers(const std::string &filename)
{
	std::ofstream output;
	output.open(filename);
	output << "Time_" << TimeParam().get_phys_units() << ", ";
	output << "MouthFlux_" << VolumeParam().get_phys_units() << '/' << TimeParam().get_phys_units() << ", ";
	output << "MouthGasFlux_" << VolumeParam().get_phys_units() << '/' << TimeParam().get_phys_units() << ", ";
	output << "LungVolume_" << VolumeParam().get_phys_units() << ", ";
	output << "LungIGVolume_" << VolumeParam().get_phys_units() << ", ";
	output << "MouthConc, ";
    output << "PleuralPressure_" << PressureParam().get_phys_units() << ", ";
    output << "Total_deposition_" << VolumeParam().get_phys_units() << ", ";
    output << "Total_diffusion_deposition_" << VolumeParam().get_phys_units() << ", ";
    output << "Total_sedimentation_deposition_" << VolumeParam().get_phys_units() << ", ";
    output << "Total_impaction_deposition_" << VolumeParam().get_phys_units() << ", ";
	output << "Total_alveolar_deposition_" << VolumeParam().get_phys_units() << ", ";
    output << "Acinar_deposition_" << VolumeParam().get_phys_units() << ", ";
    output << "Acinar_diffusion_deposition_" << VolumeParam().get_phys_units() << ", ";
    output << "Acinar_sedimentation_deposition_" << VolumeParam().get_phys_units() << ", ";
    output << "Acinar_impaction_deposition_" << VolumeParam().get_phys_units() << ", ";
	output << "Acinar_alveolar_deposition_" << VolumeParam().get_phys_units() << ", ";
    /*
    int N_gen = 1;//no. of airway generations in tree
    for (int k = 1; k < this->transport_tree->count_nodes(); ++k) {
        size_t j = this->transport_tree->get_edge_in_index(k,0);
        if (N_gen < this->transport_tree->get_weibel_order(j))
        {
            N_gen = this->transport_tree->get_weibel_order(j);
        }
    }
    for (int i = 0; i < N_gen; ++i) {
        if (i < N_gen - 1)
        {
            output << "Deposition_generation_" << i << "_" << VolumeParam().get_phys_units() << ", ";
        }
        else
        {
            output << "Deposition_generation_" << i << "_" << VolumeParam().get_phys_units() << ", ";
        }
    }
    for (int i = 0; i < N_gen; ++i) {
        if (i < N_gen - 1)
        {
            output << "Deposition_CV_generation_" << i << ", ";
        }
        else
        {
            output << "Deposition_CV_generation_" << i << ", ";
        }
    }

    //Flux and deposition for all 6th genertaion airways
    for (int J = 0; J < this->transport_tree->count_airways(); ++J) {
        if (this->transport_tree->get_weibel_order(this->transport_tree->get_edges_in_airway(J,0)) == 6)
        {
            output << "Flux_airway_" << J << ", ";
            output << "Deposition_airway_" << J << ", ";
        }
    }
*/

    output << "Zero" << '\n';

	output.close();
}

void LungSimulation::append_washout_file(const std::string &filename, const double & t)
{
	std::ofstream output;
	output.open(filename, std::ios_base::app);
	double vol_conv = 1.0/ this->params->get_conversion(VOL_KEY);
	double time_conv = 1.0 / this->params->get_conversion(TIME_KEY);
	double press_conv = 1.0 / this->params->get_conversion(PRESSURE_KEY);
	output << std::fixed << std::setprecision(OUTPUT_PRECISION_DEFAULT);
	output << t * time_conv << ", "
	       << this->transport_tree->get_edge_flux(0) * vol_conv / time_conv << ", ";
	if( this->transport_tree->get_mouth_object() != NULL)
	{
		output << (this->transport_tree->get_edge_flux(0) * this->transport_tree->get_mouth_object()->get_conc_top()) * vol_conv / time_conv << ", ";
	}
	else
	{
		output << (std::max(this->transport_tree->get_edge_flux(0), 0.) * this->transport_tree->get_node_conc(this->transport_tree->get_node_in_index(0))
			     + std::min(this->transport_tree->get_edge_flux(0), 0.) * this->transport_tree->get_node_conc(this->transport_tree->get_node_out_index(0)))
			     * vol_conv / time_conv << ", ";
	}
	output << (this->flow_tree->get_total_edge_volume() + this->flow_tree->sum_all_termnode_volumes()) *vol_conv << ", "
	       << this->transport_tree->get_tot_IG_volume(this->options->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value())*vol_conv << ", ";
	if(this->transport_tree->get_mouth_object() != NULL)
	{
		output << this->transport_tree->get_mouth_object()->get_conc_top() << ", ";
	}
	else
	{
		output << this->transport_tree->get_node_conc(0) << ", ";
	}
	output << this->flow_tree->get_Ppl()*press_conv << ", ";
    //Total depositions
    double total_dep = 0;
    double total_diff_dep = 0;
    double total_sed_dep = 0;
    double total_imp_dep = 0;
	double total_alv_dep = 0;
    for (int k = 0; k < this->transport_tree->count_nodes(); ++k) {
        total_dep += this->transport_tree->get_total_deposition(k);
        total_diff_dep += this->transport_tree->get_total_diff_deposition(k);
        total_sed_dep += this->transport_tree->get_total_sed_deposition(k);
        total_imp_dep += this->transport_tree->get_total_imp_deposition(k);
		total_alv_dep += this->transport_tree->get_total_alv_deposition(k);
    }
    output << total_dep * vol_conv << ", " << total_diff_dep * vol_conv << ", " << total_sed_dep * vol_conv << ", " 
		   << total_imp_dep * vol_conv << ", " << total_alv_dep * vol_conv << ", ";

    //Acinar depositions
    double acin_cumul_dep = 0;
    double acin_cumul_diff_dep = 0;
    double acin_cumul_sed_dep = 0;
    double acin_cumul_imp_dep = 0;
	double acin_cumul_alv_dep = 0;
    for (int kt = 0; kt < this->transport_tree->count_term_nodes(); ++kt) {
        std::vector<size_t> acin_nodes = this->transport_tree->get_acin_node_map(kt);
        for (int i = 0; i < acin_nodes.size(); ++i) {
            size_t ka = acin_nodes[i];
            acin_cumul_dep += this->transport_tree->get_total_deposition(ka);
            acin_cumul_diff_dep += this->transport_tree->get_total_diff_deposition(ka);
            acin_cumul_sed_dep += this->transport_tree->get_total_sed_deposition(ka);
            acin_cumul_imp_dep += this->transport_tree->get_total_imp_deposition(ka);
			acin_cumul_alv_dep += this->transport_tree->get_total_alv_deposition(ka);
        }
    }
    output << acin_cumul_dep * vol_conv << ", " << acin_cumul_diff_dep * vol_conv << ", " << acin_cumul_sed_dep * vol_conv 
		   << ", " << acin_cumul_imp_dep * vol_conv << ", " << acin_cumul_alv_dep * vol_conv << ", ";


    /*
    //Generational depositions
    int N_gen = 1;//no. of airway generations in tree
    for (int k = 1; k < this->transport_tree->count_nodes(); ++k) {
        size_t j = this->transport_tree->get_edge_in_index(k,0);
        if (N_gen < this->transport_tree->get_weibel_order(j) + 1)
        {
            N_gen = this->transport_tree->get_weibel_order(j) + 1;
        }
    }
    Eigen::VectorXd gen_dep = Eigen::VectorXd::Zero(N_gen);
    for (int k = 0; k < this->transport_tree->count_nodes(); ++k) {
        size_t j;
        if (k == 0) {j = this->transport_tree->get_edge_out_index(k,0);}
        else {j = this->transport_tree->get_edge_in_index(k,0);}
        size_t gen = this->transport_tree->get_weibel_order(j);
        if (k >= this->transport_tree->count_nodes())
        {
            std::cout << "Problem here. k" << k << ' ' << this->transport_tree->count_nodes() << std::endl;
        }
        if (gen >= N_gen)
        {
            std::cout << "Problem here. gen" << gen << ' ' << N_gen << std::endl;
        }
        gen_dep[gen] += this->transport_tree->get_total_deposition(k);
    }
    for (int i = 0; i < N_gen; ++i) {
        if(i < N_gen - 1) {
            output << gen_dep[i] * vol_conv << ", ";
        }
        else
        {
            output << gen_dep[i] * vol_conv << ", ";
        }
    }

    //Generational deposition standard deviation
    std::vector<size_t> n_air_gen(N_gen,0);//no of airways in each generation
    //std::cout << this->transport_tree->count_airways() << std::endl;
    for (int J = 0; J < this->transport_tree->count_airways(); ++J) {
        size_t gen = this->transport_tree->get_weibel_order(this->transport_tree->get_edges_in_airway(J,0));
        if (gen < N_gen)
        {
            n_air_gen[gen] ++;
        }
    }
    //for (int i = 0; i < N_gen; ++i) {
    //    std::cout << i << ' ' << n_air_gen[i] << std::endl;
    //}
    std::vector<double> mean_dep(N_gen);
    for(int i = 0; i < N_gen; i++)
    {
        mean_dep[i] = gen_dep[i] / n_air_gen[i];
    }

    std::vector<double> gen_sd(N_gen,0);
    Eigen::VectorXd edge_airway_dep = this->transport_tree->get_all_airway_depositions();
    for (int J = 0; J < this->transport_tree->count_airways(); ++J)
    {
        size_t gen = this->transport_tree->get_weibel_order(this->transport_tree->get_edges_in_airway(J,0));
        if (gen < N_gen)
        {
            gen_sd[gen] += std::pow(mean_dep[gen] - edge_airway_dep[this->transport_tree->get_edges_in_airway(J,0)],2);
        }
        //if (J == 1){std::cout << gen << std::endl;}
    }
    //std::cout << this->transport_tree->get_generational_deposition(this->transport_tree->get_edges_in_airway(1,0)) << ' ' << gen_dep[1] << ' ' << mean_dep[1] << ' ' << edge_airway_dep[this->transport_tree->get_edges_in_airway(1,0)] << std::endl;
    for (int i = 0; i < N_gen; ++i) {
        gen_sd[i] = std::sqrt(gen_sd[i] / n_air_gen[i]);
        //std::cout << i << ' ' << gen_sd[i] << std::endl;
    }


    for (int i = 0; i < N_gen; ++i) {
        if(i < N_gen - 1) {
            output << gen_sd[i] / mean_dep[i] << ", ";
        }
        else
        {
            output << gen_sd[i] / mean_dep[i] << ", ";
        }
    }

    //Flux for all 6th generation airways
    Eigen::VectorXd edge_airway_dep = this->transport_tree->get_all_airway_depositions();
    for (int J = 0; J < this->transport_tree->count_airways(); ++J) {
        if (this->transport_tree->get_weibel_order(this->transport_tree->get_edges_in_airway(J,0)) == 6)
        {
            output << this->transport_tree->get_edge_flux(this->transport_tree->get_edges_in_airway(J,0)) << ", ";
            output << edge_airway_dep[this->transport_tree->get_edges_in_airway(J,0)] << ", ";
        }
    }
*/
    output << 0 << '\n';

	output.close();

	//std::cout << 0 << std::endl;
}

int LungSimulation::print_acinus_vtk(const std::string & filename,
										const int & tstep, const bool & original)
{
	//complete filename
	std::stringstream ss;
	ss << filename << ".vtk";
	std::string fname = ss.str().c_str();
	//open file
	std::ofstream output;
	output.open(fname);
	//define output scale factors
	double vol_scale = 1.0 / this->params->get_conversion(VOL_KEY);
	double length_scale = 1.0 / this->params->get_conversion(LENGTH_KEY);
	output << std::fixed << std::setprecision(OUTPUT_PRECISION_DEFAULT);
	
	//vtk headers
	size_t Npts =  this->flow_tree->count_term_nodes();
	output <<  "# vtk DataFile Version 3.0\nNetwork_data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	output <<  "POINTS " << Npts << " float\n";
	//print position vectors
	for (size_t kt = 0; kt < Npts; kt++)
	{
		size_t k = this->flow_tree->get_first_term_index() + kt;
		output << this->flow_tree->get_node(k)->get_pos(0) * length_scale << ' ' <<
				  this->flow_tree->get_node(k)->get_pos(1) * length_scale << ' ' <<
				  this->flow_tree->get_node(k)->get_pos(2) * length_scale << '\n';
	}
	
	//get volumes and concs
	output <<  "\n\nPOINT_DATA " << Npts << "\n";
	output <<  "\nSCALARS Volume float\nLOOKUP_TABLE default\n";
	
	for (size_t kt = 0; kt < this->transport_tree->count_acini(); kt++)
	{
		output << this->transport_tree->get_acinus_volume(kt) * vol_scale << '\n';
	}
	
	output <<  "\nSCALARS IGVolume float\nLOOKUP_TABLE default\n";
	for (size_t kt = 0; kt < this->transport_tree->count_acini(); kt++)
	{
		output << this->transport_tree->get_acinus_IG_volume(kt, 
					this->options->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value()) * vol_scale << '\n';
	}
	
	output <<  "\nSCALARS Conc float\nLOOKUP_TABLE default\n";
	for (size_t kt = 0; kt < this->transport_tree->count_acini(); kt++)
	{
		output << this->transport_tree->get_acinus_IG_volume(kt, this->options->get_option<bool>(SIMULATE_DEPOSITION_KEY)->get_value()) 
					/ this->transport_tree->get_acinus_volume(kt) << '\n';
	}
	
	output.close();
	
	//vtkSmartPointer<vtkMultiBlockDataSet> data =
	//vtkSmartPointer<vtkMultiBlockDataSet>::New();
	//for(size_t k = this->flow_tree->get_first_term_index(); k < this->flow_tree->count_nodes(); k++)
	//{
	//	size_t kt = k - this->flow_tree->get_first_term_index();
	//	double scale = this->params->get_conversion(LENGTH_KEY);
	//	network::Position pos = this->flow_tree->get_node(k)->get_original_pos() / scale;
	//	auto sphere = vtkSmartPointer<vtkSphereSource>::New();
	//	sphere->SetCenter(pos.x[0], pos.x[1], pos.x[2]);
	//	double rad = pow(0.75*this->flow_tree->get_termnode_volume(kt)/M_PI,1.0/3.0)/scale;
 // 		sphere->SetRadius(rad);
	//	// Make the surface smooth.
	//	sphere->SetPhiResolution(6);
	//	sphere->SetThetaResolution(6);
	//	sphere->Update();
	//	//map to poly data
	//	auto sphdata = vtkSmartPointer<vtkPolyData>::New();
	//	sphdata->ShallowCopy(sphere->GetOutput());
	//	auto conc = vtkSmartPointer<vtkDoubleArray>::New();
	//	double concval[1] = {this->transport_tree->get_acinus_IG_volume(kt)
	//		/this->flow_tree->get_termnode_volume(kt)};
	//	conc->SetNumberOfComponents(1);
	//	conc->SetName("Concentration");
	//	conc->InsertNextTuple(concval);
	//	sphdata->GetFieldData()->AddArray(conc);
	//	data->SetBlock(unsigned(kt), sphdata);
	//	

	//	/*auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	//	

	//	std::stringstream ss;
	//	ss << filename << "_" << kt << "_" << tstep << ".vtp";
	//	writer->SetFileName(ss.str().c_str());
	//	writer->SetInputData(sphdata);
	//	writer->SetDataModeToBinary();
	//	writer->Update();
	//	auto test = writer->GetInput();
	//	cout << test->GetNumberOfCells() << "\n\n";
	//	writer->Write();*/
	//}


	//vtkSmartPointer<vtkXMLMultiBlockDataWriter> output = 
	//	vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
	//std::stringstream ss;
	//ss << filename << "_" << tstep << ".vtm";
	//output->SetFileName(ss.str().c_str());
	//output->SetInputData(data);
	//output->SetDataModeToBinary();
	//output->Update();
	//output->Write();

	return 0;
}
