#ifndef LUNG_SIMULATION_H
#define LUNG_SIMULATION_H

#include "list_template.h"
#include "network_3D.h"
#include "airway_flow_network.h"
#include "airway_transport_network.h"
#include "define_options.h"
#include "define_params.h"

//class for a simulation
class LungSimulation
{
private:
	double time;
	std::map<std::string,std::string> tree_files;
	std::map<std::string,std::string> acin_files;
public:
	//--------Member variables-------//
	std::shared_ptr<AirwayFlowNetwork> flow_tree;        //pointer to airway network being simulated
	std::shared_ptr<AirwayTransportNetwork> transport_tree;
	std::shared_ptr<AirwayNetwork> acin_template;
	std::shared_ptr<SimOptionList> options;        //pointer to option list object
	std::shared_ptr<SimParameterList> params;      //pointer to parameter list object

	//-------Member functions-------//
	LungSimulation(){};                     //default constructor(blank)
	LungSimulation(const int &, char**);    //constructor from command line args
	int read_acin_template();
	int initialise();
	int create_flow_network();
	int import_network();   //import network from branch, node and termnode files -- update airway dead-space based on inputs
	int simulate();

	void sort_files();
	//void create_tree_maps();
	void copy_fluxes_to_transport_network();

	void print_acinus_headers(const std::string &);
	void print_washout_file_headers(const std::string &filename);
	void append_acinus_vals(const std::string &, const double & t);
	void append_washout_file(const std::string &, const double & t);
	int print_acinus_vtk(const std::string & filename,  const int & tstep, 
		                 const bool & original = false);
};

#endif