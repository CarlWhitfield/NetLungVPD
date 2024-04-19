#ifndef AIRWAY_FLOW_NETWORK_H
#define AIRWAY_FLOW_NETWORK_H

#include<network_3D.h>
#include "airway_edges_nodes.h"
#include "fv_element_tree.h"
#include "list_template.h"
#include "define_options.h"
#include "define_params.h"
#include "elastance.h"
#include "resistance.h"
#include<map>
#include<Eigen/Sparse>

#define FN_NODE_COUNT 1  //number of node variables
//node vals
#define FN_NODE_PRESSURE 0
const std::string fn_node_strings[] = {"Pressure_cmH20"};

#define FN_EDGE_COUNT 2  //number of edge variables
//edge vals
#define FN_EDGE_RESISTANCE 0
#define FN_EDGE_FLUX 1
const std::string fn_edge_strings[] = {"Resistance_cmH20s_per_L", "Flux_L_per_s"};


#define FN_TNODE_COUNT 6  //number of term node variables
//term_node_vals
#define FN_TNODE_VOLUME 0
#define FN_TNODE_ELASTANCE 1
#define FN_TNODE_RESISTANCE 2
#define FN_TNODE_INERTANCE 3
#define FN_TNODE_PRESSGRAD 4
#define FN_TNODE_INTRINSIC_VOL 5

inline void output_Eigen_error_message(const int & info)
{
	switch(info)
	{
	case Eigen::Success:
		{
			std::cout << "Success\n";
		} break;

	case Eigen::NumericalIssue:
		{
			std::cout << "Numerical Issue\n";
		} break;

	case Eigen::NoConvergence:
		{
			std::cout << "No convergence\n";
		} break;

	case Eigen::InvalidInput:
		{
			std::cout << "Invalid input\n";
		} break;

	default:
		{
			std::cout << "Did not recognise error meassage.\n";
		}
	}
}

class AirwayFlowNetwork;

class PressureVolumeSolver
{
protected:
	AirwayFlowNetwork *tree;        //pointer to airway network being simulated
	double dt;
	Eigen::VectorXd Vold;
public:
	PressureVolumeSolver(){};
	PressureVolumeSolver(AirwayFlowNetwork *t, const double & dtime);
	void update_flow(const Eigen::VectorXd & Vterm, const Eigen::VectorXd & Vtermold);
	
	virtual void setup_flow_calc(const bool & Flow_BC){};
	virtual void update_press_and_vol(const size_t &){};
	virtual double get_expansion_volume(const size_t &){ return 0; }
};

class AirwayFlowNetwork: public AirwayNetwork
{
private:
	std::vector<Eigen::VectorXd> node_vals, edge_vals, term_node_vals;
	double Ppl, mouth_resistance;
	std::shared_ptr<ResistanceObject> res_calc;
	std::shared_ptr<ElastanceObject> e_calc;
	std::shared_ptr<PressureVolumeSolver> pv_solver;
	inline void init_var_vectors()
	{
		this->node_vals.resize(FN_NODE_COUNT);
		for(size_t n = 0; n < FN_NODE_COUNT; n++)
		{
			this->node_vals[n] = Eigen::VectorXd::Zero(this->count_nodes());
		}

		this->edge_vals.resize(FN_EDGE_COUNT);
		for(size_t n = 0; n < FN_EDGE_COUNT; n++)
		{
			this->edge_vals[n] = Eigen::VectorXd::Zero(this->count_edges());
		}

		this->term_node_vals.resize(FN_TNODE_COUNT);
		for(size_t n = 0; n < FN_TNODE_COUNT; n++)
		{
			this->term_node_vals[n] = Eigen::VectorXd::Zero(this->count_term_nodes());
		}
	}
	void setup_gravity_dist(const char & option, const network::Position & grav_dir,
							const double & gradient, const double & Ntot, 
							std::vector<double> & grav_dist);
	void setup_acinus_vols_dist(const bool & option, const double & bag_vol, 
								const double & Ntot, std::vector<double> & acinus_vols);
public:

	AirwayFlowNetwork():AirwayNetwork(){};
	AirwayFlowNetwork(std::map<long int, std::shared_ptr<AirwayNetworkNode>> &node, 
		              std::map<long int, std::shared_ptr<AirwayNetworkEdge>> &edge)
		              :AirwayNetwork(node, edge){};
	AirwayFlowNetwork(std::map<long int, std::shared_ptr<AirwayNetworkNode>> &node, 
		              std::map<long int, std::shared_ptr<AirwayNetworkEdge>> &edge, 
		              SimOptionList *o, SimParameterList *p);
	AirwayFlowNetwork(const std::string & node_filename, const std::string & branch_filename, 
		              const std::string & termnode_filename,  const double & l_scale)
					  :AirwayNetwork(node_filename, branch_filename, termnode_filename, l_scale){};

	void setup_afn(SimOptionList *, SimParameterList*);  
	void assign_params_to_flow_network(SimOptionList *, SimParameterList*);                                                  //print vtk file of nodes and branches
	//int print_network_for_input(const std::string &filehead);                                  //print nodes and branches for reading back into code
	std::shared_ptr<PressureVolumeSolver> create_pv_solver(SimOptionList *, SimParameterList*);

	inline double get_Ppl() const { return (this->Ppl); };
	inline double get_mouth_resistance() const { return (this->mouth_resistance); }
	inline double get_node_pressure(const size_t & k) const { return (this->node_vals[FN_NODE_PRESSURE][k]); }
	inline double get_edge_flux(const size_t & j) const { return (this->edge_vals[FN_EDGE_FLUX][j]); }
	inline double get_edge_resistance(const size_t & j) const { return (this->edge_vals[FN_EDGE_RESISTANCE][j]); }
	inline double get_edge_resistance_for_flux(const size_t & j, const double & flux) const 
	{ 
		return (this->res_calc->calc_resistance(this->EdgeVec[j]->get_geom(), flux) / this->EdgeVec[j]->branch_count());
	}
	inline double get_edge_resistance_grad_for_flux(const size_t & j, const double & flux) const 
	{ 
		return (this->res_calc->resistance_flux_grad(this->EdgeVec[j]->get_geom(), flux) / this->EdgeVec[j]->branch_count());
	}
	inline double get_termnode_volume(const size_t & kt) const { return (this->term_node_vals[FN_TNODE_VOLUME][kt]); }
	inline double get_termnode_resistance(const size_t & kt) const { return (this->term_node_vals[FN_TNODE_RESISTANCE][kt]); }
	inline double get_termnode_inertance(const size_t & kt) const { return (this->term_node_vals[FN_TNODE_INERTANCE][kt]); }
	inline double get_termnode_elastance(const size_t & kt) const { return (this->term_node_vals[FN_TNODE_ELASTANCE][kt]); }
	inline double get_termnode_elastance_for_vol(const size_t & kt, const double & vol) const 
	{ 
		double ext = pow(vol / (0.5*this->term_node_vals[FN_TNODE_INTRINSIC_VOL][kt]), 1.0/3.0);
		return (this->e_calc->calc_elastance(ext, this->term_node_vals[FN_TNODE_ELASTANCE][kt]));
	}
	inline double get_termnode_elastance_vol_grad(const size_t & kt, const double & vol) const
	{ 
		double ext = pow(vol / (0.5*this->term_node_vals[FN_TNODE_INTRINSIC_VOL][kt]), 1.0/3.0);
		return (this->e_calc->elastance_ext_grad(ext, this->term_node_vals[FN_TNODE_ELASTANCE][kt]))/vol;
	}

	inline double get_termnode_intrinsic_volume(const size_t & kt) const { return (this->term_node_vals[FN_TNODE_INTRINSIC_VOL][kt]); }
	inline double get_termnode_relative_pressure(const size_t & kt) const { return (this->term_node_vals[FN_TNODE_PRESSGRAD][kt]); }

	//these return the Eigen vectors by const reference for efficiency
	inline const Eigen::VectorXd& get_all_node_pressures() const { return this->node_vals[FN_NODE_PRESSURE]; }
	inline const Eigen::VectorXd& get_all_edge_fluxes() const{ return this->edge_vals[FN_EDGE_FLUX]; }
	inline const Eigen::VectorXd& get_all_edge_resistances() const{ return this->edge_vals[FN_EDGE_RESISTANCE]; }
	inline const Eigen::VectorXd& get_all_termnode_volumes() const { return this->term_node_vals[FN_TNODE_VOLUME]; }
	inline const Eigen::VectorXd& get_all_termnode_resistances() const { return this->term_node_vals[FN_TNODE_RESISTANCE]; }
	inline const Eigen::VectorXd& get_all_termnode_inertances() const { return this->term_node_vals[FN_TNODE_INERTANCE]; }
	inline const Eigen::VectorXd& get_all_termnode_elastances() const { return this->term_node_vals[FN_TNODE_ELASTANCE]; }
	inline const Eigen::VectorXd& get_all_termnode_intrinsic_volumes() const { return this->term_node_vals[FN_TNODE_INTRINSIC_VOL]; }
	inline const Eigen::VectorXd& get_all_termnode_relative_pressures() const { return this->term_node_vals[FN_TNODE_PRESSGRAD]; }
	
	inline double sum_all_termnode_volumes() const { return this->term_node_vals[FN_TNODE_VOLUME].sum(); }
	inline double get_lung_expansion_volume(const size_t & it) const { return this->pv_solver->get_expansion_volume(it); }

	inline void solve_pressure_volume_update(const size_t & i){ this->pv_solver->update_press_and_vol(i); }
	inline void update_Ppl(const double & p){ this->Ppl = p; }
	inline void update_node_pressure(const size_t & k, const double & p){ this->node_vals[FN_NODE_PRESSURE][k] = p; }
	inline void update_edge_flux(const size_t & j, const double & f){ this->edge_vals[FN_EDGE_FLUX][j] = f; }
	inline void update_termnode_volume(const size_t & kt, const double & v) { this->term_node_vals[FN_TNODE_VOLUME][kt] = v; }
	inline void update_termnode_resistance(const size_t & kt, const double & r){ this->term_node_vals[FN_TNODE_RESISTANCE][kt] = r; }
	inline void update_termnode_inertance(const size_t & kt, const double & i){ this->term_node_vals[FN_TNODE_INERTANCE][kt] = i; }
	inline void update_termnode_elastance(const size_t & kt, const double & e){ this->term_node_vals[FN_TNODE_ELASTANCE][kt] = e; }
	inline void update_termnode_intrinsic_volume(const size_t & kt, const double & iv){ this->term_node_vals[FN_TNODE_INTRINSIC_VOL][kt] = iv; }
	inline void update_termnode_relative_pressure(const size_t & kt, const double & rp){ this->term_node_vals[FN_TNODE_PRESSGRAD][kt] = rp; }
	
	inline void update_all_node_pressures(const Eigen::VectorXd & p){ this->node_vals[FN_NODE_PRESSURE] = p; }
	inline void update_all_edge_fluxes(const Eigen::VectorXd & f){ this->edge_vals[FN_EDGE_FLUX] = f; }
	inline void update_all_termnode_volumes(const Eigen::VectorXd & v) { this->term_node_vals[FN_TNODE_VOLUME] = v; }
	inline void update_all_termnode_resistances(const Eigen::VectorXd & r){ this->term_node_vals[FN_TNODE_RESISTANCE] = r; }
	inline void update_all_termnode_inertances(const Eigen::VectorXd & i){ this->term_node_vals[FN_TNODE_INERTANCE] = i; }
	inline void update_all_termnode_elastances(const Eigen::VectorXd & e){ this->term_node_vals[FN_TNODE_ELASTANCE] = e; }
	inline void update_all_termnode_intrinsic_volumes(const Eigen::VectorXd & iv){ this->term_node_vals[FN_TNODE_INTRINSIC_VOL] = iv; }
	inline void update_all_termnode_relative_pressures(const Eigen::VectorXd & rp){ this->term_node_vals[FN_TNODE_PRESSGRAD] = rp; }

	inline void update_all_edge_resistances_by_rule()
	{
		for(size_t j = 0; j < this->count_edges(); j++)
		{
			this->edge_vals[FN_EDGE_RESISTANCE][j] = this->res_calc->calc_resistance(this->EdgeVec[j]->get_geom(), this->edge_vals[FN_EDGE_FLUX][j]) / this->get_edge(j)->branch_count();
		}
	}

	int print_flow_vtk(const std::string & filename, SimParameterList *p);
};

std::vector<double> read_flow_file(const std::string & , const double &dt, 
								   const double & sec_to_simtime, const double & L_to_simvol, 
								   inlist::Parameter<double> * stop_time);

class FlowObject
{
protected:
	double dt;
public:
	FlowObject(const double & dtin)
	{
		this->dt = dtin;
	}
	virtual double get_volume(const size_t & i_time){ return 0; }
};

class FileFlowObject: public FlowObject
{
private:
	std::vector<double> inflation;
	void parse_flow_input(const std::vector<double> & flow, const double & infl0 = 0);

public:
	FileFlowObject(const std::vector<double> & flow, const double & dtin):FlowObject(dtin)
	{
		this->parse_flow_input(flow);
	}

	double get_volume(const size_t & i_time);
};

class SinFlowObject: public FlowObject
{
private:
	double tau;

public:
	SinFlowObject(const double & dtin, const double & time_sig):FlowObject(dtin)
	{
		this->tau = time_sig;
	}

	double get_volume(const size_t & i_time);
};

class StepFlowObject: public FlowObject
{
private:
	double tau;

public:
	StepFlowObject(const double & dtin, const double & time_sig):FlowObject(dtin)
	{
		this->tau = time_sig;
	}

	double get_volume(const size_t & i_time);
};

class LinearSinusoidalPressureVolumeSolver: public PressureVolumeSolver
{
	protected:
		double Pp0, Ppsin, Ppcos, Q0sin, Q0cos, VT, tau, VFRC, DP;
		Eigen::VectorXd P0, Psin, Pcos, V0, Vsin, Vcos;

	public:
		//LinearSinusoidalPressureVolumeSolver():PressureVolumeSolver(){};
		//LinearSinusoidalPressureVolumeSolver(AirwayNetwork * t):PressureVolumeSolver(t){};
		LinearSinusoidalPressureVolumeSolver(AirwayFlowNetwork * t,  const double & TidalVol, 
			const double & PleuralDPressure, const double & TidalTime, const double & VFRC,
			const double & dt):PressureVolumeSolver(t, dt)
		{
			this->VFRC = VFRC;
			this->VT = TidalVol;
			this->tau = TidalTime;
			this->DP = PleuralDPressure;
		};
		void setup_flow_calc(const bool & Flow_BC);
		void update_press_and_vol(const size_t & it);
		double get_expansion_volume(const size_t & it);
};

class LinearGenericPressureVolumeSolver: public PressureVolumeSolver
{
private:
	Eigen::SparseMatrix<double> A;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> LU_solver;
protected:
	double VT, tau, V0;
	char flow_opt;
	size_t mat_dim_size;
	Eigen::VectorXd Vold_old, x;
	std::vector<Eigen::VectorXd> breath_sols;
	bool breaths_converged;

	std::shared_ptr<FlowObject> lung_volume_model;
	std::vector<double> flow_in;
	void initialise(Eigen::SparseMatrix<double> & mat);
	void run_breaths(const int & nbreaths);
	void update_from_solution();
public:
	LinearGenericPressureVolumeSolver(AirwayFlowNetwork * t, const double & TidalVol, const double & TidalTime, 
			                    const double & dt, const char & flow_input_option):PressureVolumeSolver(t, dt)
	{
		this->VT = TidalVol;
		this->tau = TidalTime;
		this->flow_opt = flow_input_option;
		this->breaths_converged = false;
	};
	LinearGenericPressureVolumeSolver(AirwayFlowNetwork * t,  const double & dt, 
			                    const std::vector<double> & mouth_flux):PressureVolumeSolver(t, dt)
	{
		this->flow_opt = FLOW_FROM_FILE;
		this->flow_in = mouth_flux;
		this->breaths_converged = false;
	}

	void parse_flow_input_option();
	void run_breaths_until_converged();
	virtual void setup_flow_calc(const bool & Flow_BC);
	virtual void update_press_and_vol(const size_t & i_time);
	void update_flow_equations(Eigen::VectorXd & func, const size_t & i_time);
	void fill_matrix(Eigen::SparseMatrix<double> & mat);
	double get_expansion_volume(const size_t & i_time);
};

class GenericPressureVolumeSolver : public LinearGenericPressureVolumeSolver
{
private:
	Eigen::SparseMatrix<double> Jacobian;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicgstab_solver;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> LU_solver;
protected:
	double tol;
public:
	GenericPressureVolumeSolver(AirwayFlowNetwork * t, const double & TidalVol, const double & TidalTime, 
			                    const double & dt, const char & flow_input_option)
								:LinearGenericPressureVolumeSolver(t, TidalVol, TidalTime, dt, flow_input_option)
	{
		this->tol = 1E-12;
	}
	GenericPressureVolumeSolver(AirwayFlowNetwork * t,  const double & dt, 
			                    const std::vector<double> & mouth_flux):LinearGenericPressureVolumeSolver(t, dt, mouth_flux)
	{
		this->tol = 1E-12;
	}

	void setup_flow_calc(const bool & Flow_BC);
	void update_press_and_vol(const size_t & i_time);
	void update_jacobian();
};

#endif
