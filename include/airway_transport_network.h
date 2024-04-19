#ifndef AIRWAY_TRANSPORT_NETWORK_H
#define AIRWAY_TRANSPORT_NETWORK_H

#include "airway_edges_nodes.h"
#include "fv_element_tree.h"
#include "acinus_models.h"
#include "compartmental.h"
#include "diffusivity.h"
#include <chrono>
#include <vector>
#include <map>
#include <Eigen/Sparse>
//#include <boost/math/quadrature/trapezoidal.hpp>

#define TN_NODE_COUNT 14
//node vals
#define TN_NODE_CONC 0
#define TN_NODE_DEPOSITION 1
#define TN_NODE_VOLUME 2
#define TN_INNER_NODE_VOLUME 3
#define TN_NODE_CUMULATIVE_DEPOSITION 4
#define TN_NODE_CUMULATIVE_DIFFUSIVE_DEPOSITION 5
#define TN_NODE_CUMULATIVE_SEDIMENTATION_DEPOSITION 6
#define TN_NODE_CUMULATIVE_IMPACTION_DEPOSITION 7
#define TN_NODE_CUMULATIVE_ALVEOLAR_DEPOSITION 8
#define TN_NODE_TOTAL_DEP_OVER_VOL 9
#define TN_NODE_BRANGLE 10
#define TN_NODE_GENERATIONAL_DEPOSITION 11
#define TN_NODE_PENDELLUFT 12
#define TN_NODE_SPECT_DEPOSITION 13
const std::string tn_node_strings[] = {"Concentration", "Deposition_Rate_L_per_s",  "Volume_L", "Inner_Volume_L",
	"Cumulative_Deposition_L", "Cumulative_Diffusion_Deposition_L", "Cumulative_Sedimentation_Deposition_L", 
	"Cumulative_Impaction_deposition_L","Cumulative_Alveolar_deposition_L", "Cumulative_Deposition_Over_Node_Volume",
	"Branching_Angle","Generational_deposition","Cumul_Pendelluft_FLux", "Spect_deposition"};
//note that currently concentration is on a dimensionless scale between 0 and 1, with c = 1 corresponding to
//the inhaled concentration. So deposition rate is Volume per second rather than mass per second. To convert
//to mass per second I think you would just have to scale everything by whatever the mass of particles in one
//unit of inhaled gas is.

#define TN_EDGE_COUNT 12
//edge vals
#define TN_EDGE_FLUX 0
#define TN_EDGE_DIFFUSIVITY 1
#define TN_EDGE_GRAVITY_ANGLE 2
#define TN_EDGE_AIRWAY_DEPOSITION 3
#define TN_EDGE_AIRWAY_IMP_DEPOSITION 4
#define TN_EDGE_AIRWAY_DIFF_DEPOSITION 5
#define TN_EDGE_AIRWAY_SED_DEPOSITION 6
#define TN_EDGE_AIRWAY_NO 7
#define TN_EDGE_GENERATION_NO 8
#define TN_EDGE_MAX_FLUX 9
#define TN_EDGE_AIRWAY_LENGTH 10
#define TN_EDGE_HORSFIELD_ORDER 11
const std::string tn_edge_strings[] = {"Flux_L_per_s", "Volumetric_diffusivity_L_per_s",
                                       "Gravity_angle_(rad)", "Total_deposition_in_whole_airway_L",
                                       "Impaction_deposition_in_whole_airway_L","Diffusion_deposition_in_whole_airway_L",
                                       "Sedimentation_deposition_in_whole_airway_L",
                                       "Airway_number", "Generation", "Max_flux_L_per_s",
                                       "Airway_length_m", "Horsfield_order"};

class AirwayTransportNetwork;

class TransportSolver
{
protected:
	AirwayTransportNetwork *tree;        //pointer to airway network being simulated
	Eigen::VectorXd fv_vol_old, inner_fv_vol_old;
    bool simulate_deposition, simulate_alveolar_deposition, imp_both_dirs_option;
    double particle_size, particle_density, grav_accel, kB_T_mu, air_viscosity, mean_free_path, air_density;
    char diff_formula_option, sed_formula_option, imp_formula_option;
public:
	TransportSolver(){};
	TransportSolver(AirwayTransportNetwork *t);

	virtual void setup_transport_calc(){};
	virtual void update_concentration(const double & dt){};
	inline void update_vol_old(const Eigen::VectorXd & Vol)
	{
		this->fv_vol_old = Vol;
	}
	inline double get_node_vol_old(const size_t & k) const
	{
		return this->fv_vol_old[k];
	}
	inline const Eigen::VectorXd& get_all_node_vols_old() const
	{
		return this->fv_vol_old;
	}
	inline void update_inner_vol_old(const Eigen::VectorXd & Vol)
	{
		this->inner_fv_vol_old = Vol;
	}
	inline double get_inner_node_vol_old(const size_t & k) const
	{
		return this->inner_fv_vol_old[k];
	}
	inline const Eigen::VectorXd& get_all_inner_node_vols_old() const
	{
		return this->inner_fv_vol_old;
	}
    inline void set_deposition_params(const double & prad, const double & pdens, const double & ga,
                                      const double & kBTmu, const double & mu_air, const double & lambda,
                                      const double & rho_air, const char & diff_option, const char & sed_option, const char & imp_option,
                                      const bool & imp_dirs_option)
    {
        this->particle_size = prad;
        this->particle_density = pdens;
        this->grav_accel = ga;
        this->kB_T_mu = kBTmu;
        this->air_viscosity = mu_air;
        this->air_density = rho_air;
        this->mean_free_path = lambda;
        this->diff_formula_option = diff_option;
        this->sed_formula_option = sed_option;
        this->imp_formula_option = imp_option;
        this->imp_both_dirs_option = imp_dirs_option;
    }
};

class AirwayTransportNetwork: public FVElementNetwork
{
private:
	std::vector<Eigen::VectorXd> node_vals, edge_vals;
	std::vector<size_t> conducting_nodes, conducting_edges, edge_airway_no, no_airway_edges;
	std::vector<std::shared_ptr<Acinus>> acinus_vector;
	std::vector<std::vector<size_t>> acin_node_map, acin_edge_map, airway_to_edge_map;
	bool washin_mode, source_mode;
	std::shared_ptr<DSVolume> Mouth;                         //stores mouth dead-space object (see mouth_volume.h)
	std::shared_ptr<DiffusivityObject> conducting_diffusivity, acinar_diffusivity;
	std::shared_ptr<TransportSolver> transport_solver;
    std::map<std::pair<size_t,size_t>,double> branch_angles;  //maps edge pair to angle (needs lower number first)
    size_t airway_count;
	bool print_acinar_airways;

	void init_var_vectors();

public:	
	AirwayTransportNetwork():FVElementNetwork()
	{
		washin_mode = true;
		source_mode = true;
	};
	AirwayTransportNetwork(AirwayNetwork *, const Eigen::VectorXd &, SimOptionList*, SimParameterList*,
                           AirwayNetwork *acin_template=NULL);
	inline bool get_washin_mode(){ return washin_mode; }
	inline bool get_source_mode(){ return source_mode; }
	inline void set_washin_mode(const bool & m){ washin_mode = m; }
	inline void set_source_mode(const bool & m){ source_mode = m; }
	inline double mouth_bc()
	{
		if(washin_mode && source_mode) return 1.0;
		else return 0.0;
	}

	bool iterate_node_map(Eigen::SparseMatrix<double>::InnerIterator &);
	bool iterate_edge_map(Eigen::SparseMatrix<int>::InnerIterator &);

	inline size_t count_acini(){ return this->acinus_vector.size(); }
	inline DSVolume* get_mouth_object(){ return Mouth.get(); }
	inline double get_node_conc(const size_t & k) const { return (this->node_vals[TN_NODE_CONC][k]); }
	inline double get_edge_flux(const size_t & j) const { return (this->edge_vals[TN_EDGE_FLUX][j]); }
    inline double get_edge_gravity_angle(const size_t & j) const { return (this->edge_vals[TN_EDGE_GRAVITY_ANGLE][j]); }
	inline double get_edge_diffusivity(const size_t & j) const { return (this->edge_vals[TN_EDGE_DIFFUSIVITY][j]); }
    inline double get_node_fv_volume(const size_t & k) const { return (this->node_vals[TN_NODE_VOLUME][k]); }
	inline double get_node_inner_fv_volume(const size_t & k) const { return (this->node_vals[TN_INNER_NODE_VOLUME][k]); }
    inline double get_molecular_diffusivity() const {return this->conducting_diffusivity->get_molecular_diffusivity(); }
    inline double get_total_deposition(const size_t & k) const { return (this->node_vals[TN_NODE_CUMULATIVE_DEPOSITION][k]); }
    inline double get_total_diff_deposition(const size_t & k) const { return (this->node_vals[TN_NODE_CUMULATIVE_DIFFUSIVE_DEPOSITION][k]); }
    inline double get_total_sed_deposition(const size_t & k) const { return (this->node_vals[TN_NODE_CUMULATIVE_SEDIMENTATION_DEPOSITION][k]); }
    inline double get_total_imp_deposition(const size_t & k) const { return (this->node_vals[TN_NODE_CUMULATIVE_IMPACTION_DEPOSITION][k]); }
	inline double get_total_alv_deposition(const size_t & k) const { return (this->node_vals[TN_NODE_CUMULATIVE_ALVEOLAR_DEPOSITION][k]); }
    inline double get_generational_deposition(const size_t & k) const {return (this->node_vals[TN_NODE_GENERATIONAL_DEPOSITION][k]);}
    inline const Eigen::VectorXd& get_all_node_pendelluft() const { return (this->node_vals[TN_NODE_PENDELLUFT]); }
    inline double get_spect_deposition(const size_t & k) const { return (this->node_vals[TN_NODE_SPECT_DEPOSITION][k]); }

	inline const Eigen::VectorXd& get_all_node_concs() const { return (this->node_vals[TN_NODE_CONC]); }
	inline const Eigen::VectorXd& get_all_edge_diffusivities() const { return (this->edge_vals[TN_EDGE_DIFFUSIVITY]); }
    inline const Eigen::VectorXd& get_all_edge_fluxes() const { return (this->edge_vals[TN_EDGE_FLUX]); }
    inline const Eigen::VectorXd& get_all_airway_depositions() const { return (this->edge_vals[TN_EDGE_AIRWAY_DEPOSITION]); }
    inline const Eigen::VectorXd& get_all_airway_imp_depositions() const { return (this->edge_vals[TN_EDGE_AIRWAY_IMP_DEPOSITION]); }
    inline const Eigen::VectorXd& get_all_airway_diff_depositions() const { return (this->edge_vals[TN_EDGE_AIRWAY_DIFF_DEPOSITION]); }
    inline const Eigen::VectorXd& get_all_airway_sed_depositions() const { return (this->edge_vals[TN_EDGE_AIRWAY_SED_DEPOSITION]); }
    inline const Eigen::VectorXd& get_all_airway_numbers() const { return (this->edge_vals[TN_EDGE_AIRWAY_NO]); }
    inline const Eigen::VectorXd& get_all_edge_max_fluxes() const { return (this->edge_vals[TN_EDGE_MAX_FLUX]); }
    inline const Eigen::VectorXd& get_all_airway_generations() const { return (this->edge_vals[TN_EDGE_GENERATION_NO]); }
    inline const Eigen::VectorXd& get_all_horsfield_orders() const { return (this->edge_vals[TN_EDGE_HORSFIELD_ORDER]); }
	inline const Eigen::VectorXd& get_all_airway_lengths() const { return (this->edge_vals[TN_EDGE_AIRWAY_LENGTH]); }
    inline const Eigen::VectorXd& get_all_node_volumes() const { return (this->node_vals[TN_NODE_VOLUME]); }
	inline const Eigen::VectorXd& get_all_inner_node_volumes() const { return (this->node_vals[TN_INNER_NODE_VOLUME]); }

    void update_conducting_diffusivities();
	void update_acinar_diffusivities();
    void update_branching_angles();     //angles between branches
    void update_node_volumes();  //outer volume
	void update_inner_node_volumes();  //inner volume

	inline void update_edge_flux(const size_t & j, const double & flux){ this->edge_vals[TN_EDGE_FLUX][j] = flux; }
	inline void update_node_conc(const size_t & k, const double & conc){ this->node_vals[TN_NODE_CONC][k] = conc; }
	void update_acinar_volume(const size_t & kt, const double & v);

	void compute_all_fluxes(const double & dt);
    void compute_all_edge_gravity_angles(const pos::Position<double> grav_dir);
	inline void update_all_edge_fluxes(const Eigen::VectorXd & flux){
        this->edge_vals[TN_EDGE_FLUX] = flux;
    }
	inline void update_all_node_concs(const Eigen::VectorXd & conc){
        this->node_vals[TN_NODE_CONC] = conc;
    }
    inline void update_all_node_depositions(const Eigen::VectorXd & dep){
        this->node_vals[TN_NODE_DEPOSITION] = dep;
    }
    inline void update_all_volumes(const Eigen::VectorXd & vol){
        this->node_vals[TN_NODE_VOLUME] = vol;
    }
	inline void update_all_inner_volumes(const Eigen::VectorXd & vol){
        this->node_vals[TN_INNER_NODE_VOLUME] = vol;
    }
    inline void update_all_cumul_depositions(const Eigen::VectorXd & dep){
        this->node_vals[TN_NODE_CUMULATIVE_DEPOSITION] = this->node_vals[TN_NODE_CUMULATIVE_DEPOSITION] + dep;
    }
    inline void update_all_cumul_diff_depositions(const Eigen::VectorXd & diff_dep){
        this->node_vals[TN_NODE_CUMULATIVE_DIFFUSIVE_DEPOSITION] = this->node_vals[TN_NODE_CUMULATIVE_DIFFUSIVE_DEPOSITION] + diff_dep;
    }
    inline void update_all_cumul_sed_depositions(const Eigen::VectorXd & sed_dep){
        this->node_vals[TN_NODE_CUMULATIVE_SEDIMENTATION_DEPOSITION] = this->node_vals[TN_NODE_CUMULATIVE_SEDIMENTATION_DEPOSITION] + sed_dep;
    }
    inline void update_all_cumul_imp_depositions(const Eigen::VectorXd & imp_dep){
        this->node_vals[TN_NODE_CUMULATIVE_IMPACTION_DEPOSITION] = this->node_vals[TN_NODE_CUMULATIVE_IMPACTION_DEPOSITION] + imp_dep;
    }
	inline void update_all_cumul_alv_depositions(const Eigen::VectorXd & alv_dep){
        this->node_vals[TN_NODE_CUMULATIVE_ALVEOLAR_DEPOSITION] = this->node_vals[TN_NODE_CUMULATIVE_ALVEOLAR_DEPOSITION] + alv_dep;
    }
    inline void update_all_relative_cumul_depositions(const Eigen::VectorXd & rel_dep){
        this->node_vals[TN_NODE_TOTAL_DEP_OVER_VOL] = this->node_vals[TN_NODE_TOTAL_DEP_OVER_VOL] + rel_dep;
    }
    inline void update_all_branching_angles(const Eigen::VectorXd & brangle){
        this->node_vals[TN_NODE_BRANGLE] = brangle;
    }
    //inline void update_total_depositions(const Eigen::VectorXd & dep_sum){
        //this->node_vals[TN_NODE_TOTAL_DEPOSITION] = this->node_vals[TN_NODE_TOTAL_DEPOSITION] + dep_sum;
    //}
    inline void update_generational_depositions(const Eigen::VectorXd & dep_sum){
        this->node_vals[TN_NODE_GENERATIONAL_DEPOSITION] = this->node_vals[TN_NODE_GENERATIONAL_DEPOSITION] + dep_sum;
    }
    inline void update_all_spect_depositions(const Eigen::VectorXd & dep){
        this->node_vals[TN_NODE_SPECT_DEPOSITION] = this->node_vals[TN_NODE_SPECT_DEPOSITION] + dep;
    }
    inline void update_all_airway_depositions(const Eigen::VectorXd & dep){
        this->edge_vals[TN_EDGE_AIRWAY_DEPOSITION] = this->edge_vals[TN_EDGE_AIRWAY_DEPOSITION] + dep;
    }
    inline void update_all_airway_imp_depositions(const Eigen::VectorXd & dep){
        this->edge_vals[TN_EDGE_AIRWAY_IMP_DEPOSITION] = this->edge_vals[TN_EDGE_AIRWAY_IMP_DEPOSITION] + dep;
    }
    inline void update_all_airway_diff_depositions(const Eigen::VectorXd & dep){
        this->edge_vals[TN_EDGE_AIRWAY_DIFF_DEPOSITION] = this->edge_vals[TN_EDGE_AIRWAY_DIFF_DEPOSITION] + dep;
    }
    inline void update_all_airway_sed_depositions(const Eigen::VectorXd & dep){
        this->edge_vals[TN_EDGE_AIRWAY_SED_DEPOSITION] = this->edge_vals[TN_EDGE_AIRWAY_SED_DEPOSITION] + dep;
    }
    inline void update_all_airway_numbers(const Eigen::VectorXd & airway){
        this->edge_vals[TN_EDGE_AIRWAY_NO] = airway;
    }
    inline void update_all_airway_generations(const Eigen::VectorXd & gen){
        this->edge_vals[TN_EDGE_GENERATION_NO] = gen;
    }
    inline void update_all_horsfield_orders(const Eigen::VectorXd & ord){
        this->edge_vals[TN_EDGE_HORSFIELD_ORDER] = ord;
    }
    inline void update_all_max_fluxes(const Eigen::VectorXd & vel){
        this->edge_vals[TN_EDGE_MAX_FLUX] = this->edge_vals[TN_EDGE_MAX_FLUX].cwiseMax(vel);
    }
    inline void update_node_pendelluft(const Eigen::VectorXd & pendelluft){
        this->node_vals[TN_NODE_PENDELLUFT] = this->node_vals[TN_NODE_PENDELLUFT] + pendelluft;
    }
    inline void update_all_airway_lengths(const Eigen::VectorXd & lengths){
        this->edge_vals[TN_EDGE_AIRWAY_LENGTH] = lengths;
    }
	inline void solve_concentration_update(const double & dt){ this->transport_solver->update_concentration(dt); }

    inline void set_branching_angle(size_t j1, size_t j2, double angle)
    {
        //ensures index is consistently defined so that smallest index comes first

        this->branch_angles[std::minmax(j1,j2)] = angle;
		//std::cout << this->get_horsfield_order(std::min(j1, j2)) << ' ' << this->branch_angles[std::minmax(j1, j2)] << std::endl;
        //std::cout << j2 << '\n' << branch_angles[std::minmax(j1,j2)] << std::endl;
    }
    inline double get_branching_angle(size_t j1, size_t j2)
    {
        //ensures index is consistently defined so that smallest index comes first
        return this->branch_angles[std::minmax(j1,j2)];
    }

    inline size_t get_airway_from_edge(size_t j)
    {
        return this->edge_airway_no[j];
    }
    inline size_t get_edges_in_airway(size_t J, size_t count)
    {
        return this->airway_to_edge_map[J][count];
    }
    inline size_t count_airways(){return airway_count;}//total number of airways
    inline size_t count_airway_edges(size_t J){return no_airway_edges[J];}//return number of edges in airway J

	double get_tot_IG_volume(const bool & UseInnerVol);
	inline double get_acinus_volume(const size_t kt){ return this->acinus_vector[kt]->get_volume(); }
	double get_acinus_IG_volume(const size_t & kt, const bool & UseInnerVol);
	inline void get_acinus_pos(const size_t & kt, network::Position & pos)
	{
		this->acinus_vector[kt]->get_position(pos);
	}
	int print_transport_vtk(const std::string & filename, SimParameterList *p);
    void print_end_transport_csv(const std::string & filename, SimParameterList *p);
    void print_end_transport_headers(const std::string & filename);

    inline std::vector<size_t> get_acin_node_map(const size_t kt){return acin_node_map[kt]; }
};

class FirstOrderTransportSolver: public TransportSolver
{
private:
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicgstab_solver;
	std::vector<Eigen::Triplet<double>> A_triplet;
	std::chrono::time_point<std::chrono::system_clock> start, end;
    Eigen::VectorXd L_sed, L_imp, L_diff, lam_sed, lam_diff, lam_diff_ingham, L_diff_ingham, L_alv, lam_sed_wang, L_sed_wang;
public:
	FirstOrderTransportSolver():TransportSolver(){};
	FirstOrderTransportSolver(AirwayTransportNetwork *t, const bool & deposition, const bool & alveolar_dep):TransportSolver(t)
	{
		setup_transport_calc(deposition, alveolar_dep);
	}

	void setup_transport_calc(const bool & deposition, const bool & alv_deposition);
	virtual void update_concentration(const double & dt);

    void update_diff_deposition_rate(const double & beta, const double & mol_diff);
    void update_diff_deposition_rate_ingham(const double & mol_diff);
    void update_diff_deposition_rate_yu();
    void update_diff_deposition_rate_ingham91(const std::vector<double> & alv_frac);
    void update_sed_deposition_rate(const double & u_sed, const std::vector<double> & alv_frac);
    void update_sed_deposition_rate_pich(const double & u_sed, const std::vector<double> & alv_frac);
    void update_sed_deposition_rate_wang(const double & u_sed);
    void update_imp_deposition_rate_yeh(const double & eff_imp, const double & d_p, const double & rho_p, const double mu_air);
    void update_imp_deposition_rate_zhang();
	void update_alv_deposition_rate(const double & dt);

};

inline std::shared_ptr<TransportSolver> create_transport_solver(const char & c, AirwayTransportNetwork *t,
																const bool & deposition, const bool & alv_deposition)  //no options for this yet
{
	return (std::make_shared<FirstOrderTransportSolver>(t,deposition,alv_deposition));
}

#endif
