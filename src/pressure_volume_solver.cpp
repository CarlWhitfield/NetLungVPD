#if defined(_WIN32) || defined(_WIN64)
	#ifndef _USE_MATH_DEFINES
		#define _USE_MATH_DEFINES
	#endif
#endif

#include "airway_flow_network.h"
#include "airway_edges_nodes.h"
#include "acinus_models.h"
#include <vector>
#include <Eigen/Sparse>
#include <math.h>
#include <iostream>
#include <chrono>
#include "globals.h"

PressureVolumeSolver::PressureVolumeSolver(AirwayFlowNetwork *t, const double & dtime)
{
	this->tree = t;
	this->dt = dtime;
}

void PressureVolumeSolver::update_flow(const Eigen::VectorXd & V, const Eigen::VectorXd & Vold)
{
	Eigen::VectorXd End_flux = (V - Vold)/this->dt;
	for(size_t order = 0; order < this->tree->count_horsfield_orders(); order++)
	{
		for(size_t index = 0; index < this->tree->count_edges_in_horsfield_order(order); index++)
		{
			size_t j = this->tree->get_edge_index_from_horsfield_order(order, index);
			size_t k_out = this->tree->get_node_out_index(j);
			double tot_flux_out;
			if(this->tree->count_edges_out(k_out) == 0)
			{
				tot_flux_out = End_flux[k_out - this->tree->get_first_term_index()];
			}
			else
			{
				tot_flux_out = 0;
				for (size_t jo = 0; jo < this->tree->count_edges_out(k_out); jo++)
				{
					tot_flux_out += this->tree->get_edge_flux(this->tree->get_edge_out_index(k_out, jo));
				}
			}
			this->tree->update_edge_flux(j,tot_flux_out);
		}
	}
}

void LinearSinusoidalPressureVolumeSolver::setup_flow_calc(const bool & Flow_BC)
{
	//the sinusoidal solver assumes that the periodic steady state pressures and volumes are of the form
	//p = p0 + ps*sin(2*pi*t/tau) + pc*sin(2*pi*t/tau) and
	//v = v0 + vs*sin(2*pi*t/tau) + vc*sin(2*pi*t/tau)
	//so we are solving for p0 and v0, ps and vs, and pc and vc separately

	this->tree->update_all_edge_resistances_by_rule();
	size_t Nnodes = this->tree->count_nodes();
	size_t Nedges = this->tree->count_edges();
	size_t Ntermnodes = this->tree->count_term_nodes();
	size_t term_start = this->tree->get_first_term_index();

	//solve 0 part first
	Eigen::VectorXd Kinv(Ntermnodes), PG(Ntermnodes), PGoverK(Ntermnodes), Vrest(Ntermnodes);
	double Vrestsum = 0, PGdotKinv = 0, Kinvsum = 0;
	for(size_t k = term_start; k < Nnodes; k++) //loop over terminal nodes
	{
		size_t k_rel = k - term_start;
		
		Vrest[k_rel] = this->tree->get_termnode_intrinsic_volume(k_rel);  //sum up 0 tension volume (all same if no grav)
		Kinv[k_rel] = 1.0 / this->tree->get_termnode_elastance(k_rel);   //1 / elastance
		PG[k_rel] = this->tree->get_termnode_relative_pressure(k_rel);   //change in Ppl due to grav (0 if no gravity)
		PGoverK[k_rel] = PG[k_rel] * Kinv[k_rel]; 

		PGdotKinv += PGoverK[k_rel];
		Vrestsum += Vrest[k_rel];
		Kinvsum += Kinv[k_rel];
	}
	//so the pleural pressure (which can vary in space, but assumed separable is)
	//Pleural pressure: Ppl = Ppl0 + PG + PPls*sin(2*pi*t/tau) + Pplc*cos(2*pi*t/tau)
	//Total acinar volume: Vtot = VFRC + VT/2 + Vs*sin(2*pi*t/tau) + Vc*cos(2*pi*t/tau)
	//Assume acinar tension is proportional to deviation from FRC
	//Intra acinus pressure: P_k = Ppl_k + K_k (V_k - Vrest_k) + R_k dV_k/dt + I_k d2V_k/dt2
	//the constant part of the pressure P_k must be same as mouth (P_k = 0)
	//because the average flux must be 0, so by definition 
	this->P0 = Eigen::VectorXd::Zero(Nnodes);

	//define matrix to solve sine and cosine parts
	size_t Nvars = 2 * (Nnodes + Ntermnodes + 2);
	Eigen::SparseMatrix<double> Asin(Nvars, Nvars);

	//define vector for RHS of equation
	Eigen::VectorXd Bsin = Eigen::VectorXd::Zero(Nvars);

	//fill matrix with triplets (faster)
	std::vector<Eigen::Triplet<double>> Asin_triplet;
	Asin_triplet.reserve(4 * Nedges + 2 * Nnodes + 10 * Ntermnodes + 6); 
	//not sure if this is exactly the number of non zeros but it only needs to be close

	/*  
	We have the following equations:
	The first 2*Nnodes rows of Asin are the Conductance Laplacian for the pressure drop: 
	Lcond Psin = Qsin and Lcond Pcos = Qcos
	If flow BC: Vtot = V0 + VT(1 - cos(2 pi t / tau))/2
	[Qsin, Qcos] = [(VT pi / tau), 0] on entry node
	If pressure BC
	[Qsin, Qcos] = to be solved for on entry node
	For both cases:
	[Qsin, Qcos] = [0, 0] on internal nodes
	[Qsin, Qcos] = [-(2 Vc pi / tau), (2 Vs pi / tau)] on terminal nodes

	Then there are 2*Nterm rows of Asin for the terminal nodes 
	Psin = PPlsin + K Vsin - R (2 pi / tau) Vcos - I (2 pi / tau)^2 Vsin  on terminal nodes
	Pcos = PPlcos + K Vcos + R (2 pi / tau) Vsin - I (2 pi / tau)^2 Vcos  on terminal nodes
	If flow BC 
	[Pplsin, PPlcos] = to be solved for
	If pressure BC
	[Pplsin, PPlcos] = [0,DP/2]

	For both cases:
	[Psin, Pcos] = 0 on entry node

	So define a common solution vector
	[Psin, Pcos, Vsin, Vcos, PPlsin, PPlcos, Q0sin, Q0cos]
	which has length 2*(Nnodes + Nterm + 2)
	*/

	//make note of starting points of different parts of the solution vector
	size_t Pcos_start = Nnodes;
	size_t Vsin_start = 2 * Nnodes;
	size_t Vcos_start = 2 * Nnodes + Ntermnodes;
	size_t Pplsin_pos = 2 * (Nnodes + Ntermnodes);
	size_t Pplcos_pos = 2 * (Nnodes + Ntermnodes) + 1;
	size_t Q0sin_pos = 2 * (Nnodes + Ntermnodes + 1);
	size_t Q0cos_pos = 2 * (Nnodes + Ntermnodes + 1) + 1;

	//first fill the laplacian in the first two diagonal blocks of Asin
	for(size_t k = 0; k < Nnodes; k++)
	{
		double w_degree = 0;
		size_t n_in = this->tree->count_edges_in(k);
		size_t n_out = this->tree->count_edges_out(k);
		//loop over all connected nodes -> fill off diagonal laplacian entries
		for(size_t counter = 0; counter < n_in + n_out; counter++)
		{
			double res;
			size_t ki;
			if(counter < n_in)  //is an edge in
			{
				res = this->tree->get_edge_resistance(this->tree->get_edge_in_index(k, counter));
				if(this->tree->get_edge_in_index(k, counter) == this->tree->get_edge_out_index(0,0))
				{
					res += this->tree->get_mouth_resistance();   //add mouth resistance to first edge
				}
				ki = this->tree->get_node_in_index(this->tree->get_edge_in_index(k, counter));
			} 
			else     //is an edge out
			{
				res = this->tree->get_edge_resistance(this->tree->get_edge_out_index(k,counter-n_in));
				if(this->tree->get_edge_out_index(k,counter-n_in) == this->tree->get_edge_out_index(0,0))
				{
					res += this->tree->get_mouth_resistance();  //add mouth resistance to first edge
				}
				ki =  this->tree->get_node_out_index(this->tree->get_edge_out_index(k,counter-n_in));
			}

			double weight = 1.0 / res; 
			w_degree += weight;
			Asin_triplet.push_back(Eigen::Triplet<double>(((int) k), ((int) ki), -weight));
			Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Pcos_start + k)), ((int)(Pcos_start + ki)), -weight));
		}
		//diagonal term
		Asin_triplet.push_back(Eigen::Triplet<double>(((int) k), ((int) k), w_degree));
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Pcos_start + k)), ((int)(Pcos_start + k)), w_degree));
	}
	//now add the terms for Qtermsin and Qcossin
	double tfactor = 2 * M_PI / (this->tau);
	for (size_t k = term_start; k < Nnodes; k++) //loop over terminal nodes
	{
		size_t k_rel = k - term_start;
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)k), ((int)(Vcos_start + k_rel)), -tfactor));
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Pcos_start + k)), ((int)(Vsin_start + k_rel)), tfactor));
	}
	//now add the terms for Q0sin and Q0cos
	Asin_triplet.push_back(Eigen::Triplet<double>(((int) 0), ((int) Q0sin_pos), -1.0));
	Asin_triplet.push_back(Eigen::Triplet<double>(((int) Nnodes), ((int) Q0cos_pos), -1.0));

	//now the first 2*Nnodes rows of Asin are filled

	//now fill the entries for the next 2*Nterm rows
	for(size_t k = term_start; k < Nnodes; k++) //loop over terminal nodes
	{
		size_t k_rel = k - term_start;
		AirwayNetworkNode *node = this->tree->get_node(k);
		double Kh = this->tree->get_termnode_elastance(k_rel);
		double Ih = this->tree->get_termnode_inertance(k_rel);
		double Rbh = this->tree->get_termnode_resistance(k_rel);

		//Identity matrices (P terms)
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Vsin_start + k_rel)), ((int) k), 1.0));
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Vcos_start + k_rel)), ((int)(Pcos_start + k)), 1.0));
		//Elastance and Inertance matrix (V terms)
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Vsin_start + k_rel)), ((int)(Vsin_start + k_rel)), -Kh + tfactor*tfactor*Ih));
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Vcos_start + k_rel)), ((int)(Vcos_start + k_rel)), -Kh + tfactor*tfactor*Ih));
		//Resistance matrix (Vterms)
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Vsin_start + k_rel)), ((int)(Vcos_start + k_rel)), -tfactor * Rbh));
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Vcos_start + k_rel)), ((int)(Vsin_start + k_rel)), tfactor * Rbh));
		//Pleural pressure term
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Vsin_start + k_rel)), ((int)(Pplsin_pos)), -1.0));
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Vcos_start + k_rel)), ((int)(Pplcos_pos)), -1.0));
	}

	//set boundary conditions in the last 4 rows
	//pressure at inlet = 0 
	Asin_triplet.push_back(Eigen::Triplet<double>(((int) (Pplsin_pos)), 0, 1.0));
	Asin_triplet.push_back(Eigen::Triplet<double>(((int) (Pplcos_pos)), ((int)Pcos_start), 1.0));
	//RHS is already 0 so no need to modify B
	if (Flow_BC){
		//set Q0sin and Q0cos
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Q0sin_pos)), ((int)(Q0sin_pos)), 1.0));
		Bsin[Q0sin_pos] = 0.5*tfactor*(this->VT);
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Q0cos_pos)), ((int)(Q0cos_pos)), 1.0));
		//cos part is zero so no need to modify B
	}
	else {   //pleural pressure BC
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Q0sin_pos)), ((int)(Pplsin_pos)), 1.0));
		//sin part is zero so no need to modify B
		Asin_triplet.push_back(Eigen::Triplet<double>(((int)(Q0cos_pos)), ((int)(Pplcos_pos)), 1.0));
		Bsin[Q0cos_pos] = 0.5*this->DP;
	}

	Asin.setFromTriplets(Asin_triplet.begin(), Asin_triplet.end());
	Eigen::VectorXd Xsin = Eigen::VectorXd::Zero(Nvars);
	Eigen::SparseLU<Eigen::SparseMatrix<double>> LU_solver;
	std::cout << "Solving flow problem...\n";
	auto start = std::chrono::system_clock::now();
	LU_solver.compute(Asin);
	Xsin = LU_solver.solve(Bsin);
	if(LU_solver.info() != Eigen::Success)
	{
		std::cout << "Solver_unsuccessful.\n";
		std::cout << LU_solver.lastErrorMessage() << '\n';
		abort_on_failure();
	}
	auto end = std::chrono::system_clock::now();
	std::cout << "Flow problem solution took " << (std::chrono::duration<double>(end-start)).count() <<  "\n";

	this->Psin = Xsin.segment(0,Nnodes);
	this->Pcos = Xsin.segment(Pcos_start,Nnodes);
	this->Vsin = Xsin.segment(Vsin_start,Ntermnodes);
	this->Vcos = Xsin.segment(Vcos_start, Ntermnodes);
	this->Ppsin = Xsin[Pplsin_pos];
	this->Ppcos = Xsin[Pplcos_pos];
	this->Q0sin = Xsin[Q0sin_pos];
	this->Q0cos = Xsin[Q0cos_pos];
	//if we have a pressure bc, the simulation starts at max PPl, not start inhalation
	//so we should transform to start inhalation at t=0
	if(!Flow_BC) {
		//shift to sine flow with t-t0 transformation
		double Q0s_new = sqrt(this->Q0sin*this->Q0sin + this->Q0cos*this->Q0cos);
		double Cos_t0 = this->Q0sin / Q0s_new;   //cos(2*pi*t0/tau)
		double Sin_t0 = this->Q0cos / Q0s_new;   //sin(2*pi*t0/tau)
		//new sine parts
		Eigen::VectorXd Ps_new = this->Psin*Cos_t0 + this->Pcos*Sin_t0;
		Eigen::VectorXd Vs_new = this->Vsin*Cos_t0 + this->Vcos*Sin_t0;
		double Pps_new = this->Ppsin*Cos_t0 + this->Ppcos*Sin_t0;
		//new cosine parts
		Eigen::VectorXd Pc_new = this->Pcos*Cos_t0 - this->Psin*Sin_t0;
		Eigen::VectorXd Vc_new = this->Vcos*Cos_t0 - this->Vsin*Sin_t0;
		double Ppc_new = this->Ppcos*Cos_t0 - this->Ppsin*Sin_t0;

		//assign the transformed variables
		this->Psin = Ps_new;
		this->Pcos = Pc_new;
		this->Vsin = Vs_new;
		this->Vcos = Vc_new;
		this->Ppsin = Pps_new;
		this->Ppcos = Ppc_new;
		this->Q0sin = Q0s_new;
		this->Q0cos = 0;
	}

	if (Flow_BC) {
		this->Pp0 = (Vrestsum - this->VFRC - 0.5*this->VT - PGdotKinv) / Kinvsum;
	}
	else {
		this->Pp0 = (Vrestsum - this->VFRC + this->Vcos.sum() - PGdotKinv) / Kinvsum;
	}
	//the vector of V0 values is then
	this->V0 = Vrest - this->Pp0 * Kinv - PGoverK;

	this->tree->update_all_termnode_volumes(this->V0 + Vcos);
	this->tree->update_all_node_pressures(this->P0 + Pcos);
	this->tree->update_Ppl(Pp0 + Ppcos);
	this->Vold = V0 + Vcos*cos(-2 * M_PI * this->dt / tau) + Vsin*sin(-2 * M_PI * this->dt / tau);
	this->tree->update_all_termnode_intrinsic_volumes(this->Vold);
	this->update_flow(this->tree->get_all_termnode_volumes(), Vold);
	//std::cout << "Ppl cos: " << this->Ppcos << ". Ppl sin: " << this->Ppsin << std::endl;
	//std::cout << "V cos: " << Vsin.sum() << ". V sin: " << Vcos.sum() << std::endl;
	//std::cout << "Q tot " << sqrt(Q0sin*Q0sin + Q0cos*Q0cos) <<	". V tot" 
	//		  << sqrt(Vsin.sum()*Vsin.sum() + Vcos.sum()*Vcos.sum()) << std::endl;

}

void LinearSinusoidalPressureVolumeSolver::update_press_and_vol(const size_t & itime)
{
	double sint = sin(2 * M_PI * (itime*this->dt) / tau);
	double cost = cos(2 * M_PI * (itime*this->dt) / tau);

	this->tree->update_all_node_pressures(P0 + Psin * sint  + Pcos * cost);
	this->Vold = V0 + Vsin * sin(2 * M_PI * (int(itime)-1) * this->dt / tau) 
		                      + Vcos * cos(2 * M_PI * (int(itime)-1) * this->dt / tau);
	this->tree->update_all_termnode_volumes(V0 + Vsin * sint + Vcos * cost);
	this->tree->update_Ppl(Pp0 + Ppsin * sint + Ppcos * cost);

	this->update_flow(this->tree->get_all_termnode_volumes(), this->Vold);
}

double LinearSinusoidalPressureVolumeSolver::get_expansion_volume(const size_t & itime)
{
	double sint = sin(2 * M_PI * (itime*this->dt) / tau);
	double cost = cos(2 * M_PI * (itime*this->dt) / tau);

	return (0.5*this->VT + (this->Q0cos*sint - this->Q0sin*cost)*this->tau / (2 * M_PI));
}

void LinearGenericPressureVolumeSolver::initialise(Eigen::SparseMatrix<double> & mat)
{
	size_t Nnodes = this->tree->count_nodes();
	size_t Nedges = this->tree->count_edges();
	size_t Ntermnodes = this->tree->count_term_nodes();
	size_t term_start = this->tree->get_first_term_index();

	//define matrix
	this->mat_dim_size =  Nedges + Nnodes + Ntermnodes + 1;

	mat = Eigen::SparseMatrix<double>(this->mat_dim_size, this->mat_dim_size);
	this->x = Eigen::VectorXd::Zero(this->mat_dim_size);
	
	//fill x vector
	this->x.head(Nedges) = this->tree->get_all_edge_fluxes();
	this->x.segment(Nedges,Nnodes) = this->tree->get_all_node_pressures();
	this->x.segment(Nedges+Nnodes,Ntermnodes) = this->tree->get_all_termnode_volumes();
	this->x[x.size()-1] = this->tree->get_Ppl();
	//fill Jacobian for time 0
	this->fill_matrix(mat);

	//get FRC volume
	this->Vold = this->tree->get_all_termnode_volumes();
	this->Vold_old = this->Vold;
	this->V0 = this->tree->sum_all_termnode_volumes();
	
	//run several beraths to intialise
	if(this->flow_opt != FLOW_FROM_FILE)
	{
		this->parse_flow_input_option();
	}
	else
	{
		this->lung_volume_model = std::make_shared<SinFlowObject>(this->dt, this->tau);
	}
}

void LinearGenericPressureVolumeSolver::fill_matrix(Eigen::SparseMatrix<double> & mat)
{
	size_t Nnodes = this->tree->count_nodes();
	size_t Nedges = this->tree->count_edges();
	size_t Ntermnodes = this->tree->count_term_nodes();
	size_t term_start = this->tree->get_first_term_index();

	std::vector<Eigen::Triplet<double>> jac_fill;
	jac_fill.reserve(3*(Nedges + Nnodes + Ntermnodes));

	//pressure drop on all edges, Jacobian of Pi - Pj = rij(qij) qij
	for(size_t j=0; j<this->tree->count_edges(); j++)
	{
		size_t i_nodein = this->tree->get_node_in_index(j);
		size_t i_nodeout = this->tree->get_node_out_index(j);
		const network::TubeGeometry *geom = this->tree->get_edge(j)->get_geom();
		jac_fill.push_back(Eigen::Triplet<double>(int(j), int(j), this->tree->get_edge_resistance_for_flux(j, x[j])
			                                                    + this->tree->get_edge_resistance_grad_for_flux(j, x[j]) * x[j]));
		jac_fill.push_back(Eigen::Triplet<double>(int(j), int(Nedges + i_nodein), -1.0));
		jac_fill.push_back(Eigen::Triplet<double>(int(j), int(Nedges + i_nodeout), 1.0));
	}

	//sum of fluxes into - out of nodes = 0
	for(size_t k=0; k<Nnodes; k++)
	{
		for (size_t jin = 0; jin < this->tree->count_edges_in(k); jin++)
		{
			jac_fill.push_back(Eigen::Triplet<double>(int(Nedges+k), int(this->tree->get_edge_in_index(k,jin)), 1.0));
		}

		for (size_t jout = 0; jout < this->tree->count_edges_out(k); jout++)
		{
			jac_fill.push_back(Eigen::Triplet<double>(int(Nedges+k), int(this->tree->get_edge_out_index(k,jout)), -1.0));
		}

		//if term node, flux must equal volume change
		if ( k >= term_start)
		{
			size_t kt = k - term_start;
			jac_fill.push_back(Eigen::Triplet<double>(int(Nedges+k), int(Nedges + Nnodes + kt), -1.0 / this->dt));
		}
	}

	//constitutive equation for term nodes Pk = Ppl + K(Vk) Vk + R dVkdt + I d2Vkdt2 
	for(size_t k=term_start; k<Nnodes; k++)
	{
		size_t kt = k - term_start;
		size_t kcoeff = Nnodes + Nedges + kt;
		jac_fill.push_back(Eigen::Triplet<double>(int(kcoeff), int(kcoeff), this->tree->get_termnode_elastance_for_vol(kt,x[kcoeff])
			                                      + this->tree->get_termnode_elastance_vol_grad(kt,x[kcoeff]) * x[kcoeff] 
		                                          + this->tree->get_termnode_resistance(kt) / this->dt 
												  + this->tree->get_termnode_inertance(kt) / (( this->dt )*( this->dt ))));
		jac_fill.push_back(Eigen::Triplet<double>(int(kcoeff), int(Nedges + Nnodes + Ntermnodes), 1.0));
		jac_fill.push_back(Eigen::Triplet<double>(int(kcoeff), int(Nedges + k), -1.0));
	}

	//P0 = 0
	jac_fill.push_back(Eigen::Triplet<double>(int(Nedges + Nnodes + Ntermnodes), int(Nedges), 1.0));

	mat.setFromTriplets(jac_fill.begin(), jac_fill.end());
}

void LinearGenericPressureVolumeSolver::update_from_solution()
{
	size_t Nnodes = this->tree->count_nodes();
	size_t Nedges = this->tree->count_edges();
	size_t Ntermnodes = this->tree->count_term_nodes();

	this->tree->update_all_node_pressures( this->x.segment(Nedges,Nnodes) );
	this->tree->update_all_termnode_volumes( this->x.segment(Nedges+Nnodes,Ntermnodes) );
	this->update_flow(this->tree->get_all_termnode_volumes(), this->Vold);
	this->tree->update_Ppl( this->x[Nnodes + Nedges + Ntermnodes] );
}

void LinearGenericPressureVolumeSolver::run_breaths(const int & nbreaths)
{
	for(size_t itime = 0; itime <= ((size_t) (nbreaths*(this->tau) / (this->dt))); itime++)
	{
		this->update_press_and_vol(itime);
	}
}

void LinearGenericPressureVolumeSolver::setup_flow_calc(const bool & Flow_BC)
{
	this->initialise(this->A);
	//if linear, can factorise first
	this->LU_solver.compute(this->A);
	this->run_breaths(5);
	this->parse_flow_input_option();
}

void LinearGenericPressureVolumeSolver::update_flow_equations(Eigen::VectorXd & func, const size_t & itime)
{
	size_t Nnodes = this->tree->count_nodes();
	size_t Nedges = this->tree->count_edges();
	size_t Ntermnodes = this->tree->count_term_nodes();
	size_t term_start = this->tree->get_first_term_index();

	func = Eigen::VectorXd::Zero(this->mat_dim_size);

	double Vnext = this->V0 + this->VT*( this->lung_volume_model->get_volume(itime) );
	//get last volume, could work out from same formula but with (t-dt), however lung inflation function may have changed
	double Vlast = this->tree->sum_all_termnode_volumes();
	//fill f vector
	//Pi - Pj = rij qij
	for(size_t j=0; j<Nedges; j++)
	{
		size_t i_nodein = this->tree->get_node_in_index(j);
		size_t i_nodeout = this->tree->get_node_out_index(j);
		func[j] = (this->tree->get_edge_resistance_for_flux(j, x[j])*x[j])  
			    - (this->x[Nedges + i_nodein] - this->x[Nedges + i_nodeout]);
	}
	//flux at all internal nodes = 0
	for(size_t k=0; k<Nnodes; k++)
	{
		for (size_t jin = 0; jin < this->tree->count_edges_in(k); jin++)
		{
			func[Nedges+k] += this->x[this->tree->get_edge_in_index(k,jin)];
		}

		for (size_t jout = 0; jout < this->tree->count_edges_out(k); jout++)
		{
			func[Nedges+k] -= this->x[this->tree->get_edge_out_index(k,jout)];
		}

		//if entry node, flux is set by bc
		if ( k==0 ) 
		{
			func[Nedges] += ( Vnext - Vlast ) / this->dt;
		}

		//if term node, flux must equal volume change
		if ( k >= this->tree->get_first_term_index())
		{
			size_t kt = k - this->tree->get_first_term_index();
			func[Nedges + k] -= (this->x[Nedges + Nnodes + kt] - this->Vold[kt]) / this->dt;
		}
	}
	for(size_t k=term_start; k<Nnodes; k++)
	{
		size_t kt = k - term_start;
		size_t kcoeff = Nnodes + Nedges + kt;
		//std::cout << kt << ": " << x[Nnodes + Nedges + Ntermnodes] << ' ' << x[Nedges+k] << ' ' 
		//	      << this->tree->get_termnode_relative_pressure(kt) << ' ' 
		//		  << this->tree->get_termnode_elastance_for_vol(kt, x[kcoeff]) * x[kcoeff] << '\n';
		func[Nnodes + Nedges + kt] = x[Nnodes + Nedges + Ntermnodes] - x[Nedges+k]  
		     + this->tree->get_termnode_relative_pressure(kt)
			 + this->tree->get_termnode_elastance_for_vol(kt, x[kcoeff]) * x[kcoeff]
	         + this->tree->get_termnode_resistance(kt)*(x[kcoeff] - this->Vold[kt]) / this->dt
			 + this->tree->get_termnode_inertance(kt)*(x[kcoeff] - 2*this->Vold[kt] + this->Vold_old[kt]) 
			 / ((this->dt)*(this->dt));
	}
	func[Nnodes + Nedges + Ntermnodes] = x[Nedges];
}

void LinearGenericPressureVolumeSolver::update_press_and_vol(const size_t & itime)
{
	if(this->flow_opt == FLOW_FROM_FILE || this->breath_sols.size() == 0)
	{
		Eigen::VectorXd dx = Eigen::VectorXd::Zero(this->mat_dim_size);
		Eigen::VectorXd func = Eigen::VectorXd::Zero(this->mat_dim_size);

		this->Vold_old = this->Vold;
		this->Vold = this->tree->get_all_termnode_volumes();
		this->update_flow_equations(func, itime);
		dx = this->LU_solver.solve(-func);
		if(this->LU_solver.info() != Eigen::Success)
		{
			std::cout << "Solver_unsuccessful.\n";
			std::cout << LU_solver.lastErrorMessage() << '\n';
			abort_on_failure();
		}
		//update
		this->x = this->x + dx;
	}
	else
	{
		size_t nbreath = size_t(itime*this->dt/tau);
		size_t i_t = itime - nbreath*int_round<size_t,double>(tau/dt);
		this->x = breath_sols[i_t];
	}
	this->update_from_solution();
}

double LinearGenericPressureVolumeSolver::get_expansion_volume(const size_t & itime)
{
	return (this->VT*this->lung_volume_model->get_volume(itime));
}

void LinearGenericPressureVolumeSolver::parse_flow_input_option()
{
	switch(this->flow_opt)
	{
	case SIN_FLOW:
		{
			this->lung_volume_model = std::make_shared<SinFlowObject>(this->dt, this->tau);
		}	break;

	case STEP_FLOW:
		{
			this->lung_volume_model = std::make_shared<StepFlowObject>(this->dt, this->tau);
		}	break;

	case FLOW_FROM_FILE:
		{
			this->lung_volume_model = std::make_shared<FileFlowObject>(this->flow_in, this->dt);
		}	break;

	default:
		{
			std::cout << "Flow option not recognised.\n";
			abort_on_failure();
		}
	}

}

void LinearGenericPressureVolumeSolver::run_breaths_until_converged()
{
	double diff = 1;
	size_t Nt = int_round<size_t,double>(this->tau/this->dt);
	std::vector<Eigen::VectorXd> breath_sols_old;
	breath_sols_old.resize(Nt);
	this->breath_sols.resize(Nt);
	for(size_t n = 0; n < Nt; n++)
	{
		this->breath_sols[n] = Eigen::VectorXd::Zero(this->mat_dim_size);
	}

	size_t iter = 0;
	size_t Nnodes = this->tree->count_nodes();
	size_t Nedges = this->tree->count_edges();
	size_t Ntermnodes = this->tree->count_term_nodes();
	int Nout = Nt / 20;
	std::cout << "Running breath convergence." << std::endl;
	while(diff > 1E-06)
	{
		std::cout << "Breath number = " << iter << std::endl;
		breath_sols_old = this->breath_sols;
		diff = 0;
		for(size_t n = 0; n < Nt; n++)
		{
			this->update_press_and_vol(iter*Nt + n);
			if ((n+1)%Nout == 0 || (n+1) == Nt) std::cout << "Breath " << iter << ": " << 100.0*(n+1) / Nt << "% complete." << std::endl;
			this->breath_sols[n] = this->x;
			double norm1 = this->x.head(Nedges).norm();
			diff += (this->x.head(Nedges) - breath_sols_old[n].head(Nedges)).norm() / norm1;
			double norm2 = this->x.segment(Nedges,Nnodes).norm();
			diff += (this->x.segment(Nedges,Nnodes) - breath_sols_old[n].segment(Nedges,Nnodes)).norm() / norm2;
			double norm3 = this->x.segment(Nedges+Nnodes, Ntermnodes).norm();
			diff += (this->x.segment(Nedges+Nnodes, Ntermnodes) - breath_sols_old[n].segment(Nedges+Nnodes, Ntermnodes)).norm() / norm3;
			diff += fabs(this->x[Nnodes + Nedges + Ntermnodes] - breath_sols_old[n][Nnodes + Nedges + Ntermnodes]) / fabs(this->x[Nnodes + Nedges + Ntermnodes]);
		}
		diff /= 4*Nt;
		iter++;
	}
	this->breaths_converged = true;
}

void GenericPressureVolumeSolver::setup_flow_calc(const bool & Flow_BC)
{
	this->initialise(this->Jacobian);
	LU_solver.analyzePattern(this->Jacobian);
	if(this->flow_opt == FLOW_FROM_FILE)
	{
		this->run_breaths(5);
	}
	else
	{
		this->run_breaths_until_converged();
	}
	this->parse_flow_input_option();
}

void GenericPressureVolumeSolver::update_press_and_vol(const size_t & itime)
{
	if(this->flow_opt == FLOW_FROM_FILE || !this->breaths_converged)
	{
		Eigen::VectorXd dx = Eigen::VectorXd::Zero(this->mat_dim_size);  //x step
		Eigen::VectorXd func = Eigen::VectorXd::Zero(this->mat_dim_size);   //to be solved
		Eigen::VectorXd Vones = Eigen::VectorXd::Ones(this->mat_dim_size);  //vector for operations
		Eigen::SparseMatrix<double> Jacobian_old;
		int iter = 0;
		this->Vold_old = this->Vold;
		this->Vold = this->tree->get_all_termnode_volumes();
		this->update_flow_equations(func, itime);   //get function to be minimised
		double res = 0;
		for(size_t irow = 0; irow < this->mat_dim_size; irow++)
		{
			res += func[irow]*func[irow];    //calculate residuals
		}
		res = sqrt(res);
		while (res > this->tol)
		{
			//solve
			//iterative solver is slower so use LU
			
			this->update_jacobian();     //get the jacobian
			double dJ = 0;    //normlaised change in Jacobian
			if (iter > 0)
			{
				Eigen::SparseMatrix<double> Jdiff = this->Jacobian - Jacobian_old;  //get Jacobian change
				Eigen::VectorXd RowsumJold = Jacobian_old.cwiseAbs()*Vones;   //sum the rows of Jold
				for (size_t k = 0; k < this->mat_dim_size; k++)
				{
					if (RowsumJold[k] == 0.0) RowsumJold[k] = 1.0;   //if rowsum is zero, set to 1
				}
				dJ = (RowsumJold.cwiseInverse().transpose()*Jdiff).norm();  //use Row sum as normalisation
			}
			if ((iter % 10) == 0 || dJ > 1E-02)
			{
				//if first iteration or Jacobian has changed by more than 1% or 10 iterations have passed
				Jacobian_old = this->Jacobian;   //update last Jacobian used
				this->LU_solver.factorize(this->Jacobian);  //update Jacobian factorisation
			}
			dx = this->LU_solver.solve(-func);

			/*this->bicgstab_solver.compute(this->Jacobian);
			Eigen::VectorXd dx_old = dx;
			dx = this->bicgstab_solver.solveWithGuess(-func, dx_old);*/

			if(this->LU_solver.info() != Eigen::Success)
			{
				std::cout << "LU Solver unsuccessful.\n";
				output_Eigen_error_message(this->LU_solver.info());
				abort_on_failure();
			}
			//}

			//update
			this->x = this->x + dx;
			this->update_flow_equations(func, itime);
			res = 0;
			for(size_t irow = 0; irow < this->mat_dim_size; irow++)
			{
				res += func[irow]*func[irow];
			}
			res = sqrt(res);
			iter++;
		}
	}
	else
	{
		size_t nbreath = size_t(itime*this->dt/tau);
		size_t i_t = itime - nbreath*int_round<size_t,double>(tau/dt);
		this->x = breath_sols[i_t];
	}
	this->update_from_solution();
	this->tree->update_all_edge_resistances_by_rule();
}

//attempted to do gradient descent, might need to return and consider normalising x for it to make sense
//or use standard algorithm
//void GenericPressureVolumeSolver::update_press_and_vol(const size_t & itime)
//{
//	if(this->flow_opt == FLOW_FROM_FILE || !this->breaths_converged)
//	{
//		Eigen::VectorXd dx = Eigen::VectorXd::Zero(this->mat_dim_size);  //x step
//		Eigen::VectorXd func = Eigen::VectorXd::Zero(this->mat_dim_size);   //to be solved
//		Eigen::VectorXd func_old = func;
//		int iter = 0;
//		this->Vold_old = this->Vold;
//		this->Vold = this->tree->get_all_termnode_volumes();
//		this->update_flow_equations(func, itime);   //get function to be minimised
//		double res = 0;
//		for(size_t irow = 0; irow < this->mat_dim_size; irow++)
//		{
//			res += func[irow]*func[irow];    //calculate residuals
//		}
//		res = sqrt(res);
//		std::cout << itime << ' ' << res << ' ' << this->tol << std::endl;
//		while (res > this->tol)
//		{
//			//solve
//			//iterative solver is slower
//			std::cout << itime << ' ' << iter << std::endl;
//			//if (iter == 0)
//			//{
//			//	this->update_jacobian();     //get the jacobian
//			//	this->LU_solver.factorize(this->Jacobian);
//			//	dx = this->LU_solver.solve(-func);
//			//	if (this->LU_solver.info() != Eigen::Success)
//			//	{
//			//		std::cout << "LU Solver unsuccessful.\n";
//			//		output_Eigen_error_message(this->LU_solver.info());
//			//		abort_on_failure();
//			//	}
//			//	std::cout << "F = " << res << std::endl;
//			//	std::cout << "dx = " << dx << std::endl;
//			//}
//			//else
//			//{
//				Eigen::VectorXd dx_old = dx;
//
//				Eigen::SparseMatrix<double> jacobian_old = this->Jacobian;
//				Eigen::VectorXd DelFold = jacobian_old.transpose()*func_old / func_old.norm();
//				this->update_jacobian();     //get the jacobian
//				std::cout << "Jacobian = " << this->Jacobian << std::endl;
//				Eigen::VectorXd DelF = this->Jacobian.transpose()*func / func.norm();
//				std::cout << "DelF = " << DelF << std::endl;
//				Eigen::VectorXd dDelF = DelF - DelFold;
//				double alpha = abs(dx_old.dot(dDelF)) / dDelF.squaredNorm();
//				dx = -alpha*DelF;
//				std::cout << "f = " << res << std::endl;
//				std::cout << "dx = " << dx << std::endl;
//			//}
//
//			/*this->bicgstab_solver.compute(this->Jacobian);
//			Eigen::VectorXd dx_old = dx;
//			dx = this->bicgstab_solver.solveWithGuess(-func, dx_old);*/
//
//			
//			//}
//
//			//update
//			this->x = this->x + dx;
//			func_old = func;
//			this->update_flow_equations(func, itime);
//			res = 0;
//			for(size_t irow = 0; irow < this->mat_dim_size; irow++)
//			{
//				res += func[irow]*func[irow];
//			}
//			res = sqrt(res);
//			iter++;
//		}
//	}
//	else
//	{
//		size_t nbreath = size_t(itime*this->dt/tau);
//		size_t i_t = itime - nbreath*int_round<size_t,double>(tau/dt);
//		this->x = breath_sols[i_t];
//	}
//	this->update_from_solution();
//	this->tree->update_all_edge_resistances_by_rule();
//}

void GenericPressureVolumeSolver::update_jacobian()
{
	for(size_t j=0; j<this->tree->count_edges(); j++)
	{
		const network::TubeGeometry *geom = this->tree->get_edge(j)->get_geom();
		this->Jacobian.coeffRef(j,j) = this->tree->get_edge_resistance_for_flux(j, x[j])
			                          + this->tree->get_edge_resistance_grad_for_flux(j,x[j]) * x[j];
	}

	for(size_t k=this->tree->get_first_term_index(); k<this->tree->count_nodes(); k++)
	{
		size_t kt = k - this->tree->get_first_term_index();
		size_t kcoeff = this->tree->count_nodes() + this->tree->count_edges() + kt;
		this->Jacobian.coeffRef(kcoeff, kcoeff) = this->tree->get_termnode_elastance_for_vol(kt, x[kcoeff])
			                                      + this->tree->get_termnode_elastance_vol_grad(kt, x[kcoeff]) * x[kcoeff]
		                                          + this->tree->get_termnode_resistance(kt) / this->dt 
												  + this->tree->get_termnode_inertance(kt) / (( this->dt )*( this->dt ));
	}
}

void FileFlowObject::parse_flow_input(const std::vector<double> & flow, const double & infl0)
{
	this->inflation.resize(flow.size() + 1);
	this->inflation[0] = infl0;

	for(size_t n = 0; n < flow.size(); n++)
	{
		this->inflation[n+1] = this->inflation[n] + (this->dt) * (flow[n]);
	}
}

double FileFlowObject::get_volume(const size_t & i_time)
{
	if (i_time < this->inflation.size()) return (this->inflation[i_time]);
	else 
	{
		std::cout << "Reached end of flow input.\n";
		abort_on_failure();
		return 0;
	}
}

double SinFlowObject::get_volume(const size_t & i_time)
{
	return ( 0.5*(1 - cos(2*M_PI*i_time*(this->dt) / (this->tau))) );
}

double StepFlowObject::get_volume(const size_t & i_time)
{
	int bn = int(2*i_time*(this->dt) / this->tau);
	double lin_func = 2*i_time*(this->dt)/this->tau - ((double) bn);
	if(bn%2)  //odd:exhale
	{
		return (1 - lin_func);
	}
	else    //even:inhale
	{
		return (lin_func);
	}
}


