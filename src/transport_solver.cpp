#include "airway_transport_network.h"
#include <vector>
#include <fstream>
#include <iostream>

TransportSolver::TransportSolver(AirwayTransportNetwork *t)
{
	this->tree = t;//comment 2
	this->fv_vol_old = this->tree->get_all_node_volumes();
	this->inner_fv_vol_old = this->tree->get_all_inner_node_volumes();
}

void FirstOrderTransportSolver::update_diff_deposition_rate(const double &beta, const double &mol_diff)
{
    double D_B, C_c, d_p, lambda;
    d_p = this->particle_size;
    lambda = this->mean_free_path;
    C_c = 1 + (2*lambda/d_p) * (1.257 + 0.4*std::exp(-0.55 * d_p / lambda));
    D_B = this->kB_T_mu * C_c / (3 * M_PI * d_p);
    this->lam_diff = Eigen::VectorXd::Zero(this->tree->count_edges());
    for (int j = 0; j < this->tree->count_edges(); ++j) {
        if (this->tree->get_edge(j)->get_geom()->get_inner_radius() > 0) {
            double a1 = 2.4048;//first zero of the 0th bessel function of 1st kind
            lam_diff[j] = D_B * this->tree->get_edge(j)->get_inner_volume()
                          / (std::pow(this->tree->get_edge(j)->get_geom()->get_inner_radius() / a1,2));

        }
    }
    Eigen::SparseMatrix<double> Inc;
    this->tree->get_incidence_matrix(Inc);
    Eigen::SparseMatrix<double> Inc_u = Inc.cwiseAbs();

    this->L_diff = 0.5 * (Inc_u.transpose() * lam_diff);

}

void FirstOrderTransportSolver::update_sed_deposition_rate(const double & u_sed, const std::vector<double> & alv_frac)
{
    this->lam_sed = Eigen::VectorXd::Zero(this->tree->count_edges());
    for (int j = 0; j < this->tree->count_edges(); ++j) {
        if (this->tree->get_edge(j)->get_geom()->get_inner_radius() > 0) {
            int Nbranches = this->tree->get_edge(j)->branch_count();
            lam_sed[j] = u_sed * sin(this->tree->get_edge_gravity_angle(j)) * this->tree->get_edge(j)->get_inner_volume()
                    / (this->tree->get_edge(j)->get_geom()->get_inner_radius());
            for (int i = 1; i < 9; ++i) {
                if (Nbranches == pow(2,i)){
                    lam_sed[j] *= 1 - alv_frac[i-1];
                }
            }
        }
    }
    Eigen::SparseMatrix<double> Inc;
    this->tree->get_incidence_matrix(Inc);
    Eigen::SparseMatrix<double> Inc_u = Inc.cwiseAbs();

    this->L_sed = 0.5 * (Inc_u.transpose() * lam_sed);
}


void FirstOrderTransportSolver::update_sed_deposition_rate_pich(const double & u_sed, const std::vector<double> & alv_frac)
//Update sedimentation deposition term using formula from Pich (1972)
{
    this->lam_sed = Eigen::VectorXd::Zero(this->tree->count_edges());
    Eigen::VectorXd edge_airway_length = Eigen::VectorXd::Zero(this->tree->count_edges());
    edge_airway_length = this->tree->get_all_airway_lengths();
    for (int j = 0; j < this->tree->count_edges(); ++j) {
        double Nbranches = this->tree->get_edge(j)->branch_count();
        double R = this->tree->get_edge(j)->get_geom()->get_inner_radius();
        double U = 0;
        if (R > 0) {U = std::abs(this->tree->get_edge_flux(j)) / (M_PI * Nbranches * std::pow(R,2));}
        if (R * U > 0){
            //define deposition efficiency term, eff
            double eps = 3 * u_sed * edge_airway_length(j) * sin(this->tree->get_edge_gravity_angle(j)) /
                    (8 * R * U);
            eps = std::min(eps,1.0);
            double eff;
            eff = (2/M_PI) * (2 * eps * sqrt(1 - pow(eps,2.0/3.0)) - pow(eps,1.0/3.0)* sqrt(1 - pow(eps,2.0/3.0)) +
                    asin(pow(eps,1.0/3.0)));
            double Nels = edge_airway_length(j) / this->tree->get_edge(j)->get_geom()->get_length();
            lam_sed[j] = (eff/Nels) * std::abs(this->tree->get_edge_flux(j));
        }
    }

    Eigen::SparseMatrix<double> Inc;
    this->tree->get_incidence_matrix(Inc);
    //Define unsigned incidence matrix
    Eigen::SparseMatrix<double> Inc_u = Inc.cwiseAbs();

    //Define a modified incidence matrix
    Eigen::VectorXd Q = this->tree->get_all_edge_fluxes();//vector of all edge fluxes
    Eigen::VectorXd sgnQ = Q.cwiseSign();
    Eigen::SparseMatrix<double> Inc_v(this->tree->count_edges(), this->tree->count_nodes());
    Inc_v = sgnQ.asDiagonal() * Inc;
    Inc_v = (Inc_v + Inc_v.cwiseAbs()) / 2;//set negative entries to zero

    //Assign half the loss term from each connecting edge to the node
    this->L_sed = 0.5 * (Inc_u.transpose() * lam_sed);
    //Instead, can assign the whole loss term from the upstream connected edge(s), and none from the downstream edge(s)
    //this->L_diff = Inc_v.transpose() * lam_diff;
}

void FirstOrderTransportSolver::update_imp_deposition_rate_yeh(const double & eff_imp, const double & d_p, const double & rho_p, const double mu_air)
//Update impaction loss term based on formula from Yeh & Schum (1980)
{
    this->L_imp = Eigen::VectorXd::Zero(this->tree->count_nodes());
    //Define vector of Stokes numbers
    Eigen::VectorXd Stk_vector = Eigen::VectorXd::Zero(this->tree->count_edges());
    double lambda = this->mean_free_path;
    double C_c = 1 + (2*lambda/d_p) * (1.257 + 0.4*std::exp(-0.55 * d_p / lambda));
    for(int j = 0; j < this->tree->count_edges(); ++j) {
        int NbranchesStk = this->tree->get_edge(j)->branch_count();
        Stk_vector[j] = std::abs(this->tree->get_edge_flux(j)) * rho_p * std::pow(d_p,2) * C_c
                / (NbranchesStk * M_PI * std::pow(this->tree->get_edge(j)->get_geom()->get_inner_radius(),3) * 36 * mu_air);
    }

    for(size_t k = 0; k < this->tree->count_nodes(); k++)
    {
        //sort edges into whether they flow into or out of the given node

        //define vectors of edge indices of edges flowing into/out of node k
        std::vector<size_t> flux_into_node_edge_index, flux_outof_node_edge_index;
        for (size_t ji = 0; ji < this->tree->count_edges_in(k); ji++)
        {
            size_t j = this->tree->get_edge_in_index(k,ji);
            //positive flux means flowing into node
            if(this->tree->get_edge_flux(j) > 0) flux_into_node_edge_index.push_back(j);
            else  flux_outof_node_edge_index.push_back(j);
        }
        for (size_t jo = 0; jo < this->tree->count_edges_out(k); jo++)
        {
            size_t j = this->tree->get_edge_out_index(k,jo);
            //positive flux means flowing out of node
            if(this->tree->get_edge_flux(j) > 0) flux_outof_node_edge_index.push_back(j);
            else  flux_into_node_edge_index.push_back(j);
        }

        for (int j_in = 0; j_in < flux_into_node_edge_index.size(); ++j_in) {
            for (int j_out = 0; j_out < flux_outof_node_edge_index.size(); ++j_out) {
                //define a deposition efficiency, eff, associated with each pair of edges connected to node k
                double eff = 0;
                double b_angle = this->tree->get_branching_angle(flux_into_node_edge_index[j_in],
                                                                 flux_outof_node_edge_index[j_out]);
                double Stk = Stk_vector[flux_into_node_edge_index[j_in]];

                //maximum no. of branches of the two edges
                int Nbranches = std::max(this->tree->get_edge(flux_into_node_edge_index[j_in])->branch_count(),
                                         this->tree->get_edge(flux_outof_node_edge_index[j_out])->branch_count());
                if (Nbranches > 1)
                {
                    //If the number of branches in one of the edges is larger than 1, i.e. if the node is part of an
                    // acinus, then assign the efficiency to be its mean given the branching angle is a uniformly
                    // distributed random variable in [0,pi/2];
                    double u = M_PI * Stk / 2.0;
                    if (u > 1)
                    {
                        eff = 1 - 4 / (3 * M_PI * u);
                    }
                        else
                    {
                        eff = 1 - (2 / M_PI) * acos(u) + (2 / (3*M_PI*u)) * (3 * std::sqrt(1 - std::pow(u,2)) - 2.0 - std::pow(1 - std::pow(u,2),1.5));
                    }
                }
                    else
                {
                    //for conducting airways, use original formula from Yeh & Schum (1980)
                    if (b_angle * Stk < 1)
                    {
                        eff = 1 - (2 / M_PI) * acos(b_angle * Stk) + (1 / M_PI) * sin(2 * acos(b_angle * Stk));
                    }
                    else
                    {
                    eff = 1.0;}
                }
                eff = std::min(eff,1.0);
                eff = std::max(eff,0.0);

                //assign deposition term
                if(this->imp_both_dirs_option) {
                    this->L_imp[k] += eff *
                                      std::min(std::abs(this->tree->get_edge_flux(flux_into_node_edge_index[j_in])),
                                               std::abs(this->tree->get_edge_flux(flux_outof_node_edge_index[j_out])));
                }
                else
                {
                    //in this case, only assign deposition if flow is towards distal end of the airways
                    if (this->tree->get_edge_flux(flux_into_node_edge_index[j_in]) > 0 &&
                        this->tree->get_edge_flux(flux_outof_node_edge_index[j_out]) > 0)
                    {
                        this->L_imp[k] += eff *
                                          std::min(std::abs(this->tree->get_edge_flux(flux_into_node_edge_index[j_in])),
                                                   std::abs(this->tree->get_edge_flux(flux_outof_node_edge_index[j_out])));
                    }
                }
            }
        }
    }
}

void FirstOrderTransportSolver::update_imp_deposition_rate_zhang()
//Update impaction loss term based on formula from Zhang et al. (1997)
{
    double d_p, rho_p, mu_air, nu_air, C_c, lambda;
    d_p = this->particle_size;
    rho_p = this->particle_density;
    mu_air = this->air_viscosity;
    nu_air = mu_air / this->air_density;
    lambda = this->mean_free_path;
    C_c = 1 + (2*lambda/d_p) * (1.257 + 0.4*std::exp(-0.55 * d_p / lambda));//Cunningham slip correction factor

    this->L_imp = Eigen::VectorXd::Zero(this->tree->count_nodes());

    //vectors of Stokes and Reynolds numbers
    Eigen::VectorXd Stk_vector = Eigen::VectorXd::Zero(this->tree->count_edges());
    Eigen::VectorXd Re_vector = Eigen::VectorXd::Zero(this->tree->count_edges());
    for(int j = 0; j < this->tree->count_edges(); ++j) {
        int Nbranches = this->tree->get_edge(j)->branch_count();
        double R, U;
        R = this->tree->get_edge(j)->get_geom()->get_inner_radius();
        U = 0;
        if (R>0){U = std::abs(this->tree->get_edge_flux(j)) / (M_PI * Nbranches * std::pow(R,2));}
        Stk_vector[j] = std::abs(this->tree->get_edge_flux(j)) * rho_p * std::pow(d_p,2) * C_c
                        / (Nbranches * M_PI * std::pow(R,3) * 36 * mu_air);
        Re_vector[j] = 2.0 * U * R / nu_air;
    }

    for(size_t k = 0; k < this->tree->count_nodes(); k++)
    {
        //sort edges into whether they flow into or out of node

        //vectors of edge indices of edges flowing into/out of node k
        std::vector<size_t> flux_into_node_edge_index, flux_outof_node_edge_index;
        for (size_t ji = 0; ji < this->tree->count_edges_in(k); ji++)
        {
            size_t j = this->tree->get_edge_in_index(k,ji);
            //positive flux means flowing into node
            if(this->tree->get_edge_flux(j) > 0) flux_into_node_edge_index.push_back(j);
            else  flux_outof_node_edge_index.push_back(j);
        }
        for (size_t jo = 0; jo < this->tree->count_edges_out(k); jo++)
        {
            size_t j = this->tree->get_edge_out_index(k,jo);
            //positive flux means flowing out of node
            if(this->tree->get_edge_flux(j) > 0) flux_outof_node_edge_index.push_back(j);
            else  flux_into_node_edge_index.push_back(j);
        }


        for (int j_in = 0; j_in < flux_into_node_edge_index.size(); ++j_in) {
            for (int j_out = 0; j_out < flux_outof_node_edge_index.size(); ++j_out) {
                //define deposition efficiency, eff, for each pair of edges connected to node k
                double eff = 0;
                double b_angle = this->tree->get_branching_angle(flux_into_node_edge_index[j_in],
                                                                 flux_outof_node_edge_index[j_out]);
                double Stk = Stk_vector[flux_into_node_edge_index[j_in]];
                double Re = Re_vector[flux_into_node_edge_index[j_in]];

                //maximum no. of branches of the two edges
                int Nbranches = std::max(this->tree->get_edge(flux_into_node_edge_index[j_in])->branch_count(),
                                         this->tree->get_edge(flux_outof_node_edge_index[j_out])->branch_count());

                if (Nbranches > 1)
                {
                    //If the number of branches in one of the edges is larger than 1, i.e. if the node is part of a
                    // trumpet, then assign this fixed branching angle
                    b_angle = asin(2.0/M_PI);
                }

                //Empirical parameters from Zhang et al. (1997)
                double a, b, c, d, f3;
                if (Stk < 0.04){
                    a = 0.0;
                    b = 0.000654;
                    c = 55.7;
                    d = 0.954;
                }
                else{
                    a = 0.19;
                    b = -0.193;
                    c = -9.5;
                    d = 1.565;
                }
                f3 = a + b * exp(c * pow(Stk,d));
                eff = pow(Re,1.0/3.0) * sin(b_angle) * f3;
                size_t Nbr_in = this->tree->get_edge(flux_into_node_edge_index[j_in])->branch_count();

                eff = std::min(eff,1.0);
                eff = std::max(eff,0.0);

                //assign deposition terms
                if(this->imp_both_dirs_option) {
                    this->L_imp[k] += eff *
                                      std::min(std::abs(this->tree->get_edge_flux(flux_into_node_edge_index[j_in])),
                                               std::abs(this->tree->get_edge_flux(flux_outof_node_edge_index[j_out])));
                }
                else
                {
                    //in this case, only assign deposition if flow is towards distal end of the airways
                    if (this->tree->get_edge_flux(flux_into_node_edge_index[j_in]) > 0 &&
                                this->tree->get_edge_flux(flux_outof_node_edge_index[j_out]) > 0)
                    {
                        this->L_imp[k] += eff *
                                          std::min(std::abs(this->tree->get_edge_flux(flux_into_node_edge_index[j_in])),
                                                   std::abs(this->tree->get_edge_flux(flux_outof_node_edge_index[j_out])));
                    }
                }
            }
        }
    }
}

void FirstOrderTransportSolver::update_alv_deposition_rate(const double & dt)
//Alveolar Loss Term
{
	this->L_alv = Eigen::VectorXd::Zero(this->tree->count_nodes());
	for(size_t kt = 0; kt < this->tree->count_acini(); kt++)
	{
		for(size_t ka = 0; ka < this->tree->get_acin_node_map(kt).size(); ka++)
		{
			size_t k = this->tree->get_acin_node_map(kt)[ka];
			double alv_vol_change_rate = (this->tree->get_node_fv_volume(k) - this->tree->get_node_inner_fv_volume(k) -
									this->fv_vol_old[k] + this->inner_fv_vol_old[k])/dt;
			if(alv_vol_change_rate > 0) //inhaling
			{
				this->L_alv[k] = alv_vol_change_rate;
			}
		}
	}
}


void FirstOrderTransportSolver::update_diff_deposition_rate_ingham(const double &mol_diff)
//Update diffusion deposition rate based on formula from Ingham (1975)
{
    this->lam_diff_ingham = Eigen::VectorXd::Zero(this->tree->count_edges());
    Eigen::VectorXd edge_airway_length = Eigen::VectorXd::Zero(this->tree->count_edges());
    edge_airway_length = this->tree->get_all_airway_lengths();

    double D_B, C_c, d_p, lambda;
    d_p = this->particle_size;
    lambda = this->mean_free_path;
    C_c = 1 + (2*lambda/d_p) * (1.257 + 0.4*std::exp(-0.55 * d_p / lambda));
    D_B = this->kB_T_mu * C_c / (3 * M_PI * d_p);

    for (int j = 0; j < this->tree->count_edges(); ++j) {
        double airway_length = edge_airway_length(j);

        if (this->tree->get_edge(j)->get_geom()->get_inner_radius() > 0) {
            //define deposition efficiency, eff, using formula from Ingham (1975)
            double Del = 0;
            if (std::abs(this->tree->get_edge_flux(j))>0) {
                int Nbranches = this->tree->get_edge(j)->branch_count();
                Del = M_PI * D_B * airway_length * Nbranches / (4 * std::abs(this->tree->get_edge_flux(j)));
            }
            double eff = 1.0 - 0.819 * std::exp(-14.63 * Del) - 0.0976 * std::exp(-89.22 * Del) - 0.0325 * std::exp(-228.0 * Del) - 0.0509 * std::exp(-125.9 * std::pow(Del,2.0/3.0));
            double Nels = airway_length / this->tree->get_edge(j)->get_geom()->get_length();
            eff = std::min(eff,1.0);
            lam_diff[j] = (eff/Nels) * std::abs(this->tree->get_edge_flux(j));
        }
    }

    this->lam_diff = Eigen::VectorXd::Zero(this->tree->count_edges());
    double a1 = 2.4048;//first zero of the 0th bessel function of 1st kind
    for (int j = 0; j < this->tree->count_edges(); ++j) {
        if (this->tree->get_edge(j)->get_geom()->get_inner_radius() > 0) {
            lam_diff[j] = D_B * this->tree->get_edge(j)->get_inner_volume()
                          / (std::pow(this->tree->get_edge(j)->get_geom()->get_inner_radius() / a1,2));
            if (lam_diff_ingham(j) < lam_diff(j))
            {
                //Uncomment this to add a simple diffusion term when it is larger than Ingham term.
                //lam_diff_ingham[j] = lam_diff(j);
            }
        }
    }

    //Define unsigned incidence matrix
    Eigen::SparseMatrix<double> Inc;
    this->tree->get_incidence_matrix(Inc);
    Eigen::SparseMatrix<double> Inc_u = Inc.cwiseAbs();

    //Define a modified incidence matrix
    Eigen::VectorXd Q = this->tree->get_all_edge_fluxes();//vector of all edge fluxes
    Eigen::VectorXd sgnQ = Q.cwiseSign();
    Eigen::SparseMatrix<double> Inc_v(this->tree->count_edges(), this->tree->count_nodes());
    Inc_v = sgnQ.asDiagonal() * Inc;
    Inc_v = (Inc_v + Inc_v.cwiseAbs()) / 2;//set negative entries to zero

    //Assign half the loss term from each connecting edge to the node
    this->L_diff_ingham = 0.5 * (Inc_u.transpose() * lam_diff_ingham);
    //Instead, can assign the whole loss term from the upstream connected edge(s), and none from the downstream edge(s)
    //this->L_diff = Inc_v.transpose() * lam_diff_ingham;
}


void FirstOrderTransportSolver::update_diff_deposition_rate_yu()
//update diffusion deposition rate based on formula from Yu & Cohen (1994)
{
    this->lam_diff = Eigen::VectorXd::Zero(this->tree->count_edges());
    Eigen::VectorXd edge_airway_length = Eigen::VectorXd::Zero(this->tree->count_edges());
    edge_airway_length = this->tree->get_all_airway_lengths();

    double D_B, d_p, nu, Sc, lambda, C_c;
    d_p = this->particle_size;
    lambda = this->mean_free_path;
    C_c = 1 + (2*lambda/d_p) * (1.257 + 0.4*std::exp(-0.55 * d_p / lambda));//Cunningham slip correction factor
    D_B = this->kB_T_mu * C_c / (3 * M_PI * d_p);
    nu = this->air_viscosity / this->air_density;
    Sc = nu / D_B;//Schmidt number

    //Empirical parameters from Yu & Cohen (1994)
    double a, b, c, d;
    a = 1.2027;
    b = -0.6067;
    c = -0.5108;
    d = 0.5081;

    for (int j = 0; j < this->tree->count_edges(); ++j) {
        int Nbranches = this->tree->get_edge(j)->branch_count();
        double R, U, Re, L, eff, Nels;
        R = this->tree->get_edge(j)->get_geom()->get_inner_radius();
        U = 0;

        if (R > 0) {U = std::abs(this->tree->get_edge_flux(j)) / (M_PI * Nbranches * std::pow(R,2));}
        Re = 2 * U * R / nu;
        L = edge_airway_length(j);
        Nels = L / this->tree->get_edge(j)->get_geom()->get_length();

        //Efficiency formula
        eff = 0.0;
        if (Re > 0) {
            eff = a * std::pow(Re, b) * std::pow(Sc, c) * std::pow(L / R, d);
        }
        eff = std::min(eff,1.0);
        lam_diff[j] = (eff/Nels) * std::abs(this->tree->get_edge_flux(j));
    }

    //Define unsigned incidence matrix
    Eigen::SparseMatrix<double> Inc;
    this->tree->get_incidence_matrix(Inc);
    Eigen::SparseMatrix<double> Inc_u = Inc.cwiseAbs();

    //Define modified incidence matrix
    Eigen::VectorXd Q = this->tree->get_all_edge_fluxes();//vector of all edge fluxes
    Eigen::VectorXd sgnQ = Q.cwiseSign();
    Eigen::SparseMatrix<double> Inc_v(this->tree->count_edges(), this->tree->count_nodes());
    Inc_v = sgnQ.asDiagonal() * Inc;
    Inc_v = (Inc_v + Inc_v.cwiseAbs()) / 2;//set negative entries to zero

    //Assign half the loss term from each connecting edge to the node
    this->L_diff = 0.5 * (Inc_u.transpose() * lam_diff);
    //Instead, can assign the whole loss term from the upstream connected edge(s), and none from the downstream edge(s)
    //this->L_diff = Inc_v.transpose() * lam_diff;

}


void FirstOrderTransportSolver::update_diff_deposition_rate_ingham91(const std::vector<double> & alv_frac)
//Update diffusion deposition rate based on formula from Ingham (1991)
{
    this->lam_diff = Eigen::VectorXd::Zero(this->tree->count_edges());
    Eigen::VectorXd edge_airway_length = Eigen::VectorXd::Zero(this->tree->count_edges());
    edge_airway_length = this->tree->get_all_airway_lengths();

    double D_B, d_p, nu, Sc, lambda, C_c;
    d_p = this->particle_size;
    lambda = this->mean_free_path;
    C_c = 1 + (2*lambda/d_p) * (1.257 + 0.4*std::exp(-0.55 * d_p / lambda));//Cunningham slip correction factor
    D_B = this->kB_T_mu * C_c / (3 * M_PI * d_p);
    nu = this->air_viscosity / this->air_density;
    Sc = nu / D_B;//Schmidt number

    //Parameters from Ingham (1991)
    double a, b, c, d;
    a = 3.033;
    b = -0.5556;
    c = -0.6667;
    d = 0.5556;

    for (int j = 0; j < this->tree->count_edges(); ++j) {
        int Nbranches = this->tree->get_edge(j)->branch_count();
        double R, U, Re, L, eff, Nels;
        R = this->tree->get_edge(j)->get_geom()->get_inner_radius();
        U = 0;

        if (R > 0) {U = std::abs(this->tree->get_edge_flux(j)) / (M_PI * Nbranches * std::pow(R,2));}
        Re = 2 * U * R / nu;
        L = edge_airway_length(j);
        Nels = L / this->tree->get_edge(j)->get_geom()->get_length();

        //Efficiency formula
        eff = 0.0;
        if (Re > 0) {
            eff = a * std::pow(Re, b) * std::pow(Sc, c) * std::pow(L / R, d);
        }
        eff = std::min(eff,1.0);
        lam_diff[j] = (eff/Nels) * std::abs(this->tree->get_edge_flux(j));

        //Adjust deposition rate based on alveolated fraction of acinar ducts, if using this adjustment
        if (Nbranches > 1) {
            for (int i = 1; i < 9; ++i) {
                if (Nbranches == pow(2, i)) {
                    lam_diff[j] *= 1 - alv_frac[i - 1];
                }
            }
        }
    }

    //Define unsigned incidence matrix
    Eigen::SparseMatrix<double> Inc;
    this->tree->get_incidence_matrix(Inc);
    Eigen::SparseMatrix<double> Inc_u = Inc.cwiseAbs();

    //Define modified incidence matrix
    Eigen::VectorXd Q = this->tree->get_all_edge_fluxes();//vector of all edge fluxes
    Eigen::VectorXd sgnQ = Q.cwiseSign();
    Eigen::SparseMatrix<double> Inc_v(this->tree->count_edges(), this->tree->count_nodes());
    Inc_v = sgnQ.asDiagonal() * Inc;
    Inc_v = (Inc_v + Inc_v.cwiseAbs()) / 2;//set negative entries to zero

    //Assign half the loss term from each connecting edge to the node
    this->L_diff = 0.5 * (Inc_u.transpose() * lam_diff);
    //Instead, can assign the whole loss term from the upstream connected edge(s), and none from the downstream edge(s)
    //this->L_diff = Inc_v.transpose() * lam_diff;
}

void FirstOrderTransportSolver::update_sed_deposition_rate_wang(const double & u_sed)
//update sedimentation deposition rate based on formula from Wang (1975)
{
    this->lam_sed = Eigen::VectorXd::Zero(this->tree->count_edges());
    Eigen::VectorXd edge_airway_length = Eigen::VectorXd::Zero(this->tree->count_edges());
    Eigen::VectorXd edge_airway_length2 = this->tree->get_all_airway_lengths();

    Eigen::VectorXd airway_length2 = Eigen::VectorXd::Zero(this->tree->count_airways());
    for (int J = 0; J < this->tree->count_airways(); ++J) {
        for (int i = 0; i < this->tree->count_airway_edges(J); ++i) {
            airway_length2[J] += this->tree->get_edge(this->tree->get_edges_in_airway(J,i))->get_geom()->get_length();
        }
    }
    //Eigen::VectorXd edge_airway_length = Eigen::VectorXd::Zero(this->tree->count_edges());
    for (int j = 0; j < this->tree->count_edges(); ++j) {
        size_t J = this->tree->get_airway_from_edge(j);
        //assign airway lengths
        edge_airway_length[j] = airway_length2[J];
    }

    for (int j = 0; j < this->tree->count_edges(); ++j) {
        double eff = 0.0;
        if (this->tree->get_edge(j)->get_geom()->get_inner_radius() > 0) {
            if (std::abs(this->tree->get_edge_flux(j))>0) {
                //Define deposition efficiency, eff, based on formula from Wang (1975)
                double R = this->tree->get_edge(j)->get_geom()->get_inner_radius();
                double U = std::abs(this->tree->get_edge_flux(j)) /
                        ( this->tree->get_edge(j)->branch_count() * M_PI * std::pow(R, 2));
                double L = edge_airway_length(j);
                double grangle = this->tree->get_edge_gravity_angle(j);
                double alpha = std::abs((M_PI / 2.0) - grangle);
                double eps = 3.0*u_sed*L* cos(alpha)/(8.0*U*R);
                if (grangle > (M_PI/2.0)){
                    //Uphill
                    double alph1, alph2;
                    alph1 = (M_PI/2.0) - (2.0 * R /(9.0*L)) * pow(2.0*u_sed / (3.0*U),0.5);
                    alph2 = (M_PI/2.0) - R * u_sed / (8.0 * U * L);
                    double sig = (u_sed * sin(alpha)/(6.0*U)) / (1 - (u_sed * sin(alpha)/(2.0 * U)));
                    double gam = pow(eps,2.0/3.0) / (1 - (u_sed * sin(alpha)/(2.0 * U)));
                    gam = std::min(gam,1.0);
                    sig = std::min(sig,1.0);
                    if (alpha < alph1){
                        eff = (1.0/M_PI) * (3.0* pow(sig*(1.0-sig),0.5) + asin(pow(1.0-sig,0.5)) +
                                (1.0-9.0* pow(sig,2.0))* asin(pow((1.0-sig)/(1+3.0*sig),0.5))) -
                                (2.0/M_PI) * (pow(gam*(1.0-gam),0.5) * (1.0 - 2.0*gam) + asin(pow(1.0-gam,0.5)));
                    } else if (alpha < alph2){
                        double zeta0;
                        zeta0 = R*u_sed* pow(sin(alpha),2.0) / (8.0*U*L* cos(alpha))
                                - pow(u_sed* sin(alpha)/(6.0*U),0.5)/16.0
                                + 7.0*L* cos(alpha)/(8.0*R* sin(alpha));
                        zeta0 = std::min(zeta0,1.0);
                        eff = (1.0/M_PI)*(1.0 + 3.0*sig)* sqrt(1.0-pow(zeta0,2.0))*
                                (zeta0 + (4.0* pow(gam,3.0/2.0)/ sqrt(1.0+3.0*sig)) -
                                sqrt(pow(zeta0,2.0) - 3.0*sig/(1.0+3.0*sig))) +
                                (1.0 - 9.0* pow(sig,2.0))* asin(sqrt(1.0 - pow(zeta0,2.0))) / M_PI -
                                (1.0/M_PI) * asin(sqrt((1.0 + 3.0*sig) * (1.0 - pow(zeta0,2.0))));
                    } else{eff = 0.0;}
                } else{
                    //Downhill
                    double zeta1;
                    zeta1 = pow(eps / (1 + (3.0*u_sed)/(4.0*u_sed*sin(alpha))),1.0/3.0);
                    zeta1 = std::min(zeta1,1.0);
                    eff = (2.0/M_PI) * asin(zeta1) + (sqrt(1- pow(zeta1,2.0)) / (M_PI * (1 + u_sed * sin(alpha)/U))) *
                                                             (4.0*eps - (2.0 + (u_sed* sin(alpha))/U) * zeta1);
                }
            }
        }
        eff = std::min(1.0,eff);
        eff = std::max(0.0,eff);
        double Nels = edge_airway_length(j) / this->tree->get_edge(j)->get_geom()->get_length();
        lam_sed[j] = (eff/Nels) * std::abs(this->tree->get_edge_flux(j));
    }

    //Define unsigned incidence matrix
    Eigen::SparseMatrix<double> Inc;
    this->tree->get_incidence_matrix(Inc);
    Eigen::SparseMatrix<double> Inc_u = Inc.cwiseAbs();

    //Assign half the loss term from each connecting edge to the node
    this->L_sed = 0.5 * (Inc_u.transpose() * this->lam_sed);

    //Modified incidence matrix
    Eigen::VectorXd Q = this->tree->get_all_edge_fluxes();//vector of all edge fluxes
    Eigen::VectorXd sgnQ = Q.cwiseSign();
    Eigen::SparseMatrix<double> Inc_v(this->tree->count_edges(), this->tree->count_nodes());
    Inc_v = sgnQ.asDiagonal() * Inc;
    Inc_v = (Inc_v + Inc_v.cwiseAbs()) / 2;//set negative entries to zero
    //Instead, can assign the whole loss term from the upstream connected edge(s), and none from the downstream edge(s)
    //this->L_sed = Inc_v.transpose() * this->lam_sed;
}




void FirstOrderTransportSolver::setup_transport_calc(const bool & deposition = false, 
													 const bool & alveolar_deposition = false) //just setup matrix shape
{

    //start = std::chrono::system_clock::now();
	size_t Nnodes = tree->count_nodes();
	
	//count entries
	size_t count = 0;
	for(size_t k = 0; k < Nnodes; k++)
	{
		size_t n_in = tree->count_edges_in(k);
		size_t n_out = tree->count_edges_out(k);
		count += n_in + n_out + 1;
	}
	//set up triplet vector for future use
	A_triplet.resize(count);  //same size as laplacian matrix
	count = 0;
	for(size_t k = 0; k < Nnodes; k++)
	{
		//loop over all connected nodes
		size_t n_in = tree->count_edges_in(k);
		size_t n_out = tree->count_edges_out(k);
		for(size_t counter = 0; counter <  n_in + n_out; counter++)
		{
			size_t ki;
			if(counter < n_in)  //get index of node in
			{
				ki = tree->get_node_in_index(tree->get_edge_in_index(k, counter));
			} 
			else     //get index of node out
			{
				ki = tree->get_node_out_index(tree->get_edge_out_index(k, counter-n_in));
			}
			A_triplet[count] = Eigen::Triplet<double>(((int) k), ((int) ki), 1.0);
			count++;
		}

		A_triplet[count] = Eigen::Triplet<double>(((int) k), ((int) k), 1.0);
		count++;
	}
	Eigen::SparseMatrix<double> A(Nnodes, Nnodes);
	A.setFromTriplets(A_triplet.begin(), A_triplet.end());
	bicgstab_solver.analyzePattern(A);
    //end = std::chrono::system_clock::now();

    if (deposition) {
        this->L_sed = Eigen::VectorXd::Zero(Nnodes); //initialisation
        this->L_imp = Eigen::VectorXd::Zero(Nnodes); // initialisation
        this->L_diff_ingham = Eigen::VectorXd::Zero(Nnodes); // initialisation
        this->L_diff = Eigen::VectorXd::Zero(Nnodes); // initialisation
		this->L_alv = Eigen::VectorXd::Zero(Nnodes); // initialisation
    }
    this->simulate_deposition = deposition;
	this->simulate_alveolar_deposition = alveolar_deposition;
}

void FirstOrderTransportSolver::update_concentration(const double & dt) {
//	start = std::chrono::system_clock::now();
    size_t Nnodes = this->tree->count_nodes();
    size_t Nedges = this->tree->count_edges();
    Eigen::SparseMatrix<double> A(Nnodes, Nnodes);//initialise A
    Eigen::VectorXd b = Eigen::VectorXd::Zero(Nnodes);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(Nnodes);
    DSVolume *mouth = this->tree->get_mouth_object();

    //Get matrices and vectors to contruct A in A*x = b
    Eigen::VectorXd D = this->tree->get_all_edge_diffusivities();//vector of all edge diffusivities
    Eigen::VectorXd Q = this->tree->get_all_edge_fluxes();//vector of all edge fluxes
    Eigen::VectorXd Vol = Eigen::VectorXd::Zero(this->tree->count_nodes());
    Eigen::VectorXd Vol_old = Eigen::VectorXd::Zero(this->tree->count_nodes());
    if (this->simulate_deposition) {
        Vol = this->tree->get_all_inner_node_volumes();
        Vol_old = this->inner_fv_vol_old;
    } else {
        Vol = this->tree->get_all_node_volumes();//vector of all edge fluxes
        Vol_old = this->fv_vol_old;
    }
    Eigen::VectorXd sgnQ = Q.cwiseSign();

    Eigen::SparseMatrix<double> Inc;
    this->tree->get_incidence_matrix(Inc);
    Eigen::SparseMatrix<double> Inc_v(Nedges, Nnodes);
    Inc_v = sgnQ.asDiagonal() * Inc;
    Inc_v = (Inc_v + Inc_v.cwiseAbs()) / 2;//set negative entries to zero

    A = dt * Inc.transpose() * (D.asDiagonal() * Inc + Q.asDiagonal() * Inc_v);//form matrix A
    A += Vol.asDiagonal();   //for some reason this is the only way to add to diagonal

    //Boundary Conditions
    if (this->tree->get_edge_flux(0) > 0) {
        //update mouth (inhalation)
        if (mouth != NULL) {
            std::vector<std::shared_ptr<FlexibleVolumeElement>> injected(1);
            injected[0] = std::make_shared<FlexibleVolumeElement>(this->tree->get_edge_flux(0) * dt,
                                                                  this->tree->mouth_bc(), mouth->get_conc_top());
            std::vector<std::shared_ptr<FlexibleVolumeElement>> ejected;
            mouth->shunt_down(injected, ejected);
			//get inhaled gas or particle volume
			double Vcin = 0;
			for (size_t n = 0; n < ejected.size(); n++)
			{
				Vcin += ejected[n]->get_igvol();
			}
            b[0] = Vcin/dt;   //set c[0] to be equal to bottom of mouth vol
            x[0] = Vcin/(dt*Q[0]);  //initial guess
        } else   //no mouth volume
        {
			b[0] = Q[0]*this->tree->mouth_bc();   //set c[0] to be equal to boundary condition
            x[0] = this->tree->mouth_bc();    //initial guess
        }
		A.coeffRef(0, 0) = Q[0] + D[0]; 
        A.coeffRef(0, this->tree->get_node_out_index(0)) = -D[0];
    } else {
		A.coeffRef(0, 0) = 1;
        A.coeffRef(0, this->tree->get_node_out_index(0)) = -1;   //no diffusive flux
        x[0] = this->tree->get_node_conc(tree->get_node_out_index(0));  //initial guess
    }
    //use last concentration as initial guess
    x.tail(this->tree->count_nodes() - 1) = this->tree->get_all_node_concs().tail(this->tree->count_nodes() - 1);
    //b is equal to v_old*c_old (element wise multiplication)
    b.tail(this->tree->count_nodes() - 1) = this->tree->get_all_node_concs().tail(this->tree->count_nodes() - 1).array()
                                            * Vol_old.tail(this->tree->count_nodes() - 1).array();

//    end = std::chrono::system_clock::now();
//    std::cout << "Transport Matrix fill took " << (std::chrono::duration<double>(end-start)).count() <<  "\n";
    //TODO: Put all deposition terms in solver child object and select object at initialisation...
    //based on input options

	//std::cout << A << '\n' << std::endl;
	//std::cout << D << '\n' << std::endl;
	//std::cout << Q << '\n' << std::endl;

    //Calculate the length of each airway
    Eigen::VectorXd airway_length = Eigen::VectorXd::Zero(this->tree->count_airways());
    for (int j = 0; j < this->tree->count_edges(); ++j) {
        size_t J = this->tree->get_airway_from_edge(j);
        airway_length[J] += this->tree->get_edge(j)->get_geom()->get_length();
    }
    Eigen::VectorXd edge_airway_length = Eigen::VectorXd::Zero(this->tree->count_edges());
    for (int j = 0; j < this->tree->count_edges(); ++j) {
        size_t J = this->tree->get_airway_from_edge(j);
        //assign airway lengths
        edge_airway_length[j] = airway_length[J];
    }
    this->tree->update_all_airway_lengths(edge_airway_length);

    Eigen::VectorXd airway_length2 = Eigen::VectorXd::Zero(this->tree->count_airways());
    for (int J = 0; J < this->tree->count_airways(); ++J) {
        for (int i = 0; i < this->tree->count_airway_edges(J); ++i) {
            airway_length2[J] += this->tree->get_edge(this->tree->get_edges_in_airway(J,i))->get_geom()->get_length();
        }
    }

    Eigen::VectorXd airway_length3 = this->tree->get_all_airway_lengths();
    if (this->simulate_deposition) {
        //Approximate alveolated fractions for acinar airways
        //Note that this is not currently being used in any of the deposition terms
        int Nb[8] = {2, 4, 8, 16, 32, 64, 128, 256};
        double wlength[8] = {1.33, 1.12, 0.93, 0.83, 0.7, 0.7, 0.7, 0.7}; //mm, will be rescaled
        double wrad[8] = {0.25, 0.245, 0.2, 0.19, 0.18, 0.17, 0.155, 0.145}; //mm, will be rescaled
        double walv_frac[8] = {0.4, 0.7, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        double duct_surf_areas[8];
        for (unsigned j = 1; j <= 8; j++) {
            //surface area of a single duct in each generation (mm^2)
            duct_surf_areas[j-1] = 2 * M_PI * wrad[j-1] * wlength[j-1];
        }
        double r_alv = 0.14;//mm (alveolar radius)
        double h_alv = 0.23;//mm (alveolar height)
        double duct_surf_ext[8];
        for (unsigned j = 1; j <= 8; j++) {
            //surface area of a single duct in each generation (mm^2)
            duct_surf_ext[j-1] = 2 * M_PI * (wrad[j-1]+h_alv-r_alv) * wlength[j-1];
        }
        double max_alv_area = M_PI * std::pow(r_alv,2);
        double N_alv[8];
        for (unsigned j = 1; j <= 8; j++) {
            //surface area of a single duct in each generation (mm^2)
            N_alv[j-1] = duct_surf_areas[j-1] / max_alv_area;
            N_alv[j-1] = std::floor(N_alv[j-1]);
        }
        double s_op = 0.561 * std::pow(2*r_alv,2);//alveolar opening area
        double alv_area[8];
        std::vector<double> alv_frac(8); //alveolated fraction for each generation of acinar airways
        for (unsigned j = 1; j <= 8; j++) {
            alv_area[j-1] = N_alv[j-1] * s_op;
            alv_frac[j-1] = 0;//Uncomment if not using alveolated fraction modification
            //Approximate alveolar surface fraction by ratio of s_op to max_alv_area
            //alv_frac[j-1] = walv_frac[j-1] * s_op / max_alv_area;
            //Alternative, approx maximum packing fraction
            //alv_frac[j-1] = std::min(walv_frac[j-1],0.85);
        }

        //Sedimentation
        double u_sed, d_p, g, C_c, mu_air, rho_p, lambda;
        mu_air = this->air_viscosity;//air viscosity
        d_p = this->particle_size;//particle diameter
        g = this->grav_accel;
        rho_p = this->particle_density;
        lambda = this->mean_free_path;
        C_c = 1 + (2 * lambda / d_p) * (1.257 + 0.4 * std::exp(-0.55 * d_p / lambda));
        u_sed = rho_p * pow(d_p, 2) * g * C_c / (18 * mu_air);
        pos::Position<double> grav_dir;
        grav_dir.x[0] = 0;
        grav_dir.x[1] = 0;
        grav_dir.x[2] = -1;
        this->tree->compute_all_edge_gravity_angles(grav_dir);

        //Update impaction deposition term
        switch (this->imp_formula_option) {
            case YEH_IMP_DEP: {
                double eff_imp = 1;
                this->update_imp_deposition_rate_yeh(eff_imp, d_p, rho_p, mu_air);
            }break;
            case ZHANG_IMP_DEP:
            {
                this->update_imp_deposition_rate_zhang();
            }break;
        }
        A += dt * this->L_imp.asDiagonal();

        //Update diffusion deposition term
        switch (this->diff_formula_option) {
            case INGHAM_DIFF_DEP:
            {
                this->update_diff_deposition_rate_ingham(this->tree->get_molecular_diffusivity());
            }break;
            case SIMPLE_DIFF_DEP:
            {
                this->update_diff_deposition_rate(0.0,this->tree->get_molecular_diffusivity());
            }break;
            case YU_DIFF_DEP:
            {
                this->update_diff_deposition_rate_yu();
            }break;
            case INGHAM91_DIFF_DEP:
            {
                this->update_diff_deposition_rate_ingham91(alv_frac);
            }break;
        }
        A += dt * this->L_diff.asDiagonal();

        //Update sedimentation deposition term
        switch (this->sed_formula_option) {
            case SIMPLE_SED_DEP:
            {
                this->update_sed_deposition_rate(u_sed, alv_frac);
            }break;
            case PICH_SED_DEP:
            {
                this->update_sed_deposition_rate_pich(u_sed, alv_frac);
            }break;
            case WANG_SED_DEP:
            {
                this->update_sed_deposition_rate_wang(u_sed);
            }break;
        }
        A += dt * this->L_sed.asDiagonal();


		//Alveolar deposition loss
		if(this->simulate_alveolar_deposition)
		{
			this->update_alv_deposition_rate(dt);
			A += dt * this->L_alv.asDiagonal();
		}
    }


    //Solve the linear problem that has just been set up
//	start = chrono::system_clock::now();
	bicgstab_solver.compute(A);
	Eigen::VectorXd xsol = bicgstab_solver.solveWithGuess(b,x);



//	end = chrono::system_clock::now();
//	std::cout << "Transport problem solution took " << (chrono::duration<double>(end-start)).count() <<  "\n";
	double trach_node_conc_old = this->tree->get_node_conc(0);
	this->update_vol_old(this->tree->get_all_node_volumes());    //store volume for next timestep
	this->update_inner_vol_old(this->tree->get_all_inner_node_volumes());    //store volume for next timestep
	this->tree->update_all_node_concs(xsol);
    //std::cout << xsol.maxCoeff() << std::endl;

    //Update things for outputs
    if (this->simulate_deposition) {
        //Node Depositions
        this->tree->update_all_node_depositions(
                (this->L_diff + this->L_sed + this->L_imp + this->L_alv).asDiagonal() * xsol);
        this->tree->update_all_cumul_depositions(
                dt * (this->L_diff + this->L_sed + this->L_imp + this->L_alv).asDiagonal() * xsol);
        this->tree->update_all_cumul_diff_depositions(dt * this->L_diff.asDiagonal() * xsol);
        this->tree->update_all_cumul_sed_depositions(dt * this->L_sed.asDiagonal() * xsol);
        this->tree->update_all_cumul_imp_depositions(dt * this->L_imp.asDiagonal() * xsol);
		this->tree->update_all_cumul_alv_depositions(dt * this->L_alv.asDiagonal() * xsol);

        //total depositions
        //Eigen::VectorXd dep_vector = dt * (this->L_diff + this->L_sed + this->L_imp + this->L_alv).asDiagonal() * xsol;
        Eigen::VectorXd dep_vector = dt * (this->L_diff + this->L_sed + this->L_imp + this->L_alv).asDiagonal() * xsol;
        Eigen::VectorXd dep_imp_vector = dt * (this->L_imp).asDiagonal() * xsol;
        //Eigen::VectorXd dep_diff_vector = dt * (this->L_diff).asDiagonal() * xsol;
        Eigen::VectorXd dep_diff_vector = dt * (this->L_diff).asDiagonal() * xsol;
        Eigen::VectorXd dep_grav_vector = dt * (this->L_sed).asDiagonal() * xsol;
		Eigen::VectorXd dep_alv_vector = dt * (this->L_alv).asDiagonal() * xsol;

        //Node deposition with acinar deposition assigned to the associated terminal conducting node (used to generate
        //gamma scintigraphy images)
        Eigen::VectorXd dep_spect = Eigen::VectorXd::Zero(this->tree->count_nodes());
        dep_spect[0] = dep_vector[0];
        for (int j = 0; j < this->tree->count_edges(); ++j) {
            size_t ki = this->tree->get_node_in_index(j);
            size_t ko = this->tree->get_node_out_index(j);
            size_t Nbr = this->tree->get_edge(j)->branch_count();
            if (Nbr == 1){dep_spect[ko] += dep_vector[ko];}//Conducting airway nodes same as original
            else {//Acinar edges - find terminal conducting node and add deposition to it
                double depo = dep_vector[ko];
                size_t kt;
                size_t jnew = j;
                size_t Nbr_new = Nbr;
                while (Nbr_new > 1){
                    kt = this->tree->get_node_in_index(jnew);
                    jnew = this->tree->get_edge_in_index(kt,0);
                    Nbr_new = this->tree->get_edge(jnew)->branch_count();
                }
                dep_spect[kt] += dep_vector[ko];
            }
        }
        this->tree->update_all_spect_depositions(dep_spect);

        //Deposition by generation
        Eigen::VectorXd len_gen = Eigen::VectorXd::Zero(35);//total length of airways in each generation
        Eigen::VectorXd surf_gen = Eigen::VectorXd::Zero(35);//total surface area of airways in each generation
        Eigen::VectorXd dep_gen = Eigen::VectorXd::Zero(35);//total deposition in each generation
        for (int k = 1; k < this->tree->count_nodes(); ++k) {
            size_t j = this->tree->get_edge_in_index(k, 0);
            double gen = this->tree->get_weibel_order(j);
            len_gen[gen] += this->tree->get_edge(j)->get_geom()->get_length();
            //n_gen[gen] += 1;
            surf_gen[gen] += this->tree->get_edge(j)->get_geom()->inner_area();
            dep_gen[gen] += dep_vector[k];
        }

        Eigen::VectorXd dep_gen_vector = Eigen::VectorXd::Zero(Nnodes);
        for (int k = 0; k < this->tree->count_nodes(); ++k) {
            size_t j;
            if (k == 0)
            {
                j = this->tree->get_edge_out_index(k,0);
            }
            else {
                j = this->tree->get_edge_in_index(k, 0);
            }
            double gen = this->tree->get_weibel_order(j);
            dep_gen_vector[k] = dep_gen[gen];
        }
        this->tree->update_generational_depositions(dep_gen_vector);

        //Save branching angles
        Eigen::VectorXd brangles = Eigen::VectorXd::Zero(this->tree->count_nodes());
        for (int i = 1; i < this->tree->count_term_nodes(); ++i) {
            brangles[i] = this->tree->get_branching_angle(this->tree->get_edge_in_index(i, 0),
                                                          this->tree->get_edge_out_index(i, 0));
        }
        this->tree->update_all_branching_angles(brangles);

        //Save deposition over edge volume
        Eigen::VectorXd Vol_edge = Eigen::VectorXd::Zero(this->tree->count_edges());
        for (int j = 0; j < this->tree->count_edges(); ++j) {
            Vol_edge[j] = this->tree->get_edge(j)->get_inner_volume();
        }
        Eigen::VectorXd Vol_node = Eigen::VectorXd::Zero(this->tree->count_nodes());
        Eigen::SparseMatrix<double> Inc_u = Inc.cwiseAbs();
        Vol_node = 0.5 * (Inc_u.transpose() * Vol_edge);
        Eigen::VectorXd rel_cum_dep =
                (dt * (this->L_diff + this->L_sed + this->L_imp + this->L_alv).asDiagonal() * xsol).array() /
                Vol_node.array();

        //.array() objects use element-wise operations
        this->tree->update_all_relative_cumul_depositions(rel_cum_dep);

        //Calculate total deposition in each airway (assigning bifurcation deposition to the airway above)
        Eigen::VectorXd airway_dep = Eigen::VectorXd::Zero(this->tree->count_airways());
        Eigen::VectorXd airway_imp_dep = Eigen::VectorXd::Zero(this->tree->count_airways());
        Eigen::VectorXd airway_diff_dep = Eigen::VectorXd::Zero(this->tree->count_airways());
        Eigen::VectorXd airway_grav_dep = Eigen::VectorXd::Zero(this->tree->count_airways());
		Eigen::VectorXd airway_alv_dep = Eigen::VectorXd::Zero(this->tree->count_airways());
        for (int k = 0; k < this->tree->count_nodes(); ++k) {
            size_t j;
            if (k == 0){j = this->tree->get_edge_out_index(k,0);}
            else{j = this->tree->get_edge_in_index(k,0);}
            size_t J = this->tree->get_airway_from_edge(j);
            airway_dep[J] += dep_vector[k];
            airway_imp_dep[J] += dep_imp_vector[k];
            airway_diff_dep[J] += dep_diff_vector[k];
            airway_grav_dep[J] += dep_grav_vector[k];
			airway_alv_dep[J] += dep_alv_vector[k];
        }

        //Deposition per whole airway
        Eigen::VectorXd edge_airway_dep = Eigen::VectorXd::Zero(this->tree->count_edges());
        Eigen::VectorXd edge_airway_imp_dep = Eigen::VectorXd::Zero(this->tree->count_edges());
        Eigen::VectorXd edge_airway_diff_dep = Eigen::VectorXd::Zero(this->tree->count_edges());
        Eigen::VectorXd edge_airway_sed_dep = Eigen::VectorXd::Zero(this->tree->count_edges());
		Eigen::VectorXd edge_airway_alv_dep = Eigen::VectorXd::Zero(this->tree->count_edges());
        for (int j = 0; j < this->tree->count_edges(); ++j) {
            size_t J = this->tree->get_airway_from_edge(j);
            edge_airway_dep[j] = airway_dep[J];
            edge_airway_imp_dep[j] = airway_imp_dep[J];
            edge_airway_diff_dep[j] = airway_diff_dep[J];
            edge_airway_sed_dep[j] = airway_grav_dep[J];
			edge_airway_alv_dep[j] = airway_alv_dep[J];
        }
        this->tree->update_all_airway_depositions(edge_airway_dep);
        this->tree->update_all_airway_imp_depositions(edge_airway_imp_dep);
        this->tree->update_all_airway_diff_depositions(edge_airway_diff_dep);
        this->tree->update_all_airway_sed_depositions(edge_airway_sed_dep);

        //Save airway generations and airway numbers
        Eigen::VectorXd edge_airway_no = Eigen::VectorXd::Zero(this->tree->count_edges());
        Eigen::VectorXd edge_generation_no = Eigen::VectorXd::Zero(this->tree->count_edges());
        Eigen::VectorXd edge_horder = Eigen::VectorXd::Zero(this->tree->count_edges());
        for (int j = 0; j < this->tree->count_edges(); ++j) {
            edge_airway_no[j] = this->tree->get_airway_from_edge(j);
            edge_generation_no[j] = this->tree->get_weibel_order(j);
            edge_horder[j] = this->tree->get_horsfield_order(j);
        }
        this->tree->update_all_airway_numbers(edge_airway_no);
        this->tree->update_all_airway_generations(edge_generation_no);
        this->tree->update_all_horsfield_orders(edge_horder);

        //Save maximum edge velocity
        Eigen::VectorXd edge_flux = Eigen::VectorXd::Zero(this->tree->count_edges());
        for (int j = 0; j < this->tree->count_edges(); ++j) {
            edge_flux[j] = this->tree->get_edge_flux(j);
        }
        this->tree->update_all_max_fluxes(edge_flux);

        //Pendelluft vector (not currently implemented)
        Eigen::VectorXd pendelluft = Eigen::VectorXd::Zero(this->tree->count_nodes());
        /*
        for (int k = 0; k < this->tree->count_nodes(); ++k) {
            for (size_t m = 0; m < this->tree->count_edges_out(k); m++) {
                for (int n = m + 1; n < this->tree->count_edges_out(k); ++n) {
                    double flux_m = this->tree->get_edge_flux(this->tree->get_edge_out_index(k, m));
                    double flux_n = this->tree->get_edge_flux(this->tree->get_edge_out_index(k, n));
                    if (flux_m * flux_n < 0.0) {
                        pendelluft(k) += dt * std::min(std::abs(flux_m), std::abs(flux_n));
                    }
                }
            }
        }
         */
        this->tree->update_node_pendelluft(pendelluft);
        Eigen::VectorXd pendelluft2 = this->tree->get_all_node_pendelluft();
        //std::cout << pendelluft.maxCoeff() << ' ' << pendelluft2.maxCoeff() << std::endl;

    }
	//after calc exhalation
	if(this->tree->get_edge_flux(0) < 0 && mouth != NULL)
	{
		std::vector<std::shared_ptr<FlexibleVolumeElement>> injected(1);
		injected[0] = std::make_shared<FlexibleVolumeElement>(-this->tree->get_edge_flux(0)*dt,
			                                     trach_node_conc_old, this->tree->get_node_conc(0));
		std::vector<std::shared_ptr<FlexibleVolumeElement>> ejected;
		mouth->shunt_up(injected, ejected);
	}
}

