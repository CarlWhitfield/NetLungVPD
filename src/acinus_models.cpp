#if defined(_WIN32) || defined(_WIN64)
	#ifndef _USE_MATH_DEFINES
		#define _USE_MATH_DEFINES
	#endif
#endif

#include "acinus_models.h"
#include "define_codes.h"
#include "globals.h"
#include <iostream>

using namespace network;

void WellMixedBagAcinus::build_model()
{
	std::vector<std::shared_ptr<AirwayNetworkNode>> nodes;
	nodes.push_back(node_in);
	std::vector<std::shared_ptr<AirwayNetworkEdge>> edges;

	double v3 = pow(this->vol / this->term_count(), 1.0/3.0);
	nodes[0]->set_intrinsic_vol(this->vol);

	this->tree = std::make_shared<FVElementNetwork>(nodes, edges);
}

void WellMixedBagAcinus::update_volume(const double & v)
{
//	double vscale = v/this->vol;
	this->vol = v;
    this->node_in->set_intrinsic_vol(v);
//	//std::cout << this->node_in->get_fv_volume() << '\n';
//
//	double r = this->node_in->get_fv_geom().get_outer_radius();
//	this->node_in->update_outer_radius(sqrt(vscale)*r);
//	//std::cout << this->node_in->get_fv_volume() << '\n';
}

void TreeAcinus::update_volume( const double & v) // const double & dt )    //distributes new acinar volume over tree
{
	double vol_old = this->vol;
	this->vol = v;

	//resize first results in faster code
//	std::vector<double> node_vol_old;
//	node_vol_old.resize(this->tree->count_nodes());
//	//store volume before update
//	for(size_t k = 0; k < this->tree->count_nodes(); k++)
//	{
//		node_vol_old[k] = this->tree->get_node(k)->get_fv_volume();
//	}
	double inner_vol = this->tree->get_total_inner_edge_volume();
	if(this->vol - inner_vol >= 0)   //so long as new volume is larger than duct volume
	{
		this->distribute_alveolar_volume(this->vol - inner_vol);
	}
	else
	{
		this->rescale_airways(this->vol / inner_vol);     //otherwise airway volume has to change
		this->distribute_alveolar_volume(0);    //sac volume becomes 0
	}
}

void TreeAcinus::distribute_alveolar_volume( const double & v )
{
	for(size_t j=0; j < this->tree->count_edges(); j++)   //loop over all nodes and rescale outer geom
	{
		if(this->tree->get_edge(j)->get_geom()->get_length() > 0)
		{
			double orad = sqrt(this->tree->get_edge(j)->get_geom()->get_inner_radius()
			                 * this->tree->get_edge(j)->get_geom()->get_inner_radius()
			                 + this->Nalv_edges[j] * v/(this->tree->get_edge(j)->get_geom()->get_length()
						                        * this->tree->get_edge(j)->branch_count() * M_PI));
			this->tree->get_edge(j)->update_outer_radius(orad);
		}
	}
}

void TreeAcinus::rescale_airways( const double & r )
{
	double sqrt_r = sqrt(r);
	for(size_t j=0; j < this->tree->count_edges(); j++)   //loop over all edges and rescale rad
	{
		double rad_old = this->tree->get_edge(j)->get_geom()->get_inner_radius();
		this->tree->get_edge(j)->update_inner_radius(sqrt_r * rad_old);
	}
}

void TreeAcinus::build_weibel_model(const double & s_r, const double & scale_factor, 
									  const double & lr_ratio, const double & sl_ratio)
{
	//this is a scaled version of the weibel symmetric model acinus
	using namespace network;
	Position dir(0, 0., -1);

	//Dimensions of acinar ducts from Hsia 2016
	int Nbin = this->node_in->point_count();
	int Nb[8] = {2,4,8,16,32,64,128,256};
	double wlength[8] = {1.33, 1.12, 0.93, 0.83, 0.7, 0.7, 0.7, 0.7}; //mm, will be rescaled
	double wrad[8] = {0.25, 0.245, 0.2, 0.19, 0.18, 0.17, 0.155, 0.145}; //mm, will be rescaled
	//double wrad[8] = {0.347, 0.343, 0.35, 0.356, 0.342, 0.347, 0.354, 0.359};//*Outer radius* from Haefeli-Bleuer (1988)
    //double wrad[8] = {0.27, 0.25, 0.235, 0.225, 0.22, 0.22, 0.215, 0.215};//Yeh & Schum (1980) typical-path model
    double orig_volume = 187.0;  //*(5.0/9.0) //mm3, scaled to FRC, will be rescaled
	double tot_alv_coverage[8] = {30.0, 67.0, 129.0, 219.0, 349.0, 661.0, 1204.0, 2720.0};
	double length_scale_factor = pow(this->vol/(Nbin*orig_volume), 1.0/3.0);

	std::vector<std::shared_ptr<AirwayNetworkNode>> nodes;
	nodes.push_back(this->node_in);
	std::vector<std::shared_ptr<AirwayNetworkEdge>> edges;
	double tot_Nalv = 0;
	this->Nalv_edges.clear();
	this->Nalv_edges.resize(8);
	//double Nb = this->node_in->point_count();
	/*double rad = pow(0.5, scale_factor)*s_r;
	double length = lr_ratio * rad;*/
	for (unsigned j = 1; j <= 8; j++)
	{
		double length = length_scale_factor*wlength[j-1];
		double rad = length_scale_factor*wrad[j-1];
		Position here = nodes[j-1]->get_pos();
		nodes.push_back(std::make_shared<AirwayNetworkNode>(here.x[0] + length*dir.x[0], 
			              here.x[1] + length*dir.x[1], here.x[2] + length*dir.x[2]));
		edges.push_back(std::make_shared<AirwayNetworkEdge>(nodes[j-1], nodes[j], Nb[j-1]*Nbin, rad));
		//count up total alveolar surface area
		this->Nalv_edges[j-1] = tot_alv_coverage[j-1];
		tot_Nalv += this->Nalv_edges[j-1];
	}
	this->tree = std::make_shared<FVElementNetwork>(nodes, edges);
	//this->set_alveolar_distribution();
	//define alveolar distribution
	for (unsigned j = 1; j <= 8; j++)
	{
		this->Nalv_edges[j-1] = this->Nalv_edges[j-1]/tot_Nalv;
	}


	double inner_vol = this->tree->get_total_inner_edge_volume();
	double v_left = this->vol - inner_vol;
	if(v_left <= 0)
	{
		this->rescale_airways( this->vol / inner_vol ); //updates inner vol
		v_left = 0;
	}
	this->distribute_alveolar_volume(v_left );
	//----checks----//
	//std::cout << inner_vol << ", ";
	//std::cout << inner_vol << '\n';
	//std::cout << return_tot_airway_vol() << ", " << vol << '\n';
	//std::cout << nodes.size() << ", " << edges.size() << '\n';
}

void TreeAcinus::build_dutrieue_model(const double & s_r, const double & asymmetry, const double & scale_factor, 
									  const double & lr_ratio, const double & sl_ratio)
{
	//TODO
	std::cout << "Error: this acinus model has not been implemnented.\n";
	abort_on_failure();
}

void TreeAcinus::build_henry_model(const double & s_r, const double & asymmetry, const double & scale_factor, 
									  const double & lr_ratio, const double & sl_ratio)
{
	//adapted from 
	using namespace network;

	double cutoff = 0.5;
	unsigned max_gen = 12;

	//initialise
	std::vector<std::shared_ptr<AirwayNetworkNode>> nodes;
	nodes.push_back(this->node_in);
	unsigned count = 1;   //counts node index
	std::vector<std::shared_ptr<AirwayNetworkEdge>> edges;

	std::vector<size_t> next_seed_node_indices;
	std::vector<double> ideal_theta_vals;
	next_seed_node_indices.push_back(0);
	ideal_theta_vals.push_back(0);
	std::vector<double> next_seed_radii;
	next_seed_radii.push_back(s_r);
	unsigned gen = 1;

	//loop while there are nodes and edges left to build
	while(next_seed_radii.size()>0)
	{
		//set up vectors for next iteration
		std::vector<size_t> new_seed_node_indices;
		new_seed_node_indices.reserve(2*next_seed_node_indices.size());
		std::vector<double> new_idtv;
		new_idtv.reserve(2*next_seed_node_indices.size());
		std::vector<double> new_seed_radii;
		new_seed_radii.reserve(2*next_seed_radii.size());
		double htheta = 2*M_PI/(pow(2,gen));
		//loop over seed nodes
		for(size_t n = 0; n < next_seed_node_indices.size(); n++)
		{
			size_t k_in = next_seed_node_indices[n];
			double Nb = nodes[k_in]->point_count();
			double rad, length;
			Position dir;
			Position pos0 = nodes[k_in]->get_pos();
			Position rel_pos = pos0 - node_in->get_pos();
			double r0 = rel_pos.magnitude();
			double theta0 = 0;
			if(rel_pos.x[0] != 0)
			{
				theta0 = atan2(rel_pos.x[1],rel_pos.x[0]);
			}
			for(int m = 0; m < 2; m++)
			{
				if (m==0)
				{
					rad = pow(0.5 * (1 - asymmetry), scale_factor) * next_seed_radii[n]; //min branch
				}
				else
				{
					rad = pow(0.5 * (1 + asymmetry), scale_factor)  * next_seed_radii[n]; //maj branch
				}
				length = lr_ratio * rad;
				if(rad < cutoff * s_r || gen == max_gen) length = sl_ratio * length;  //terminal sacs are longer than ducts
				
				double dtheta = (m-0.5)*htheta + ideal_theta_vals[n] - theta0;
				double r1;
				if (length < r0*sqrt(1.0 - cos(dtheta)*cos(dtheta)))
				{
					r1 = r0*cos(dtheta);
				}
				else
				{
					r1 = r0*cos(dtheta) + sqrt(length*length - r0*r0*(1-cos(dtheta)*cos(dtheta)));
				}
				dir = Position(r1*cos(theta0 + dtheta) - r0*cos(theta0), r1*sin(theta0 + dtheta) - r0*sin(theta0), 0);
				dir.normalise();

				size_t k_last = k_in;
				nodes.push_back(std::make_shared<AirwayNetworkNode>(pos0 + dir*length));
				nodes[count]->set_point_count(Nb);
				edges.push_back(std::make_shared<AirwayNetworkEdge>(nodes[k_last], nodes[count], Nb, rad));  //outer rad is set same as inner initially
				k_last = count;
				count++;

				if(rad >= cutoff * s_r && gen < max_gen)
				{
					new_seed_node_indices.push_back(k_last);
					new_seed_radii.push_back(rad);
					new_idtv.push_back(ideal_theta_vals[n] + (m-0.5)*htheta);
				}
			}
		}
		next_seed_radii = new_seed_radii;
		next_seed_node_indices = new_seed_node_indices;
		ideal_theta_vals = new_idtv;
		gen++;
	}
	this->tree = std::make_shared<FVElementNetwork>(nodes, edges);
	this->set_alveolar_distribution();

	double inner_vol = this->tree->get_total_inner_edge_volume();
	double v_left = this->vol - inner_vol;
	if(v_left <= 0)
	{
		this->rescale_airways( this->vol / inner_vol ); //updates inner vol
		v_left = 0;
	}
	this->distribute_alveolar_volume(v_left);
	//----checks----//
	//std::cout << inner_vol << ", ";
	//std::cout << inner_vol << '\n';
	//std::cout << return_tot_airway_vol() << ", " << vol << '\n';
	//std::cout << nodes.size() << ", " << edges.size() << '\n';
}

void TreeAcinus::set_alveolar_distribution()  //for branches in raw (not discretised model) set distribution of alveoli
{
	//initialise
	this->Nalv_edges.clear();
	this->Nalv_edges.resize(this->tree->count_edges());
	double Nalv_tot = 0;
	for(size_t wo = 0; wo < this->tree->count_weibel_orders(); wo++)   //loop over weibel orders
	{
		double factor;
		switch(wo)   //fraction of duct covered in alveoli as function of weibel order
		{
			case 0: factor = 0.4;
				break;

			case 1: factor = 0.7;
				break;

			default: factor = 1.0;
		}
		for(size_t jwo = 0; jwo < this->tree->count_edges_in_weibel_order(wo); jwo++)
		{
			size_t j= this->tree->get_edge_index_from_weibel_order(wo, jwo);  //node out index
			this->Nalv_edges[j] = factor*(this->tree->get_edge(j)->get_geom()->get_inner_radius())*
								  (this->tree->get_edge(j)->branch_count())*(this->tree->get_edge(j)->get_geom()->get_length()); 
									//No of alveoli proportional to surf area
			Nalv_tot += this->Nalv_edges[j];   
			//Number of alveoli in this generation of edges is proportional to surface area (i.e. radius * length * branches in generation)
			//multiplied by a factor <= 1 to allow for partial coverage in transitional bronchioles
			//However I do not understand why this factor is currently squared, need to check original source (I think Dutrieue model)
		}
	}

	if(Nalv_tot == 0 || Nalv_tot != Nalv_tot)
	{
		std::cout << "Error, Nalv_tot = " << Nalv_tot << std::endl;
		system("pause");
	}
	//normalise
	for (size_t j = 0; j < this->tree->count_edges(); j++)
	{
		this->Nalv_edges[j] = this->Nalv_edges[j] / Nalv_tot;   //fraction of alveoli in this edge
	}
}

void TreeAcinus::calc_alveolar_distribution()  //for branches in discretised model, calculate distribution of alveoli
{
	//initialise
	this->Nalv_edges.clear();
	this->Nalv_edges.resize(this->tree->count_edges());
	double Nalv_tot = 0;
	for(size_t j = 0; j < this->tree->count_edges(); j++)
	{
		this->Nalv_edges[j] = this->tree->get_edge(j)->get_outer_volume() - this->tree->get_edge(j)->get_inner_volume();
		Nalv_tot += this->Nalv_edges[j];
	}

	//normalise
	//if(Nalv_tot == 0 || Nalv_tot != Nalv_tot)
	//{
	//	std::cout << "Error: Nalv_tot = " << Nalv_tot << std::endl;
	//	std::cout << "Number of edges = " << this->tree->count_edges() << std::endl;
	//	for(size_t j = 0; j < this->tree->count_edges(); j++)
	//	{
	//		std::cout << "Edge " << j << ", outer vol = " << this->tree->get_edge(j)->get_outer_volume() << 
	//			         ", inner vol = "<< this->tree->get_edge(j)->get_inner_volume() << std::endl;
	//	}
	//	system("pause");
	//}
	for (size_t j = 0; j < this->tree->count_edges(); j++)
	{
		this->Nalv_edges[j] = this->Nalv_edges[j] / Nalv_tot;   //fraction of alveoli in this edge
	}
}

void TreeAcinus::import_and_scale(const AirwayNetwork & n, const double & sf)
{
	using namespace network;
	std::vector<std::shared_ptr<AirwayNetworkNode>>nodes;
	nodes.resize(n.count_nodes());
	std::vector<std::shared_ptr<AirwayNetworkEdge>>edges;
	edges.resize(n.count_edges());

	nodes[0] = this->node_in;
	for(size_t k=0; k < n.count_nodes(); k++)
	{
		if(k > 0) 
		{
			nodes[k] = std::make_shared<AirwayNetworkNode>(*(n.get_node(k)));
			nodes[k]->set_pos(nodes[0]->get_pos() + (n.get_node(k)->get_pos() - n.get_entry_node()->get_pos())*sf);
			nodes[k]->set_point_count(this->node_in->point_count() * n.get_node(k)->point_count());
		}
	}

	for(size_t j = 0; j < n.count_edges(); j++)
	{
		size_t k_in = n.get_node_in_index(j);
		size_t k_out = n.get_node_out_index(j);
		edges[j] = std::make_shared<AirwayNetworkEdge>(nodes[k_in], nodes[k_out], n.get_edge(j)->branch_count()*this->node_in->point_count(), 
			                               sf*n.get_edge(j)->get_geom()->get_inner_radius());
		edges[j]->update_outer_radius( sf*n.get_edge(j)->get_geom()->get_outer_radius());
	}
	this->tree = std::make_shared<FVElementNetwork>(nodes, edges);
	this->calc_alveolar_distribution();
	this->vol = this->tree->get_total_edge_volume();
}

//TODO: re-code this but every terminal node is well-mixed bag acinus

void DetailedTreeAcinus::assign_extra_nodes()
{
	//size_t Nedges = edges.size();
	//size_t Nnodes = nodes.size();
	//size_t count_new = 0;
	//for(size_t k = 0; k < Nnodes; k++)
	//{
	//	Position dir(0., 1., 0.);
	//	if(nodes[k]->get_volume() > 0)  //not just a branching point
	//	{
	//		Position new_pos = nodes[k]->get_pos() + dir * nodes[k]->get_fv_geom()->get_inner_radius();

	//		nodes.push_back(new AirwayNetworkNode(new_pos.x[0], new_pos.x[1], new_pos.x[2]));
	//		nodes.back()->set_point_count(nodes[k]->point_count());
	//		nodes.back()->set_fv_geom(new network::TubeGeometry(nodes[k]->get_fv_geom()->get_inner_radius(), 0));
	//		edges.push_back(new AirwayTransportEdge(nodes[k], nodes[Nnodes + count_new],
	//			            nodes[k]->point_count(), edges[k-1]->get_geom()->get_inner_radius(), d));
	//		nodes[Nnodes+count_new]->set_point_count(nodes[k]->point_count());
	//		count_new++;
	//	}
	//}

}

void DetailedTreeAcinus::distribute_alveolar_volume(const double & v)
{
	//check if extra nodes have been assigned

	//for(size_t k=tree->get_first_term_index(); k < tree->count_nodes(); k++)   //loop over all edges and rescale outer geom
	//{
	//	size_t j = this->tree->get_edge_in_index(k,0);
	//	size_t k_in = this->tree->get_node_in_index(j);   //parent node is fv_element

	//	//update geometry of interface (changes with expansion)
	//	double irad = this->tree->get_edge(k_in)->get_geom()->get_inner_radius();  //radius of fv tube
	//	double orad = sqrt(irad*irad + (Nalv[j]*this->vol)/(M_PI*this->tree->get_node(k_in)->get_fv_geom()->get_length()
	//		                                                    *this->tree->get_node(k_in)->point_count()));
	//	double connecting_rad;
	//	if(orad > irad)
	//	{
	//		connecting_rad = (irad / orad) * sqrt(3.0*this->tree->get_node(k_in)->get_fv_geom()->get_length()
	//			                                     *(orad * orad - irad * irad) / (orad - irad));
	//	}
	//	else
	//	{
	//		connecting_rad = sqrt(6.0*irad*this->tree->get_node(k_in)->get_fv_geom()->get_length());
	//	}
	//	this->tree->get_edge(j)->update_inner_radius(connecting_rad);
	//	this->tree->get_edge(j)->update_outer_radius(connecting_rad);

	//	//update bag geom
	//	this->tree->get_node(k)->get_fv_geom()->update_inner_radius(connecting_rad);
	//	this->tree->get_node(k)->get_fv_geom()->update_outer_radius(connecting_rad);
	//	this->tree->get_node(k)->get_fv_geom()->update_length(Nalv[k] * v / (this->tree->get_node(k)->point_count()
	//		                                                * this->tree->get_node(k)->get_fv_geom()->inner_area()));

	//}
}

std::shared_ptr<Acinus> parse_acinus_option(const char & c, std::shared_ptr<AirwayNetworkNode> nodein, 
											const double & s_r, const double & v0,
											const char & model_option, DiscretisationOperator *disc_op, 
											SimParameterList *p) //const double & s_r, const int & fv = 0, const double & df = 0, const double & asymm = 0)
{
	switch(c)
	{
	case WELL_MIXED_BAG:
		{
			return (std::make_shared<WellMixedBagAcinus>(nodein, v0));
		} break;

	case TREE_ACINUS:
		{
			std::shared_ptr<TreeAcinus> ptr = std::make_shared<TreeAcinus>(nodein, v0);
			ptr->construct_from_scratch(disc_op, s_r,
									   p->get_param<double>(ACIN_ASYMM_FACTOR_KEY)->get_value(),
									   p->get_param<double>(ACIN_LR_RATIO_KEY)->get_value(),
									   p->get_param<double>(ACIN_SCALE_FACTOR_KEY)->get_value(), 1.0,
									   model_option);
			return ptr;
		} break;

	case DETAILED_ACINUS:
		{
			std::cout << "Detailed acinus option needs re-writing. Aborting...\n";
			abort_on_failure();
			return NULL;
			//return (new DetailedTreeAcinus(nodein, p->get_param<int>(ACIN_FV_PER_GEN_KEY)->get_value(), 
			//	                               p->get_param<double>(MAX_FV_LENGTH_KEY)->get_value(), s_r,
			//						           p->get_param<double>(ACIN_ASYMM_FACTOR_KEY)->get_value(),
			//						           p->get_param<double>(ACIN_LR_RATIO_KEY)->get_value(),
			//						           p->get_param<double>(ACIN_SCALE_FACTOR_KEY)->get_value(), 
			//								   v0, model_option));
		} break;

	default:
		{
			std::cout << "This Acinus option has not been coded yet.\n";
			abort_on_failure();
			return NULL;
		}
	}
}

std::shared_ptr<Acinus> import_acin_template(const char & c, const double & sf, DiscretisationOperator *disc_op,
							 std::shared_ptr<AirwayNetworkNode> nodein, const AirwayNetwork & acin_template)
{
	switch(c)
	{
	case DETAILED_ACINUS:
		{
			return (std::make_shared<DetailedTreeAcinus>(acin_template, nodein, disc_op, sf));
		} break;
	default:
		{
			return (std::make_shared<TreeAcinus>(acin_template, nodein, disc_op, sf));
		} break;
	}
}


