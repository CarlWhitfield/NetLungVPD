#ifndef FV_ELEMENT_NETWORK_H
#define FV_ELEMENT_NETWORK_H
#include <network_3D.h>
#include "airway_edges_nodes.h"

//turn this into a class where updating edge radii results in updating node radii
//may need to redefine geom as private for original edge classs


class DiscretisationOperator
{
public:
	DiscretisationOperator(){};
	virtual int return_n_elements(const double & length, const double & rad){ return 0; };
};

class MaxMinDiscretisationOperator : public DiscretisationOperator
{
protected:
	int minfv;
	double dx;
public:
	MaxMinDiscretisationOperator(const int & minfv, const double & dx) : DiscretisationOperator()
	{
		this->minfv = minfv;
		this->dx = dx;
	}
	int return_n_elements(const double & length, const double & rad);
};

class SchererPeDiscretisationOperator : public DiscretisationOperator
{
protected:
	double Pemax;
public:
	SchererPeDiscretisationOperator(const double & Pemax) : DiscretisationOperator()
	{
		this->Pemax = Pemax;
	}
	int return_n_elements(const double & length, const double & rad);
};


class FVElementNetwork: public AirwayNetwork
{
public:
	FVElementNetwork():AirwayNetwork(){};
	FVElementNetwork(std::vector< std::shared_ptr<AirwayNetworkNode>> nodes, 
		             std::vector< std::shared_ptr<AirwayNetworkEdge>> edges):AirwayNetwork(nodes,edges){};
	FVElementNetwork(const AirwayNetwork & t, DiscretisationOperator *disc)
	{
		this->discretise_for_transport(t, disc);
	}
	FVElementNetwork(const AirwayNetwork & t, DiscretisationOperator *disc,
		             std::vector<std::vector<std::shared_ptr<AirwayNetworkEdge>>> & old_edges_to_new)
	{
		this->discretise_for_transport(t, disc, old_edges_to_new);
	}

	void convert_network_in(AirwayNetwork * t);

	void discretise_for_transport(const AirwayNetwork & t, DiscretisationOperator *disc);
	void discretise_for_transport(const AirwayNetwork & t, DiscretisationOperator *disc,
		                          std::vector<std::vector<std::shared_ptr<AirwayNetworkEdge>>> & old_edges_to_new);


	//inline void update_all_fv_element_rads()
	//{
	//	for(size_t k = 0; k < this->count_nodes(); k++)
	//	{
	//		if(this->get_node(k)->get_fv_geom().get_length() != 0)
	//		{
	//			size_t j_in = this->get_edge_in_index(k,0);
	//			this->get_node(k)->update_inner_radius(this->get_edge(j_in)->get_geom()->get_inner_radius());
	//		}
	//	}
	//}
};

#endif