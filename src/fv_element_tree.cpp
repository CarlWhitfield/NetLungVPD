#include "fv_element_tree.h"
#include <vector>

int MaxMinDiscretisationOperator::return_n_elements(const double & length, const double & rad){

	return std::max(this->minfv, int(ceil(length / this->dx)));
}

int SchererPeDiscretisationOperator::return_n_elements(const double & length, const double & rad){
	double Pe = length / (0.52*rad);
	return ceil(Pe / this->Pemax);
}

void FVElementNetwork::discretise_for_transport(const AirwayNetwork & t, DiscretisationOperator *disc_op)
{
    std::vector<std::vector<std::shared_ptr<AirwayNetworkEdge>>> default_vec =
            std::vector<std::vector<std::shared_ptr<AirwayNetworkEdge>>>();
    this->discretise_for_transport(t, disc_op, default_vec);
}

void FVElementNetwork::discretise_for_transport(const AirwayNetwork & tree_in, DiscretisationOperator *disc,
		                             std::vector<std::vector<std::shared_ptr<AirwayNetworkEdge>>> & old_edges_to_new)
{
	using namespace network;

	//copy nodes
	this->NodeVec.clear();
	this->EdgeVec.clear();
	this->NodeVec.resize(tree_in.count_nodes());
	for(size_t k = 0; k < tree_in.count_nodes(); k++)
	{
		this->NodeVec[k] = tree_in.get_node_smart_ptr(k);
	}

	//re-discretize network to make fv nodes
	//store old network
	size_t NNodes0 = tree_in.count_nodes();
	size_t NEdges0 = tree_in.count_edges();
	old_edges_to_new.resize(NEdges0);

	//vector to count no. of finite volumes in each branch
	std::vector<int> nfv;
	size_t totfv = 0;

	for(size_t j = 0; j < NEdges0; j++)
	{
		int nfvh = disc->return_n_elements(tree_in.get_edge(j)->get_geom()->get_length(),
			tree_in.get_edge(j)->get_geom()->get_inner_radius());
		//record how many elements each edge will be split into
		nfv.push_back(nfvh);
		totfv += nfvh;
		//record nodes at start and end of transport branch
		double rad = tree_in.get_edge(j)->get_geom()->get_inner_radius();
		double orad = tree_in.get_edge(j)->get_geom()->get_outer_radius();
        //nfvh includes one of the bifurcation nodes
		size_t k_start = tree_in.get_node_in_index(j);
		size_t k_end = tree_in.get_node_out_index(j); //this needs to be assigned a volume
		Position direction = tree_in.get_edge(j)->get_direction_vector();
		std::shared_ptr<AirwayNetworkNode> last_node = this->get_node_smart_ptr(k_start);
		for(int el = 0; el < nfvh-1; el++)  //loop over fv elements (nodes)
		{
			//create extra intermediate nodes and edges
			double frac = ((double) el + 1.0) / ((double) nfvh);
			std::shared_ptr<AirwayNetworkNode> this_node = std::make_shared<AirwayNetworkNode>
				                    (tree_in.get_node(k_start)->get_pos(0) + frac*direction.x[0],
									 tree_in.get_node(k_start)->get_pos(1) + frac*direction.x[1],
									 tree_in.get_node(k_start)->get_pos(2) + frac*direction.x[2]);
			this_node->set_point_count(tree_in.get_edge(j)->branch_count());
				
			this->NodeVec.push_back(this_node);
			this->EdgeVec.push_back(std::make_shared<AirwayNetworkEdge>(last_node, this_node, 
				                             tree_in.get_edge(j)->branch_count(), rad, orad));
			old_edges_to_new[j].push_back(this->EdgeVec.back());

			last_node = this_node;
		}
		//final edge connecting to k_end       //(zero volume node at bifurcation)
		this->EdgeVec.push_back(std::make_shared<AirwayNetworkEdge>(last_node, this->get_node_smart_ptr(k_end),
			                                   tree_in.get_edge(j)->branch_count(), rad, orad));
		old_edges_to_new[j].push_back(this->EdgeVec.back());
	}
	//reorder and stuff
	this->update_node_edge_maps();
	this->reorder_network();
    this->fill_incidence_matrix();
}