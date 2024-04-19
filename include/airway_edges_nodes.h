#ifndef AIRWAY_EDGES_NODES_H
#define AIRWAY_EDGES_NODES_H

#include "define_codes.h"
#include "network_3D.h"
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <Eigen/Sparse>

class AirwayNetworkNode: public network::Node 
{
	private:
    double intrinsic_vol;
	inline void blank_init()
	{
		this->seed_rad = -1;
		this->original_pos = &(this->pos);
		this->original_Npts = &(this->Npts);
		this->intrinsic_vol = 0;
	}

protected:
	network::Position alt_pos;
	network::Position *original_pos;
	double *original_Npts;
	double alt_npts;
	double seed_rad;
	//these point to pos and Npts unless these values have been substituted by averaging 
public:
	AirwayNetworkNode():Node()
	{
		this->blank_init();
	}
	AirwayNetworkNode(const double &x, const double &y, const double &z):Node(x,y,z)
	{
		this->blank_init();
	}
	AirwayNetworkNode(const network::Position & p):Node(p)
	{
		this->blank_init();
	}

	inline network::Position get_original_pos() const { return (*(this->original_pos)); }
	inline network::Position get_alt_pos() const { return (this->alt_pos); }
	inline bool has_pos_changed() const { return (this->original_pos != &(this->pos)); }
	inline double get_original_npts() const { return (*(this->original_Npts)); }
	inline double get_alt_npts() const { return (this->alt_npts); }
	inline bool has_npts_changed() const { return (this->original_Npts != &(this->Npts)); }
	inline bool has_seed_rad_changed() const { return (this->seed_rad > 0); }
	inline double get_seed_rad() const { return (this->seed_rad); }
//	inline const network::TubeGeometry& get_fv_geom() const { return this->fv_geom; }
//	inline double get_fv_volume() const { return (this->Npts * this->fv_geom.outer_volume()); }
//	inline double get_fv_inner_volume() const { return (this->Npts * this->fv_geom.inner_volume()); }

	inline void copy_all_vals(AirwayNetworkNode* n)
	{
		Node::copy_node_vals(n);
		this->seed_rad = n->get_seed_rad();
		this->alt_pos = n->get_alt_pos();
		this->alt_npts = n->get_alt_npts();
		if(n->has_pos_changed())
		{
			original_pos = &(this->alt_pos);
		}
		else
		{
			original_pos = &(this->pos);
		}
		if(n->has_npts_changed())
		{
			original_Npts = &(this->alt_npts);
		}
		else
		{
			original_Npts = &(this->alt_npts);
		}
	}

	inline void set_original_pos(const network::Position & p)
	{ 
		this->alt_pos = p; 
		this->original_pos = &(this->alt_pos);
	}
	inline void set_original_npts(const double & c)
	{ 
		this->alt_npts = c;
		this->original_Npts = &(this->alt_npts);
	}
	inline void set_seed_rad(const double & r){ this->seed_rad = r; }
	
	inline void set_intrinsic_vol(const double & vol){ this->intrinsic_vol = vol; }
    inline double get_intrinsic_vol(){ return this->intrinsic_vol; }

//	inline void set_fv_geom(const double & rad, const double & len){ this->fv_geom = network::TubeGeometry(rad,len); }
//	inline void set_fv_geom(const double & rad, const double & len, const double & orad)
//	                                                           { this->fv_geom = network::TubeGeometry(rad,len,orad); }
//	inline void update_inner_radius(const double & rad){ this->fv_geom.update_inner_radius(rad); }
//	inline void update_outer_radius(const double & orad){ this->fv_geom.update_outer_radius(orad); }
};

class AirwayNetworkEdge: public network::Edge<AirwayNetworkNode>
{
private:
	inline void init_blank()
	{
		this->original_geom = this->geom;
	}
protected:
	std::shared_ptr<network::TubeGeometry> original_geom;
public:
	AirwayNetworkEdge(std::shared_ptr<AirwayNetworkNode> ni, std::shared_ptr<AirwayNetworkNode> no)
		:Edge<AirwayNetworkNode>(ni,no)
	{
		this->init_blank();
	} 
	AirwayNetworkEdge(std::shared_ptr<AirwayNetworkNode> ni, std::shared_ptr<AirwayNetworkNode> no, 
		              const double &Nb, const double &rad)
		:Edge<AirwayNetworkNode>(ni,no,Nb,rad)
	{
		this->init_blank();
	}
	AirwayNetworkEdge(std::shared_ptr<AirwayNetworkNode> ni, std::shared_ptr<AirwayNetworkNode> no, 
		              const double &Nb, const double &rad, const double & orad)
		:Edge<AirwayNetworkNode>(ni,no,Nb,rad,orad)
	{
		this->init_blank();
	}

	inline void copy_all_vals(AirwayNetworkEdge* e)
	{
		network::Edge<AirwayNetworkNode>::copy_edge_vals(e);
		if(this->has_geom_changed())
		{
			original_geom = std::make_shared<network::TubeGeometry>(*(e->get_original_geom()));
		}
	}

	inline network::TubeGeometry* get_original_geom() const { return (this->original_geom.get()); }

	inline bool has_geom_changed() const { return (this->original_geom.get() != this->geom.get()); }
	inline void set_original_geom(const std::shared_ptr<network::TubeGeometry> g){ this->original_geom = g; }

    inline pos::Position<double> get_direction_vector() const
    {
        return (this->node_out->get_pos() - this->node_in->get_pos());
    }

	inline void update_inner_radius(const double & rad)
	{ 
		this->geom->update_inner_radius(rad); 
		//this->node_out->update_inner_radius(rad);
	}
	inline void update_outer_radius(const double & orad)
	{ 
		this->geom->update_outer_radius(orad); 
		//this->node_out->update_outer_radius(orad);
	}
};

template<class N, class E> class AirwayNetworkPlus: public network::Network<N,E>
{
private:

protected:
    void fill_incidence_matrix();//Declare function to fill matrix (J)
	void process_extra_inputs(const double & length_scale);
	std::vector<std::vector<network::Position>> original_tnode_map;
public:
    Eigen::SparseMatrix<double> Incidence;//Adding Incidence matrix (J)

	AirwayNetworkPlus():network::Network<AirwayNetworkNode, AirwayNetworkEdge>(){};
	AirwayNetworkPlus(std::vector<std::shared_ptr<N>> &node, std::vector<std::shared_ptr<E>> &edge)
			     :network::Network<N,E>(node, edge){
        this->fill_incidence_matrix();
    };
	AirwayNetworkPlus(std::map<long int, std::shared_ptr<N>> &node, std::map<long int,  std::shared_ptr<E>> &edge)
			     :network::Network<N,E>(node, edge){
        this->fill_incidence_matrix();
    };
	AirwayNetworkPlus(const std::string & node_filename, const std::string & branch_filename, 
					const std::string & termnode_filename, const double & l_scale)
					:network::Network<N,E>(node_filename, branch_filename, termnode_filename, l_scale)
	{
        this->process_extra_inputs(l_scale);
        this->fill_incidence_matrix();
	}

    int print_files_for_input(const std::string &fhead, const double & length_scale) const;
	int print_files_for_input(const std::string &fhead, const double & length_scale, const std::map<char, std::vector<std::vector<double>>>
								& extra_vals) const;


	void import_terminal_node_map(const std::string & fhead, const double & length_scale);
	void print_terminal_node_map(const std::string & fhead, const double & length_scale);
    inline void get_incidence_matrix(Eigen::SparseMatrix<double> & incidence){incidence = this->Incidence;}//Function to get incidence matrix (J)

};

template<class N, class E> int AirwayNetworkPlus<N,E>::print_files_for_input(const std::string &fhead,
                                    const double & length_scale) const
{
    std::map<char, std::vector<std::vector<double>>> extra_vals =
                    std::map<char, std::vector<std::vector<double>>>();
    return this->print_files_for_input(fhead, length_scale, extra_vals);
}


template<class N, class E> int AirwayNetworkPlus<N,E>::print_files_for_input(const std::string &fhead, 
						const double & length_scale, const std::map<char, std::vector<std::vector<double>>>
						 & extra_vals) const
{
	std::map<char, std::vector<std::vector<double>>> extra_val_final = extra_vals;

	//get extra values
	bool rad_check = false;
	std::vector<std::vector<double>> originial_edge_rads;
	originial_edge_rads.resize(this->count_edges());
	for(size_t j = 0; j < this->count_edges(); j++)
	{
		if(this->get_edge(j)->has_geom_changed())
		{
			originial_edge_rads[j].push_back(
					            this->get_edge(j)->get_original_geom()->get_inner_radius());
			rad_check = true;
		}
	}
	if(rad_check) extra_val_final['r'] = originial_edge_rads;

	bool npts_check = false, seed_rad_check = false, pos_check = false;
	std::vector<std::vector<double>> original_node_pos, original_node_npts, seed_node_rad;
	original_node_pos.resize(this->count_nodes());
	original_node_npts.resize(this->count_nodes());
	seed_node_rad.resize(this->count_nodes());
	for(size_t k = 0; k < this->count_nodes(); k++)
	{
		if(this->get_node(k)->has_npts_changed())
		{
			original_node_npts[k].push_back(this->get_node(k)->get_original_npts());
			npts_check = true;
		}
		if(this->get_node(k)->has_pos_changed())
		{
			for(size_t n = 0; n < 3; n++)
			{
				original_node_pos[k].push_back(this->get_node(k)->get_original_pos().x[n]);
			}
			pos_check = true;
		}
		if(this->get_node(k)->has_seed_rad_changed())
		{
			seed_node_rad[k].push_back(this->get_node(k)->get_seed_rad()*length_scale);
			seed_rad_check = true;
		}
	}
	if(npts_check) extra_val_final['n'] = original_node_npts;
	if(pos_check) extra_val_final['p'] = original_node_pos;
	if(seed_rad_check) extra_val_final['r'] = seed_node_rad;
	
	//output
	return network::Network<N,E>::print_files_for_input(fhead, length_scale, extra_val_final);
}

template<class N, class E> void AirwayNetworkPlus<N,E>::process_extra_inputs(const double & length_scale)
{
	std::vector<char> extra_node_args = this->get_extra_node_args();
	for(size_t i_arg = 0; i_arg < extra_node_args.size(); i_arg++)
	{
		for(size_t k = 0; k < this->count_nodes(); k++)
		{
			std::vector<double> input_vals = this->get_extra_node_inputs(extra_node_args[i_arg], 
																			this->get_node(k));
			switch(extra_node_args[i_arg])
			{
			case 'p':   //original node position to follow
				{
					if(input_vals.size() >= 3)
					{
							this->get_node(k)->set_original_pos(
								 network::Position(length_scale*(input_vals[0]),
														length_scale*(input_vals[1]),
														length_scale*(input_vals[2])));
					}

				} break;

			case 'n':   //original npts to follow
				{
					if(input_vals.size() >= 1)
					{
							this->get_node(k)->set_original_npts(input_vals[0]);
					}

				} break;

			case 'r':   //seed rad to follow
				{
					if(input_vals.size() >= 1)
					{
							this->get_node(k)->set_seed_rad(length_scale*(input_vals[0]));
					}

				} break;
			default:
				{
					std::cout << "Error did not recognise node option: " << 
									extra_node_args[i_arg] << '\n';
				}
			}
		}
	}

	std::vector<char> extra_edge_args = this->get_extra_edge_args();
	for(size_t i_arg = 0; i_arg < extra_edge_args.size(); i_arg++)
	{
		for(size_t j = 0; j < this->count_edges(); j++)
		{
			std::vector<double> input_vals = this->get_extra_edge_inputs(extra_edge_args[i_arg], 
																			this->get_edge(j));
			switch(extra_edge_args[i_arg])
			{
			case 'r':   //orig rad to follow
				{
					if(input_vals.size() >= 1)
					{
						double orig_rad = length_scale*(input_vals[0]);
						this->get_edge(j)->set_original_geom(std::make_shared<network::TubeGeometry>(orig_rad, 
														this->get_edge(j)->get_geom()->get_length()));
					}

				} break;

			default:
				{
					std::cout << "Error did not recognise edge option: " << extra_edge_args[i_arg] << '\n';
				}
			}
		}
	}
}


template<class N, class E> void AirwayNetworkPlus<N,E>::import_terminal_node_map(const std::string & fhead, 
																				  const double & length_scale)
{
	std::stringstream ss;
	//file name
	ss << fhead << "." << TNODE_MAP_FILE_EXT;
	std::string filename = ss.str();
	if(!check_infile(filename))
	{
		std::vector<std::vector<double>> data = parse_csv_file<double>(filename);
		this->original_tnode_map.resize(0);
		this->original_tnode_map.resize(this->count_term_nodes());
		for(size_t l = 0; l < data.size(); l++)
		{
			//get node index
			int tnode_index = int_round<int,double>(data[l][0]);
			//add position to map
			this->original_tnode_map(this->tnode_index_map[tnode_index]).push_back
				(network::Position(data[l][1],data[l][2],data[l][3]));
		}
		//check that number of termnodes matches Nt
		for(size_t k = this->get_first_term_index(); k < this->count_nodes(); k++)
		{
			size_t kt = k - this->get_first_term_index();
			if(this->original_tnode_map[kt].size() != 
				int_round<size_t,double>(this->get_node(k)->point_count()))
			{
				std::cout << "Warning, point count does not equal original number of nodes.\n";
			}
		}
	}
	
}
template<class N, class E> void AirwayNetworkPlus<N,E>::print_terminal_node_map(const std::string & fhead, 
																				const double & length_scale)
{
	std::stringstream ss;
	//file name
	ss << fhead << "." << TNODE_MAP_FILE_EXT;
	std::string filename = ss.str().c_str();
	std::ofstream outfile;
	if(!check_outfile(filename))
	{
		outfile.open(filename,std::ios_base::out);
		for(size_t tn = 0; tn < this->original_tnode_map.size(); tn++)
		{
			for(size_t it = 0; it < this->original_tnode_map[tn].size(); it++)
			{
				outfile << tn;
				for(int k = 0; k < 3; k++) outfile << ", " << this->original_tnode_map[tn][it].x[k];
				outfile << '\n';
			}
		}
	}
	outfile.close();
}

//Fill incidence matrix (J)
template<class N, class E> void AirwayNetworkPlus<N,E>::fill_incidence_matrix()
{
    //initialise matrix
    this->Incidence = Eigen::SparseMatrix<double>(this->count_edges(),this->count_nodes());
    std::vector<Eigen::Triplet<double>> Inc_triplet(2*this->count_edges());   //for building incidence matrix (move to object)
    //fill matrix - loop over all nodes
    for(size_t j = 0; j < this->count_edges(); j++)
    {
        int node_in = this->get_node_in_index(j);
        Inc_triplet[2*j] = Eigen::Triplet<double>(((int) j), node_in, 1.0);
        int node_out = this->get_node_out_index(j);
        Inc_triplet[2*j+1] = Eigen::Triplet<double>(((int) j), node_out, -1.0);
    }

    //build matrix
    this->Incidence.setFromTriplets(Inc_triplet.begin(),Inc_triplet.end());
}

typedef AirwayNetworkPlus<AirwayNetworkNode, AirwayNetworkEdge> AirwayNetwork;

#endif