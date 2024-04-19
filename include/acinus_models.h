#ifndef ACINUS_MODELS_H
#define ACINUS_MODELS_H

#include <network_3D.h>
#include "airway_edges_nodes.h"
#include "fv_element_tree.h"
#include "define_params.h"

class Acinus {
protected:
	std::vector<double> blank_vec;
	std::shared_ptr<AirwayNetworkNode> node_in;
	double Nterm, vol;  //no of airways terminating into this (i.e. number of acini)
	inline void node_setup(std::shared_ptr<AirwayNetworkNode> n_in)
	{
		this->node_in = n_in;
		this->tree = NULL;
		this->Nterm = n_in->point_count();
	}
public:
	Acinus(std::shared_ptr<AirwayNetworkNode> n_in)
	{
		this->node_setup(n_in);
	}
	Acinus(std::shared_ptr<AirwayNetworkNode> n_in, const double & v0)
	{ 
		this->node_setup(n_in);
		this->vol = v0;
	}

	std::shared_ptr<FVElementNetwork> tree;
	void copy_structure(Acinus *);
	int copy_vals(Acinus *);

	inline double term_count() const { return (this->Nterm); }
	inline double get_volume() const { return (this->vol); }
	inline void get_position(network::Position & pos) { pos = node_in->get_original_pos(); }
	//this gets overwritten for more complex acinus models
	virtual void build_model(){};
	virtual void update_volume(const double & volume)
	{ 
		this->vol = volume; 
	}
};

class WellMixedBagAcinus: public Acinus
{
public:
	WellMixedBagAcinus(std::shared_ptr<AirwayNetworkNode> n_in, const double & v0):Acinus(n_in, v0)
	{ 
		this->build_model(); 
	}
	WellMixedBagAcinus(const Acinus & a):Acinus(a)
	{ 
		this->build_model(); 
	};

	void build_model(); 
	void update_volume(const double & volume);
};

class TreeAcinus: public Acinus
{
protected:
	std::vector<double> Nalv_edges, edge_fluxes;
	void import_and_scale(const AirwayNetwork & n, const double & sf);
	inline void setup(DiscretisationOperator *disc_op)
	{
		FVElementNetwork tree_copy = *(this->tree);
		this->tree->discretise_for_transport(tree_copy, disc_op);

		this->calc_alveolar_distribution();
		this->edge_fluxes.resize(this->tree->count_edges());
		std::fill(this->edge_fluxes.begin(), this->edge_fluxes.end(), 0);
	}
public:
	TreeAcinus(std::shared_ptr<AirwayNetworkNode> n, const double & v0):Acinus(n, v0){};   //constructor, just set node_in 

	//constructor, import network from template using pre-existing acinus
	TreeAcinus(const AirwayNetwork & n, const Acinus & a,  
				DiscretisationOperator *disc_op, const double & import_scale) :Acinus(a)
	{
		this->import_and_scale(n, import_scale);
		this->setup(disc_op);
	}

	//constructor to import from scratch
	TreeAcinus(const AirwayNetwork & n, std::shared_ptr<AirwayNetworkNode> node, 
		DiscretisationOperator *disc_op, const double & import_scale)
			   :Acinus(node)
	{
		this->import_and_scale(n, import_scale);
		this->setup(disc_op);
	}

	//constructor to build from scratch
	TreeAcinus(std::shared_ptr<AirwayNetworkNode> node, DiscretisationOperator *disc_op,
			   const double & s_r, const double & asymm, const double & lr, const double & sf, const double & sl_ratio, 
			    const double & v0, const char & model_opt):Acinus(node, v0)
	{
		this->construct_from_scratch(disc_op, s_r, asymm, lr, sf, sl_ratio, model_opt);
	}

	//constructor to build from pre-existing acinus
	TreeAcinus(const Acinus & a, DiscretisationOperator *disc_op, const double & s_r,
		       const double & asymm, const double & lr, const double & sf, const double & sl_ratio, 
			   const char & model_opt):Acinus(a)
	{
		this->build_model(s_r, asymm, sf, lr, sl_ratio, model_opt);
		this->setup(disc_op);
	};

	inline void construct_from_scratch(DiscretisationOperator *disc_op, const double & s_r,
		            const double & asymm, const double & lr, const double & sf, const double & sl_ratio, 
			        const char & model_opt)
	{
		this->build_model(s_r, asymm, sf, lr, sl_ratio, model_opt);
		this->setup(disc_op);	
	}

	void update_volume(const double & v);
	void rescale_airways(const double & ratio);
	void set_alveolar_distribution(  );
	void calc_alveolar_distribution();
	inline void build_model(const double & s_r, const double & asymm, const double & sf, 
		             const double & lr, const double & sl_ratio, const char & model_opt)
	{
		switch(model_opt)
		{
		case WEIBEL_ACIN:
			{
				return(this->build_weibel_model(s_r, sf, lr, sl_ratio));
			} break;

		case DUTRIEUE_ACIN:
			{
				return(this->build_dutrieue_model(s_r, asymm, sf, lr, sl_ratio));
			} break;

		case HENRY_ACIN:
			{
				return(this->build_henry_model(s_r, asymm, sf, lr, sl_ratio));
			} break;
			
		default:
			{
				std::cout << "Error: could not parse acinus model option.\n";
			}
		}
	}
	void build_weibel_model(const double & s_r, const double & sf, const double & lr, const double & sl_ratio);
	void build_dutrieue_model(const double & s_r, const double & asymm, const double & sf, const double & lr, const double & sl_ratio);
	void build_henry_model(const double & s_r, const double & asymm, const double & sf, const double & lr, const double & sl_ratio);

	virtual void distribute_alveolar_volume( const double & v );

    inline std::vector<double> get_Nalv_edges() const { return (this->Nalv_edges); };
};

class DetailedTreeAcinus: public TreeAcinus
{
protected:
	void assign_extra_nodes();
	//TODO: write some initialisation routine for this
	std::vector<WellMixedBagAcinus*> alveoli;
public:
	DetailedTreeAcinus(std::shared_ptr<AirwayNetworkNode> n, const double & v0):TreeAcinus(n, v0){};

	DetailedTreeAcinus(const AirwayNetwork & n, const Acinus & a,  
						DiscretisationOperator *disc_op, const double & import_scale)
					   :TreeAcinus(n,a,disc_op,import_scale){};

	DetailedTreeAcinus(const Acinus & a, DiscretisationOperator *disc_op,
		       const double & s_r, const double & asymm, const double & lr, const double & sf, const double & sl_ratio, 
			   const char & model_opt):TreeAcinus(a, disc_op, s_r, asymm, lr, sf, sl_ratio, model_opt){};

	DetailedTreeAcinus(std::shared_ptr<AirwayNetworkNode> node, DiscretisationOperator *disc_op,
			   const double & s_r, const double & asymm, const double & lr, const double & sf, const double & sl_ratio, 
			   const double & v0, const char & model_opt) 
			   :TreeAcinus(node,disc_op,s_r,asymm,lr,sf,sl_ratio,v0,model_opt){};

	DetailedTreeAcinus(const AirwayNetwork & n, std::shared_ptr<AirwayNetworkNode> node,
					DiscretisationOperator *disc_op,
					   const double & import_scale):TreeAcinus(n,node,disc_op,import_scale){};

	void distribute_alveolar_volume( const double & );
};

std::shared_ptr<Acinus> parse_acinus_option(const char & c, std::shared_ptr<AirwayNetworkNode> nodein, const double & s_r,
					const double & v0, const char & model_option, DiscretisationOperator *disc_op, 
					SimParameterList *p);
std::shared_ptr<Acinus> import_acin_template(const char & c, const double & sf, DiscretisationOperator *disc_op, std::shared_ptr<AirwayNetworkNode> nodein,
							 const AirwayNetwork & acin_template);

#endif
