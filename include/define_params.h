#ifndef	DEFINE_PARAMS_H
#define DEFINE_PARAMS_H

#include "list_template.h"
#include "define_options.h"
#include "param_defaults.h"
#include <string>

class SimParameterList: public inlist::ParameterList<int, double>
{
public:
	SimParameterList();   //constructor for simulation
	void check_and_convert(SimOptionList *o);
};


class RatioParam: public inlist::Parameter<double>
{
public:
	RatioParam():Parameter<double>(){};
	RatioParam(const double & a, const std::string & nam):Parameter<double>(a,nam){};
	inline bool isOK()
	{
		return (value >= 0 && value <=1);
	}
};

class PhysParam: public inlist::Parameter<double>
{
protected:
	double PhysValue, SIValue;
	std::string SIunits, PhysUnits;
	double *phys_to_sim_conversion;

public:
	PhysParam():Parameter<double>(){};
	PhysParam(const double & a, const std::string & nam):Parameter<double>(a, nam)
	{ 
		PhysValue = a;
		phys_to_sim_conversion = NULL; 
	}
	PhysParam(const double & a, const std::string & nam, double * conv):Parameter<double>(a, nam)
	{
		PhysValue = a;
		set_conversion(conv);
	}

	void update_value( const double & val ); //this is to be used to update the sim value parameter directly
	inline void read(const std::string &c){this->update_value(atof(c.c_str()));}  //this is common to all phys params
	
	inline double get_phys_value() const { return PhysValue; }
    inline double get_SI_value() const { return SIValue; }
	inline std::string get_phys_units() const { return PhysUnits; };
	std::string phys_value_string();
	inline void set_conversion( double *conv )
	{ 
		phys_to_sim_conversion = conv; 
		calc_sim_from_phys_value();
	}
	inline void calc_sim_from_phys_value()
	{ 
		if(phys_to_sim_conversion != NULL) value = get_phys_value() * (*phys_to_sim_conversion);
	}

	virtual void calc_SI_from_Phys_value(){ SIValue = get_phys_value(); }
	virtual void define_units()
	{ 
		PhysUnits = ""; 
		SIunits = "";
	}
};

class SpecificPhysParam: public PhysParam
{
public:
	SpecificPhysParam():PhysParam(){};
	SpecificPhysParam(const double & a, const std::string & nam):PhysParam(a,nam){};
	SpecificPhysParam(const double & a, const std::string & nam, double * conv):PhysParam(a,nam,conv){};

	void init()
	{
		this->define_units();
		this->calc_SI_from_Phys_value();
	}

	bool isOK()
	{
		return get_phys_value() >= 0;
	}

	virtual void define_units()
	{
		std::cerr << "Warning: using virtual function\n";
	}
	virtual void calc_SI_from_Phys_value()
	{
		std::cerr << "Warning: using virtual function\n";
	}
};

class VolumeParam: public SpecificPhysParam
{
public:
	VolumeParam():SpecificPhysParam()
	{
		this->init();	
	}
	VolumeParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	VolumeParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}
	
	void define_units();
	void calc_SI_from_Phys_value();
};

class TimeParam: public SpecificPhysParam
{
public:
	TimeParam():SpecificPhysParam()
	{
		this->init();	
	}
	TimeParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	TimeParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}
	
	void define_units();
	void calc_SI_from_Phys_value();
};

class LengthParam: public SpecificPhysParam
{
public:
	LengthParam():SpecificPhysParam()
	{
		this->init();	
	}
	LengthParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	LengthParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();
};

class ReciprocalLengthParam: public SpecificPhysParam
{
public:
	ReciprocalLengthParam():SpecificPhysParam()
	{
		this->init();	
	}
	ReciprocalLengthParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	ReciprocalLengthParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();
};

class PressureParam: public SpecificPhysParam
{
public:
	PressureParam():SpecificPhysParam()
	{
		this->init();	
	}
	PressureParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	PressureParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();
};

class PressurePerLengthParam: public SpecificPhysParam
{
public:
	PressurePerLengthParam():SpecificPhysParam()
	{
		this->init();	
	}
	PressurePerLengthParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	PressurePerLengthParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();
};

class ElastanceParam: public SpecificPhysParam
{
public:
	ElastanceParam():SpecificPhysParam()
	{
		this->init();	
	}
	ElastanceParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	ElastanceParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();
};

class ResistanceParam: public SpecificPhysParam
{
public:
	ResistanceParam():SpecificPhysParam()
	{
		this->init();	
	}
	ResistanceParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	ResistanceParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();
};

class InertanceParam: public SpecificPhysParam
{
public:
	InertanceParam():SpecificPhysParam()
	{
		this->init();	
	}
	InertanceParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	InertanceParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();
};


class ViscosityParam: public SpecificPhysParam
{
public:
	ViscosityParam():SpecificPhysParam()
	{
		this->init();	
	}
	ViscosityParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	ViscosityParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();
};

class DensityParam: public SpecificPhysParam
{
public:
	DensityParam():SpecificPhysParam()
	{
		this->init();	
	}
	DensityParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	DensityParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();

};

class DiffusivityParam: public SpecificPhysParam
{
public:
	DiffusivityParam():SpecificPhysParam()
	{
		this->init();	
	}
	DiffusivityParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
	{
		this->init();	
	}
	DiffusivityParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
	{
		this->init();	
	}

	void define_units();
	void calc_SI_from_Phys_value();
};

class AccelerationParam: public SpecificPhysParam
{
public:
    AccelerationParam():SpecificPhysParam()
    {
        this->init();
    }
    AccelerationParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
    {
        this->init();
    }
    AccelerationParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
    {
        this->init();
    }

    void define_units();
    void calc_SI_from_Phys_value();
};

class ParticleDensityParam: public SpecificPhysParam
{
public:
    ParticleDensityParam():SpecificPhysParam()
    {
        this->init();
    }
    ParticleDensityParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
    {
        this->init();
    }
    ParticleDensityParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
    {
        this->init();
    }

    void define_units();
    void calc_SI_from_Phys_value();
};

class VolumeFluxParam: public SpecificPhysParam
{
public:
    VolumeFluxParam():SpecificPhysParam()
    {
        this->init();
    }
    VolumeFluxParam(const double & a, const std::string & nam):SpecificPhysParam(a,nam)
    {
        this->init();
    }
    VolumeFluxParam(const double & a, const std::string & nam, double * conv):SpecificPhysParam(a,nam,conv)
    {
        this->init();
    }

    void define_units();
    void calc_SI_from_Phys_value();
};

#endif
