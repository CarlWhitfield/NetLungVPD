#if defined(_WIN32) || defined(_WIN64)
	#ifndef _USE_MATH_DEFINES
		#define _USE_MATH_DEFINES
	#endif
#endif

#include "define_codes.h"
#include "diffusivity.h"
#include "network_3D.h"
#include <math.h>

using namespace network;

double TaylorDiffusivityCalc::calc_diffusivity(const TubeGeometry * g, const double & flux) const
{
	double u_mag = fabs(flux / g->inner_area());
	return (molecular_diffusivity + u_mag * u_mag * g->inner_area() / (192 * M_PI * molecular_diffusivity));  //return taylor dispersion
}

double SchererDiffusivityCalc::calc_diffusivity(const TubeGeometry * g, const double & flux) const
{
	double u_mag = fabs(flux / g->inner_area());
	if(flux > 0)
	{
		return (this->molecular_diffusivity + (1.08*u_mag*2.0*g->get_inner_radius()));   //return scherer (1975) dispersion 
	}
	else
	{
		return (this->molecular_diffusivity + (0.37*u_mag*2.0*g->get_inner_radius()));   //return scherer (1975) dispersion 
	}
}

double ModifiedSchererDiffusivityCalc::calc_diffusivity(const TubeGeometry * g, const double & flux) const
{
	double u_mag = fabs(flux / g->inner_area());
	if (flux > 0)
	{
		return (this->molecular_diffusivity + (0.7*u_mag*2.0*g->get_inner_radius()));   //return scherer (1975) dispersion
        //return (this->molecular_diffusivity + (1.08*u_mag*2.0*g->get_inner_radius()) + 0.6*u_mag*g->get_length());   //return scherer (1975) dispersion
	}
	else
	{
		return (this->molecular_diffusivity + (0.26*u_mag*2.0*g->get_inner_radius()));   //return scherer (1975) dispersion
	}
}

double TaulbeeTrumpetDiffusivityCalc::calc_diffusivity(const TubeGeometry * g, const double & flux) const
{
    double u_mag = fabs(flux / g->inner_area());
    return (this->molecular_diffusivity + 0.6*u_mag*g->get_length());   //return Taulbee & Yu (1975) apparent diffusivity
}


std::shared_ptr<DiffusivityObject> parse_diff_option(const char & c, const double & md, const bool & particles)
{
	switch (c)
	{
	case MOLECULAR_DIFFUSIVITY:
		return (std::make_shared<MolecularDiffusivityCalc>(md));
		break;

	case TAYLOR_DISPERSION:
		return (std::make_shared<TaylorDiffusivityCalc>(md));
		break;

	case SCHERER_DISPERSION:
		if (particles) return (std::make_shared<ModifiedSchererDiffusivityCalc>(md));
		else return (std::make_shared<SchererDiffusivityCalc>(md));
		break;

        case TAULBEE_DISPERSION:
            return (std::make_shared<TaulbeeTrumpetDiffusivityCalc>(md));
            break;

	default:
		return (std::make_shared<DiffusivityObject>(md));
	}
}