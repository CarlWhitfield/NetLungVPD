#if defined(_WIN32) || defined(_WIN64)
	#ifndef _USE_MATH_DEFINES
		#define _USE_MATH_DEFINES
	#endif
#endif

#include "resistance.h"
#include "network_3D.h"
#include "math.h"
#include "define_codes.h"

ResistanceObject::ResistanceObject(const double & visc, const double & dens)
{
	this->viscosity = visc;
	this->density = dens;
}

PoiseuilleResistanceCalc::PoiseuilleResistanceCalc(const double & visc, const double & dens)
	                     :ResistanceObject(visc, dens){}

PedleyResistanceCalc::PedleyResistanceCalc(const double & visc, const double & dens)
	                 :ResistanceObject(visc, dens){}

double PoiseuilleResistanceCalc::calc_resistance(const network::TubeGeometry * geom, const double & flux=0)  const 
{
	double rad = geom->get_inner_radius();
	return 8 * viscosity * geom->get_length() / (M_PI * rad * rad * rad * rad);
}

double PedleyResistanceCalc::calc_resistance(const network::TubeGeometry * geom, const double & flux) const 
{
	double Rpois = PoiseuilleResistanceCalc(viscosity, density).calc_resistance(geom);
	double flux_mag = abs(flux);
	double Z = 0.925 *sqrt(flux_mag * density / (2.0 * M_PI * viscosity * geom->get_length()));
	if (Z > 1) return Z*Rpois;
	else return Rpois;
}

double PedleyResistanceCalc::resistance_flux_grad(const network::TubeGeometry * geom, const double & flux) const 
{
	double Rpois = PoiseuilleResistanceCalc(viscosity, density).calc_resistance(geom);
	double flux_mag = sqrt(flux * flux);
	double Z = 0.925 *sqrt(flux_mag * density / (2.0 * M_PI * viscosity * geom->get_length()));
	if (Z > 1) return 0.5 * Z * Rpois / flux;
	else return 0;
}

std::shared_ptr<ResistanceObject> parse_res_option(const char & c, const double & visc, const double & dens)
{
	switch(c)
	{
	case POISEUILLE:
		return std::make_shared<PoiseuilleResistanceCalc>(visc, dens);
		break;

	case PEDLEY:
		return std::make_shared<PedleyResistanceCalc>(visc, dens);
		break;

	default:
		return std::make_shared<ResistanceObject>(visc, dens);
		break;
	}
}