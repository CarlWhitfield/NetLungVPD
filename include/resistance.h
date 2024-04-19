#ifndef RESISTANCE_H
#define RESISTANCE_H

#include "network_3D.h"

class ResistanceObject
{
protected:
	double viscosity, density;
public:
	ResistanceObject(const double & visc, const double & dens);
	virtual double calc_resistance(const network::TubeGeometry *, const double & flux) const { return 0; }
	virtual double resistance_flux_grad(const network::TubeGeometry *, const double & flux) const { return 0; }
};

class PoiseuilleResistanceCalc: public ResistanceObject
{
public:
	PoiseuilleResistanceCalc(const double & visc, const double & dens);
	double calc_resistance(const network::TubeGeometry *, const double & flux) const ;
};

class PedleyResistanceCalc: public ResistanceObject
{
public:
	PedleyResistanceCalc(const double & visc, const double & dens);
	double calc_resistance(const network::TubeGeometry *, const double & flux) const ;
	double resistance_flux_grad(const network::TubeGeometry *, const double & flux) const ;
};

std::shared_ptr<ResistanceObject> parse_res_option(const char & c, const double & visc, const double & dens);

#endif