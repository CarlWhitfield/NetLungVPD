#ifndef DIFFUSIVITY_H
#define DIFFUSIVITY_H

#include "network_3D.h"
#include<algorithm>

class DiffusivityObject
{
protected:
	double molecular_diffusivity;
public:
	DiffusivityObject(const double & md){ this->molecular_diffusivity = md; }
	virtual double calc_diffusivity(const network::TubeGeometry *, const double & flux) const { return 0; }
    virtual double get_molecular_diffusivity() const { return this-> molecular_diffusivity; }
};

class MolecularDiffusivityCalc: public DiffusivityObject
{
public:
	 MolecularDiffusivityCalc(const double & md):DiffusivityObject(md){};
	 inline double calc_diffusivity(const network::TubeGeometry *, const double & flux) const { return this->molecular_diffusivity; }
};

class TaylorDiffusivityCalc: public DiffusivityObject
{
public:
	TaylorDiffusivityCalc(const double & md):DiffusivityObject(md){};
	double calc_diffusivity(const network::TubeGeometry *, const double & flux) const ;
};

class SchererDiffusivityCalc: public DiffusivityObject
{
public:
	SchererDiffusivityCalc(const double & md):DiffusivityObject(md){};
	double calc_diffusivity(const network::TubeGeometry *, const double & flux) const ;
};

class ModifiedSchererDiffusivityCalc : public DiffusivityObject
{
public:
	ModifiedSchererDiffusivityCalc(const double & md) :DiffusivityObject(md){};
	double calc_diffusivity(const network::TubeGeometry *, const double & flux) const;
};

class TaulbeeTrumpetDiffusivityCalc : public DiffusivityObject
{
public:
    TaulbeeTrumpetDiffusivityCalc(const double & md) :DiffusivityObject(md){};
    double calc_diffusivity(const network::TubeGeometry *, const double & flux) const;
};

std::shared_ptr<DiffusivityObject> parse_diff_option(const char & c, const double & md, const bool & particles);

#endif