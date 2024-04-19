#ifndef ELASTANCE_H
#define ELASTANCE_H

#include "network_3D.h"
#include <algorithm>

class ElastanceObject
{
public:
	ElastanceObject(){};
	virtual double calc_elastance(const double & ext, const double & K0){ return K0; }
	virtual double elastance_ext_grad(const double & ext, const double & K0){ return 0; }
};

class LinearElastance: public ElastanceObject
{
public:
	LinearElastance():ElastanceObject(){};
};

class TawhaiElastance: public ElastanceObject
{
protected:
	double a, b, ext0;

public:
	TawhaiElastance(const double & a0, const double & b0, const double & extfrc):ElastanceObject()
	{
		this->a = a0;
		this->b = b0;
		this->ext0 = extfrc;
	}
	double calc_elastance(const double & ext, const double & K0);
	double elastance_ext_grad(const double & ext, const double & K0);
};

std::shared_ptr<ElastanceObject> parse_elast_option(const char & c);

#endif