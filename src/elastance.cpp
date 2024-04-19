#if defined(_WIN32) || defined(_WIN64)
	#ifndef _USE_MATH_DEFINES
		#define _USE_MATH_DEFINES
	#endif
#endif

#include "elastance.h"
#include "define_codes.h"
#include "globals.h"

double TawhaiElastance::calc_elastance(const double & ext, const double & K0)
{
	double extsq = ext*ext;
    double extsq0 = (this->ext0)*(this->ext0);
	double z = (3*(this->a) + (this->b));
	double gamma = 0.75*z*(extsq - 1);
	double gamma0 = 0.75*z*(extsq0 - 1);

	double K = K0 * exp(gamma) * (3*z*z*(extsq-1)*(extsq-1)/extsq + z*(extsq+1)/(extsq*extsq))
		         / (exp(gamma0) *(3*z*z*(extsq0-1)*(extsq0-1)/extsq0 + z*(extsq0+1)/(extsq0*extsq0)));
	return K;
}

double TawhaiElastance::elastance_ext_grad(const double & ext, const double & K0)
{
	double extsq = ext*ext;
    double extsq0 = (this->ext0)*(this->ext0);
	double z = (3*(this->a) + (this->b));
	double gamma = 0.75*z*(extsq - 1);
	double dgdx = 1.5*z*ext;
	double gamma0 = 0.75*z*(extsq0 - 1);

	double dKdx = K0 * (dgdx * exp(gamma) * (3*z*z*(extsq-1)*(extsq-1)/extsq + z*(extsq+1)/(extsq*extsq))
		               + exp(gamma) * (12*z*z*(extsq-1) - 6*z*z*(extsq-1)*(extsq-1)/extsq)/ext
					   + exp(gamma) * (2*z - 4*(extsq+1)/extsq)/(ext*extsq))
		         / (exp(gamma0) *(3*z*z*(extsq0-1)*(extsq0-1)/extsq0 + z*(extsq0+1)/(extsq0*extsq0)));
	return dKdx;
}

std::shared_ptr<ElastanceObject> parse_elast_option(const char & c)
{
	std::shared_ptr<ElastanceObject> efunc;
	switch(c)
	{
	case LINEAR_ELASTICITY:
		{
			efunc = std::make_shared<LinearElastance>();
		} break;

	case TAWHAI_ELASTICITY:
		{
			efunc = std::make_shared<TawhaiElastance>(0.433, -0.661, pow(2.0,1.0/3.0));
		} break;

	default:
		{
			std::cout << "Did not recognise elastance option. \n";
			abort_on_failure();
		} break;
	}

	return efunc;
}


	//inline double calc_elastance_for_vol(const double & volume)
	//{
	//	double ext = pow(volume/(0.5*this->intrinsic_volume), 1.0/3.0);
	//	return (this->elastance_calc->calc_elastance(ext, this->elast0));
	//}
	//inline double calc_elastance_vol_grad(const double & volume)
	//{
	//	double ext = pow(volume/(0.5*this->intrinsic_volume), 1.0/3.0);
	//	return ( ( (1.0/3.0)*ext/volume ) * ( this->elastance_calc->elastance_ext_grad(ext, this->elast0)) );
	//}