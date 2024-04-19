#ifndef	DEFINE_OPTIONS_H
#define DEFINE_OPTIONS_H

#include "list_template.h"
#include "option_defaults.h"
#include "define_codes.h"
#include "network_3D.h"
#include <string>

class SimOptionList: public inlist::OptionList<char,bool>
{
public:
	SimOptionList();
	network::Position get_grav_direction();
};



#endif
