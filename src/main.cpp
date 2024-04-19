// main.cpp : Defines the entry point for the console application.
//

#include "lung_simulation.h"
#include <iostream>

int main(int argc, char** argv)
{
	auto start = std::chrono::system_clock::now();

	LungSimulation sim(argc, argv);   //construct simulation environment from command line arguments

	//initialise simulation
	if (sim.initialise())
	{
		std::cout << "Error initialising simulation, aborting.\n";
		return EXIT_FAILURE;
	}
	
	//run simulation
	if (sim.simulate())
	{
		std::cout << "Error running simulation, aborting.\n";
		return EXIT_FAILURE;
	}

	auto end = std::chrono::system_clock::now();

	std::cout << "Whole code took: " << (std::chrono::duration<double>(end-start)).count() << "s.\n"; 
	std::cout << "Number of flow nodes = " << sim.flow_tree->count_nodes() << '\n';
	std::cout << "Number of transport nodes = " << sim.transport_tree->count_nodes() << '\n';
	std::cout << "Number of term nodes = " << sim.flow_tree->count_term_nodes() << '\n';

	return EXIT_SUCCESS;
}
