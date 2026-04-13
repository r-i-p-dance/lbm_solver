// lbm_solver.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>

int main()
{   
    
    // D2Q9 Lattice Constants
    // The 9 direction vectors
    const int cx[9] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
    const int cy[9] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };

    // The 9 weights corresponding to the directions above
    const double w[9] = { 4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36 };

	// Opposite indices for streaming step
    const int opposite[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

    // Simulation params
    int Nx = 400;
    int Ny = 100;
    double tau = 0.52;
	double rho0 = 1.0;      // Initial density



    // Create initial distribution function
	std::vector<double> f(Nx * Ny * 9, 1.0); 

    // Initialize to equilibrium distribution
    for (int q = 0; q < 9; q++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
				f[q * Nx * Ny + i * Ny + j] = w[q] * rho0; 
            }
        }
    }

    for (int q = 0; q < 9; q++) {
        std::cout << "f[" << q << "] at (0,0): " << f[q * Nx * Ny + 0 * Ny + 0] << std::endl;
    }
    
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu
