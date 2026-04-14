// lbm_solver.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>

int main()
{   
    // -------------------------------------------------------------
    // Define constants and params
    // -------------------------------------------------------------
    
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



    // -------------------------------------------------------------
    // Create initial distribution function
    // -------------------------------------------------------------
	std::vector<double> f(Nx * Ny * 9, 1.0); 

    // Initialize to equilibrium distribution
    for (int q = 0; q < 9; q++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
				f[q * Nx * Ny + i * Ny + j] = w[q] * rho0;
            }
        }
    }
	// Test print the distribution function at (0,0)
    std::cout << "Test print the distribution function: " << std::endl;
    for (int q = 0; q < 9; q++) {
        std::cout << "f[" << q << "] at (0,0): " << f[q * Nx * Ny + 0 * Ny + 0] << std::endl;
    }
    std::cout << std::endl;



	// -------------------------------------------------------------
    // Density field
    // -------------------------------------------------------------
	std::vector<double> rho(Nx * Ny, 0.0); 

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int q = 0; q < 9; q++) {
                rho[i * Ny + j] += f[q * Nx * Ny + i * Ny + j];
            }
        }
    }
	// Test density field (should be 1 at all points)
    std::cout << "Test density field (should be 1 at all points)" << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            std::cout << "rho[" << i << ", " << j << "]: " << rho[i * Ny + j] << std::endl;
    }
    std::cout << std::endl;
    
    // -------------------------------------------------------------
    // Compute velocity u 
    // -------------------------------------------------------------
    std::vector<double> u_x(Nx * Ny, 0.0);
    std::vector<double> u_y(Nx * Ny, 0.0);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int q = 0; q < 9; q++) {
                u_x[i * Ny + j] += f[q * Nx * Ny + i * Ny + j] * cx[q];
				u_y[i * Ny + j] += f[q * Nx * Ny + i * Ny + j] * cy[q];
            }
			u_x[i * Ny + j] /= rho[i * Ny + j];
            u_y[i * Ny + j] /= rho[i * Ny + j];
        }
    }
    // Test velocity field 
    std::cout << "Test velocity field" << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << "u_x[" << i << ", " << j << "]: " << u_x[i * Ny + j] << std::endl;
            std::cout << "u_y[" << i << ", " << j << "]: " << u_y[i * Ny + j] << std::endl;
        }
    }
    std::cout << std::endl;

	// -------------------------------------------------------------
	// Compute equilibrium distribution function feq
    // -------------------------------------------------------------
    // Compute c_u
    std::vector<double> c_u(Nx * Ny * 9, 0.0);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int q = 0; q < 9; q++) {
                c_u[q * Nx * Ny + i * Ny + j] = cx[q] * u_x[i * Ny + j] + cy[q] * u_y[i * Ny + j];
            }
        }
    }

    // Compute u_sq
    std::vector<double> u_sq(Nx * Ny, 0.0);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            u_sq[i * Ny + j] = u_x[i * Ny + j] * u_x[i * Ny + j] + u_y[i * Ny + j] * u_y[i * Ny + j];
        }
    }

	// Compute feq
	std::vector<double> f_eq(Nx * Ny * 9, 0.0);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int q = 0; q < 9; q++) {
                f_eq[q * Nx * Ny + i * Ny + j] = w[q] * rho[i * Ny + j] * (1 + 3 * c_u[q * Nx * Ny + i * Ny + j] + 4.5 * c_u[q * Nx * Ny + i * Ny + j] * c_u[q * Nx * Ny + i * Ny + j] - 1.5 * u_sq[i * Ny + j]);
            }
        }
    }

    // Test f_eq
    std::cout << "Test f_eq" << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << "At(" << i << ", " << j << ") : " << std::endl;
            for (int q = 0; q < 9; q++) {
                std::cout << "f_eq[" << q << "]: " << f_eq[q * Nx * Ny + i * Ny + j] << std::endl;
            }
            std::cout << std::endl;
        }
    }

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu
