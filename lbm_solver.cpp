// lbm_solver.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>

// -------------------------------------------------------------
// Constants — global so all functions can see them
// -------------------------------------------------------------
// 
// D2Q9 Lattice Constants
    // The 9 direction vectors
const int cx[9] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
const int cy[9] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };

// The 9 weights corresponding to the directions above
const double w[9] = { 4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36 };

// Opposite indices for streaming step
const int opposite[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

// Inlet velocity
double u_inlet = 0.04;

// -------------------------------------------------------------
// Functions
// -------------------------------------------------------------
void compute_macroscopic(const std::vector<double>& f,
    std::vector<double>& rho,
    std::vector<double>& u_x,
    std::vector<double>& u_y,
    int Nx, int Ny) {

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = i * Ny + j;

            rho[idx] = 0.0;
            u_x[idx] = 0.0;
            u_y[idx] = 0.0;

            for (int q = 0; q < 9; q++) {
                double fq = f[q * Nx * Ny + idx];
                rho[idx] += fq;
                u_x[idx] += fq * cx[q];
                u_y[idx] += fq * cy[q];
            }

            u_x[idx] /= rho[idx];
            u_y[idx] /= rho[idx];
        }
    }
}

void collide(std::vector<double>& f,
    const std::vector<double>& rho,
    const std::vector<double>& u_x,
    const std::vector<double>& u_y,
    double tau, int Nx, int Ny) {

    // Do not create arrays for c_u, u_sq, f_eq, as they are not needed outside this function. Just compute them on the fly for every cell.
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = i * Ny + j;
            double usq = u_x[idx] * u_x[idx] + u_y[idx] * u_y[idx];

            for (int q = 0; q < 9; q++) {
                double cu = cx[q] * u_x[idx] + cy[q] * u_y[idx];
                double feq = w[q] * rho[idx] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * usq);

                int fidx = q * Nx * Ny + idx;
                f[fidx] -= (1.0 / tau) * (f[fidx] - feq);
            }
        }
    }
}

void stream(std::vector<double>& f,
    std::vector<double>& f_new,
    int Nx, int Ny) {

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int q = 0; q < 9; q++) {
                // Where does this particle come FROM?
                int i_src = (i - cx[q] + Nx) % Nx;
                int j_src = (j - cy[q] + Ny) % Ny;

                f_new[q * Nx * Ny + i * Ny + j] = f[q * Nx * Ny + i_src * Ny + j_src];
            }
        }
    }
    // Swap buffers
    std::swap(f, f_new);

}

void apply_outlet_bc(std::vector<double>& f, int Nx, int Ny) {
    for (int j = 0; j < Ny; j++) {
        for (int q = 0; q < 9; q++) {
            f[q * Nx * Ny + (Nx - 1) * Ny + j] = f[q * Nx * Ny + (Nx - 2) * Ny + j];
        }
    }
}

void apply_inlet_bc(std::vector<double>& f,
    std::vector<double>& rho,
    std::vector<double>& u_x,
    std::vector<double>& u_y,
    int Nx, int Ny) {
    for (int j = 0; j < Ny; j++) {

        u_x[j] = u_inlet; // Set a constant velocity at the inlet
        u_y[j] = 0.0;
        rho[j] = 1.0; // Set a constant density at the inlet

        double u_sq_inlet = u_x[j] * u_x[j];
        for (int q = 0; q < 9; q++) {
            double cu_inlet = cx[q] * u_inlet;
            f[q * Nx * Ny + j] = w[q] * rho[j] * (1.0 + 3.0 * cu_inlet + 4.5 * cu_inlet * cu_inlet - 1.5 * u_sq_inlet);
        }
    }
}

void bounce_back(std::vector<double>& f,
	const std::vector<bool>& obstacle,
	int Nx, int Ny) {

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (obstacle[i * Ny + j]) {
				for (int q = 0; q < 9; q++) {
					int fidx = q * Nx * Ny + i * Ny + j;
					int fidx_opposite = opposite[q] * Nx * Ny + i * Ny + j;
					std::swap(f[fidx], f[fidx_opposite]);
				}
			}
		}
	}
}



int main()
{   
    // Simulation params
    int Nx = 400;
    int Ny = 100;
	double rho0 = 1.0;      // Initial density
    double tau = 0.52;

    std::vector<double> f(Nx * Ny * 9, 1.0);
    std::vector<double> f_new(9 * Nx * Ny, 0.0);
    std::vector<double> rho(Nx * Ny, 0.0);
    std::vector<double> u_x(Nx * Ny, 0.0);
    std::vector<double> u_y(Nx * Ny, 0.0);   


    
    // Create initial distribution function
    for (int q = 0; q < 9; q++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                f[q * Nx * Ny + i * Ny + j] = w[q] * rho0;
            }
        }
    }

	// Initialize obstacle (circular cylinder in the middle of the domain)
    std::vector<bool> obstacle(Nx * Ny, false);
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            obstacle[i * Ny + j] = (i - Nx / 4) * (i - Nx / 4) + (j - Ny / 2) * (j - Ny / 2) < (Ny / 10) * (Ny / 10);
        }
    }
    for (int i = 0; i < Nx; i++) {
        obstacle[i * Ny + 0] = true; // Bottom wall
        obstacle[i * Ny + (Ny - 1)] = true; // Top wall
    }

	apply_outlet_bc(f, Nx, Ny);
	compute_macroscopic(f, rho, u_x, u_y, Nx, Ny);
    apply_inlet_bc(f, rho, u_x, u_y, Nx, Ny);
	collide(f, rho, u_x, u_y, tau, Nx, Ny);
	bounce_back(f, obstacle, Nx, Ny);
	stream(f, f_new, Nx, Ny);
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu
