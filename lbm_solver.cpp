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
    
    // -------------------------------------------------------------
    // Create initial distribution function
    // -------------------------------------------------------------
    for (int q = 0; q < 9; q++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                f[q * Nx * Ny + i * Ny + j] = w[q] * rho0;
            }
        }
    }

	compute_macroscopic(f, rho, u_x, u_y, Nx, Ny);
	collide(f, rho, u_x, u_y, tau, Nx, Ny);
	stream(f, f_new, Nx, Ny);
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu
