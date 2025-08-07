#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// Function declarations
double phi0(double x);
double velocity(double x);
double lambda(double x);
double godunov_flux(double g, double phiL, double phiR);
double lax_wendroff_flux(double g, double phiL, double phiR);
double calculate_limiter(const std::vector<double>& g, const std::vector<double>& phi, int k, int interface_index);
double koren_limiter(double r);
double limited_flux(const std::vector<double>& g, const std::vector<double>& phi, int k, int interface_index);
void apply_periodic_bc(std::vector<double>& phi);
void initialize_grid(std::vector<double>& x, std::vector<double>& g);
void initialize_phi(const std::vector<double>& x, std::vector<double>& phi);
void rk3_step(double dt, std::vector<double>& phi, const std::vector<double>& g,
    const std::vector<double>& x, const std::vector<double>& Lambda, int interface_index);
//double compute_exact_solution_testcase1(double x, double t, double gL, double gR, double xs);
double compute_exact_solution_testcase1(double x, double t, double gL, double gR, double xs, double lambda_l, double lambda_r);
double compute_L1_error(const std::vector<double>& phi,
                        const std::vector<double>& phi_exact,
                        double dx);

// Parameters
const int NX = 4000;
const double X_MIN = 0.0;
const double X_MAX = 1.0;
const int TESTCASE = 1;
double T_FINAL = 0.1;
double CFL = 0.1;
double G_L = 0.5;
double G_R = 1.0;
double LAMBDA_L = 0.0;
double LAMBDA_R = 0.0;

const double INTERFACE_X = 0.5;
int interface_index;

int main() {
  if (TESTCASE == 1) {
    T_FINAL = 0.4;
    CFL = 0.4;
    G_L = 0.5;
    G_R = 1.0;
    LAMBDA_L = 0.0;
    LAMBDA_R = 0.0;
  } else if (TESTCASE == 3) {
    T_FINAL = 0.0666;
    CFL = 0.02;
    G_L = 3.0;
    G_R = 1.0;
    LAMBDA_L = 1.0;
    LAMBDA_R = 0.0;
  }
  // Grid edges (Nx+1), velocities (Nx+1), and phi (Nx+2 with ghost cells)
  std::vector<double> x(NX + 1);
  std::vector<double> g(NX + 1);
  std::vector<double> phi(NX + 2, 0.0);  // phi[0] and phi[Nx+1] = ghosts
  std::vector<double> Lambda(NX + 2, 0.0); // No source, zero everywhere

  initialize_grid(x, g);
  initialize_phi(x, phi);

  apply_periodic_bc(phi);

  for (int i = 0; i <= NX + 1; ++i) {
    if (i == 0) {
      Lambda[i] = lambda(x[0] - (x[1] - x[0]) / 2.0);
    } else if (i == NX + 1) {
      Lambda[i] = lambda(x[NX] + (x[NX] - x[NX - 1]) / 2.0);
    } else {
      double xc = 0.5 * (x[i] + x[i - 1]);
      Lambda[i] = lambda(xc);
    }
  }

  double t = 0.0;
  double dt;
  double dx = x[1] - x[0];
  while (t < T_FINAL) {
    double max_g = 0.0;
    for (int k = 0; k <= NX; ++k) {
      if (std::abs(g[k]) > max_g) max_g = std::abs(g[k]);
    }
    dt = CFL * dx / max_g;
    if (t + dt > T_FINAL) dt = T_FINAL - t;

    rk3_step(dt, phi, g, x, Lambda, interface_index);
    t += dt;
  }

  // Compute and export exact solution for comparison
  std::vector<double> phi_exact(NX + 2, 0.0);
  for (int k = 1; k <= NX; ++k) {
    double xc = 0.5 * (x[k] + x[k - 1]);  // cell center
    phi_exact[k] = compute_exact_solution_testcase1(
        xc, T_FINAL, G_L, G_R, INTERFACE_X, LAMBDA_L, LAMBDA_R
    );
  }

  std::cout <<"The L1-error is " << compute_L1_error(phi, phi_exact, dx) << "\n";
  // Export results to file
  std::ofstream file("result.csv");
  file << "x,phi,phi_exact\n";
  for (int k = 1; k <= NX; ++k) {
    double xc = 0.5 * (x[k] + x[k - 1]);
    file << std::setprecision(8) << xc << "," << phi[k] << "," << phi_exact[k] << "\n";
  }
  file.close();

  std::cout << "Simulation finished. Data saved to result.csv\n";
  std::cout << "To visualise the data, please run python3 1d_plot.py\n";
  //std::cout << "Interface index at approx x = " << x[interface_index] << std::endl;

  return 0;
}

// Exact initial Gaussian condition phi0(x)
double phi0(double x) {
  const double c_x = 0.3;
  const double sigma = std::sqrt(0.002);
  double factor = 1.0 / (std::sqrt(2.0 * M_PI) * sigma);
  return factor * std::exp(-0.5 * ((x - c_x) * (x - c_x)) / (sigma * sigma));
}

double velocity(double x) {
  return (x < INTERFACE_X) ? G_L : G_R;
}

double lambda(double x) {
  return (x < INTERFACE_X) ? LAMBDA_L : LAMBDA_R;
}

// Godunov flux (low order)
double godunov_flux(double g, double phiL, double phiR) {
  return (g > 0) ? g * phiL : g * phiR;
}

// Lax-Wendroff flux (high order)
double lax_wendroff_flux(double g, double phiL, double phiR) {
  return 0.5 * g * (phiL + phiR);
}

// Calculate limiter ratio r_k with interface aware corrections
double calculate_limiter(const std::vector<double>& g, const std::vector<double>& phi, int k, int interface_index) {
  // Adjust calculation around interface to use continuous values (cf. equation (3.13))
  int N = (int)g.size();

  // Helper lambda to compute R ratio
  auto R = [](const std::vector<double>& z, const std::vector<double>& gs) -> double {
    bool all_positive = true, all_negative = true;
    for (double gi : gs) {
      if (gi < 0) all_positive = false;
      if (gi > 0) all_negative = false;
    }
    if (all_positive) {
      double denom = z[2] - z[1];
      if (std::abs(denom) < 1e-14) return 0.0;
      return (z[1] - z[0]) / denom;
    }
    if (all_negative) {
      double denom = z[1] - z[2];
      if (std::abs(denom) < 1e-14) return 0.0;
      return (z[3] - z[1]) / denom;
    }
    return 0.0;
  };

  // Normal periodic indexing with wrap-around
  auto idx = [&](int ind) {
    while (ind < 0) ind += N;
    while (ind >= N) ind -= N;
    return ind;
  };

  // Compute the ratio according to location w.r.t interface
  // Special treatment for k = interface_index-1, k=interface_index, k=interface_index+1
  if (k == interface_index - 1) {
    // r_{K-1}
    std::vector<double> z(4);
    std::vector<double> gs(4);
    z[0] = g[idx(k - 2)] * phi[idx(k - 2)];
    z[1] = g[idx(k - 1)] * phi[idx(k - 1)];
    z[2] = g[idx(k)] * phi[idx(k)];
    // Apply psi^{-1}_L to right side value (interface side)
    double phi_right = phi[idx(interface_index)];
    double g_right = g[idx(interface_index)];
    z[3] = phi_right * g_right; // For flux continuity psi_L = Id, so inverse is identity
    gs[0] = g[idx(k - 2)];
    gs[1] = g[idx(k - 1)];
    gs[2] = g[idx(k)];
    gs[3] = g_right;
    return koren_limiter(R(z, gs));
  }
  else if (k == interface_index) {
    // r_{K}
    std::vector<double> z(4);
    std::vector<double> gs(4);
    // For flux continuity psi_L=Id so no modification
    double phi_left_2 = phi[idx(interface_index - 2)];
    double phi_left_1 = phi[idx(interface_index - 1)];
    double phi_right = phi[idx(interface_index)];
    double phi_right_1 = phi[idx(interface_index + 1)];
    z[0] = phi_left_2 * g[idx(interface_index - 2)];
    z[1] = phi_left_1 * g[idx(interface_index - 1)];
    z[2] = phi_right * g[idx(interface_index)];
    z[3] = phi_right_1 * g[idx(interface_index + 1)];

    gs[0] = g[idx(interface_index - 2)];
    gs[1] = g[idx(interface_index - 1)];
    gs[2] = g[idx(interface_index)];
    gs[3] = g[idx(interface_index + 1)];
    return koren_limiter(R(z, gs));
  }
  else if (k == interface_index + 1) {
    // r_{K+1}
    std::vector<double> z(4);
    std::vector<double> gs(4);
    double phi_left_1 = phi[idx(interface_index - 1)];
    double phi_right = phi[idx(interface_index)];
    double phi_right_1 = phi[idx(interface_index + 1)];
    double phi_right_2 = phi[idx(interface_index + 2)];
    z[0] = phi_left_1 * g[idx(interface_index - 1)];
    z[1] = phi_right * g[idx(interface_index)];
    z[2] = phi_right_1 * g[idx(interface_index + 1)];
    z[3] = phi_right_2 * g[idx(interface_index + 2)];

    gs[0] = g[idx(interface_index - 1)];
    gs[1] = g[idx(interface_index)];
    gs[2] = g[idx(interface_index + 1)];
    gs[3] = g[idx(interface_index + 2)];
    return koren_limiter(R(z, gs));
  }
  else {
    // Usual case
    int km3 = idx(k - 3);
    int km2 = idx(k - 2);
    int km1 = idx(k - 1);
    int k_ = idx(k);
    double numerator = g[km1] * phi[km1] - g[km2] * phi[km2];
    double denominator = g[k_] * phi[k_] - g[km1] * phi[km1];
    double r = (std::fabs(denominator) > 1e-14) ? numerator / denominator : 0.0;
    return koren_limiter(r);
  }
} 


// Calculate Koren limiter l(r_k)
double koren_limiter(double r) {
  if (r <= 0) return 0.0;
  return std::fmax(0.0, std::fmin(std::fmin(2.0 * r, (2.0 + r) / 3.0), 2.0));
}

// Calculate limited flux F_k with special handling at interface
double limited_flux(const std::vector<double>& g, const std::vector<double>& phi, int k, int interface_index) {
  double limiter = calculate_limiter(g, phi, k, interface_index);
  double F_low = godunov_flux(g[k], phi[k - 1], phi[k]);
  double F_high = lax_wendroff_flux(g[k], phi[k - 1], phi[k]);
  return F_low + limiter * (F_high - F_low);
}

// Apply periodic boundary conditions using ghost cells (consistent with Aymard paper 3.1)
void apply_periodic_bc(std::vector<double>& phi) {
  phi[0] = phi[phi.size() - 2]; // Left ghost cell = last physical cell
  phi[phi.size() - 1] = phi[1]; // Right ghost cell = first physical cell
}

// Initialise grid cell edges and velocities at interfaces
void initialize_grid(std::vector<double>& x, std::vector<double>& g) {
  double dx = (X_MAX - X_MIN) / NX;
  for (int k = 0; k <= NX; ++k) {
    x[k] = X_MIN + k * dx; // cell edges
    g[k] = velocity(x[k]);
  }
  // Find interface index for grid edge closest to INTERFACE_X
  interface_index = 0;
  for (int k = 0; k <= NX; ++k) {
    if (std::abs(x[k] - INTERFACE_X) < dx * 0.5) {
      interface_index = k;
      break;
    }
  }
}

// Initialize phi with Gaussian, phi[k]: average in cell [x[k], x[k+1]]
void initialize_phi(const std::vector<double>& x, std::vector<double>& phi) {
  const double c_x = 0.3;
  const double sigma = std::sqrt(0.002);
  for (int k = 1; k <= NX; ++k) {
    // Cell center approx: midpoint between x[k-1], x[k]
    double xc = 0.5 * (x[k] + x[k - 1]);
    double exponent = std::exp(-0.5 * ((xc - c_x) * (xc - c_x)) / (sigma * sigma));
    double factor = 1.0 / (std::sqrt(2.0 * M_PI) * sigma);
    phi[k] = factor * exponent;
  }
}

// Single Runge-Kutta 3rd order time step with interface flux splitting
void rk3_step(double dt, std::vector<double>& phi, const std::vector<double>& g,
    const std::vector<double>& x, const std::vector<double>& Lambda, int interface_index) {
  int N = (int)phi.size();

  auto compute_rhs = [&](const std::vector<double>& phi_curr) {
    std::vector<double> F_L(N, 0.0);
    std::vector<double> F_R(N, 0.0);

    // Compute fluxes away from interface normally
    for (int k = 1; k < interface_index; ++k) {
      F_L[k] = limited_flux(g, phi_curr, k, interface_index);
      F_R[k] = F_L[k];
    }
    for (int k = interface_index + 1; k <= NX; ++k) {
      F_L[k] = limited_flux(g, phi_curr, k, interface_index);
      F_R[k] = F_L[k];
    }

    // At interface k = interface_index, compute left flux normally, apply transmission condition for right flux
    // For flux continuity psi_L = Id, so right flux = left flux
    // Use formulas for first and second order flux at interface from paper eq (3.11) and (3.12)

    // First order flux (Godunov) on left side of interface
    double g_left = g[interface_index - 1];
    double phi_left_1 = phi_curr[interface_index - 1];
    double g_right = g[interface_index];
    double phi_right = phi_curr[interface_index];

    double F_low_L = ((g_left > 0) ? g_left * phi_left_1 : 0.0) + ((g_right < 0) ? g_right * phi_right : 0.0);

    // Second order flux (Lax-Wendroff) on left side of interface
    double F_high_L = 0.5 * (g_left * phi_left_1 + g_right * phi_right);

    // Calculate limiter near interface using corrected calculation
    double limiter = calculate_limiter(g, phi_curr, interface_index, interface_index);

    // Combine to get high order flux at interface left side
    F_L[interface_index] = F_low_L + limiter * (F_high_L - F_low_L);

    // For flux continuity, right flux equals left flux
    F_R[interface_index] = F_L[interface_index];

    std::vector<double> rhs(N, 0.0);
    double dx = x[1] - x[0];

    // Update rhs using split fluxes: phi^{n+1}_k = phi^n_k - dt/dx (F_L_{k+1} - F_R_k) - dt Lambda_k phi_k
    for (int k = 1; k <= NX; ++k) {
      rhs[k] = -(F_L[k + 1] - F_R[k]) / dx - Lambda[k] * phi_curr[k];
    }
    return rhs;
  };

  std::vector<double> phi1 = phi;
  apply_periodic_bc(phi1);
  auto k1 = compute_rhs(phi1);
  for (int i = 1; i <= NX; ++i)
    phi1[i] = phi[i] + dt * k1[i];
  apply_periodic_bc(phi1);

  auto k2 = compute_rhs(phi1);
  std::vector<double> phi2 = phi;
  for (int i = 1; i <= NX; ++i)
    phi2[i] = 0.75 * phi[i] + 0.25 * (phi1[i] + dt * k2[i]);
  apply_periodic_bc(phi2);

  auto k3 = compute_rhs(phi2);
  for (int i = 1; i <= NX; ++i)
    phi[i] = (1.0 / 3.0) * phi[i] + (2.0 / 3.0) * (phi2[i] + dt * k3[i]);
  apply_periodic_bc(phi);
}

// Function to compute the exact solution phi_exact at position x and time t
double compute_exact_solution_testcase1(double x, double t, double gL, double gR, double xs, double lambda_l, double lambda_r) {
  // For flux continuity transmission: ψL = Id (identity)
  // Applies to test case 1: continuous flux, velocity discontinuity, no source.

  // Helper lambda for trace function Tr(t)
  auto Tr = [&](double tau) {
    double val = phi0(xs - gL * tau);
    return (gL / gR) * val;  // Flux continuity relation
  };

  if (x < xs) {
    // Left side: transported initial condition at speed gL
    return phi0(x - gL * t) * std::exp(-lambda_l * t);
  }
  else if (x < xs + gR * t) {
    // Right side near interface: use trace function at shifted time tau
    double tau = t - (x - xs) / gR;
    return Tr(tau) * std::exp(-lambda_l * tau);
  }
  else {
    // Far right side: transported initial condition at speed gR
    return phi0(x - gR * t) * std::exp(-lambda_r * t);
  }
}

// Computes the discrete L1 error norm between numerical and exact solutions
double compute_L1_error(const std::vector<double>& phi,
                        const std::vector<double>& phi_exact,
                        double dx) {
    double numerator = 0.0;
    double denominator = 0.0;
    size_t N = phi.size();

    for (size_t k = 0; k < N; ++k) {
        numerator += std::abs(phi[k] - phi_exact[k]) * dx;
        denominator += std::abs(phi_exact[k]) * dx;
    }

    // To avoid division by zero
    if (denominator == 0.0) {
        return numerator; // or return 0 or some indicator that exact solution is zero
    }

    return numerator / denominator;
}


