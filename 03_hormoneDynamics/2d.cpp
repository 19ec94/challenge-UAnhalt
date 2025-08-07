#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <utility>
#include <algorithm>


using GridField2D = std::vector<std::vector<double>>;

enum class CellPhase { Omega1, Omega2, Omega3 };

std::ostream& operator<<(std::ostream& os, const CellPhase& phase) {
  switch (phase) {
    case CellPhase::Omega1: os << "Omega1"; break;
    case CellPhase::Omega2: os << "Omega2"; break;
    case CellPhase::Omega3: os << "Omega3"; break;
    default: os << "Unknown Phase"; break;
  }
  return os;
}


struct DomainParams {
  int Nc; // Number of cell cycles
  double Dc; // Duration of one cell cycle in x  (width in x)
  double Ys; // Maturity threshhold (0 < Ys < 1)
             //
  int cellsPerHalfCycle; // Numer of cells per half-cycle in x
                         //
  int Nx; // Total number of cells in x (computed)
  int Ny; // Total number of cells in y (choosen s.t Ys aligns with Grid)

  double Lx; // Total domain length in x (computed)
  double Ly=1.0; // Total domain length in y (fixed)

  double dx; // Grid spacing in x (computed)
  double dy; // Grid spacing in y (computed)

  DomainParams(int Nc_, double Dc_, double Ys_, int cellsPerHalfCycle_, int Ny_)
    : Nc(Nc_), Dc(Dc_), Ys(Ys_), cellsPerHalfCycle(cellsPerHalfCycle_), Ny(Ny_)
  {
    Nx = 2 * Nc * cellsPerHalfCycle;
    Lx = Nc * Dc;
    dx = Lx / Nx;
    dy = Ly / Ny;

    double j_s = Ys / dy;
    if (std::abs(j_s - std::round(j_s)) > 1e-12) {
      throw std::runtime_error("Ys must align exactly with vertical grid lines");
    }
  }
};


class Grid {
  public:
    const DomainParams& params;
    std::vector<double> xCoords; // edges in x, size Nx+1
    std::vector<double> yCoords; // edges in y, size Ny+1
    std::vector<double> xCenters; // centers in x, size Nx
    std::vector<double> yCenters; // centers in y, size Ny

    Grid(const DomainParams& p): params(p) {
      xCoords.resize(params.Nx + 1);
      yCoords.resize(params.Ny + 1);
      xCenters.resize(params.Nx);
      yCenters.resize(params.Ny);
      initialize();
    }

    void initialize() {
      for (int i=0; i <= params.Nx; ++i){
        xCoords[i] = i * params.dx;
      }
      for (int j=0; j <= params.Ny; ++j){
        yCoords[j] = j * params.dy;
      }
      for (int i=0; i < params.Nx; ++i){
        xCenters[i] = 0.5 * (xCoords[i] + xCoords[i+1]);
      }
      for (int j=0; j < params.Ny; ++j){
        yCenters[j] = 0.5 *(yCoords[j] + yCoords[j+1]);
      }
    }

    // Utility function: get cell index p for the horizontal cell cycle
    int cellCycleIndex(double x) const {
      int p = static_cast<int> (x / params.Dc);
      if (p < 0) p = 0;
      if (p >= params.Nc) p = params.Nc - 1;
      return p;
    }
};


class PhaseMap {
  public:
    const DomainParams& params;
    const Grid& grid;
    std::vector<std::vector<CellPhase>> phases; // size Nx x Ny

    PhaseMap(const DomainParams& p, const Grid& g) : params(p), grid(g) {
      phases.resize(params.Nx,
          std::vector<CellPhase> (params.Ny, CellPhase::Omega3));
      assignPhases();
    }

    void assignPhases() {
      for (int i = 0; i < params.Nx; ++i) {
        double x = grid.xCenters[i];
        int p = grid.cellCycleIndex(x);
        for (int j = 0; j < params.Ny; ++j){
          double y = grid.yCenters[j];
          if (y < params.Ys) {
            if (x >= p * params.Dc && x < p * params.Dc + 0.5 * params.Dc)
              phases[i][j] = CellPhase::Omega1;
            else if (x >= p * params.Dc + 0.5 * params.Dc && x < (p+1) * params.Dc)
              phases[i][j] = CellPhase::Omega2;
            else
              phases[i][j] = CellPhase::Omega3; // Outside subdomain ? Edge case
          }
          else {
            phases[i][j] = CellPhase::Omega3;
          }
        }
      }
    }

    CellPhase getPhase(int i, int j) const {
      return phases[i][j];
    }
};


enum class InterfaceType {
  None, 
  FluxContinuity, //Ω1 ↔ Ω2 interface with flux continuity
  FluxDoubling, //Ω2 ↔ Ω1 next cycle interface with flux doubling (mitosis) 
  Dirichlet //Ω2 ↔ Ω3 interface with Dirichlet condition (φ=0) 
};

std::ostream& operator<<(std::ostream& os, const InterfaceType& type) {
  switch (type) {
    case InterfaceType::None: os << "None"; break;
    case InterfaceType::FluxContinuity: os << "FluxContinuity"; break;
    case InterfaceType::FluxDoubling: os << "FluxDoubling"; break;
    case InterfaceType::Dirichlet: os << "Dirichlet"; break;
    default: os << "Unknown Interface Type"; break;
  }
  return os;
}


struct Interface {
  int i; // x-index of interface (between cell i-1 and i)
  int j; // y-index of interface (between cell j-1 and j)

  InterfaceType type;
  //Aditional info (could be orientation, phase boundaries, etc)
  //For example, verticle or horizontal interface
  enum Orientation {Vertical, Horizontal } orientation;

  Interface(int i_, int j_, InterfaceType t, Orientation o) 
    : i(i_), j(j_), type(t), orientation(o) {}
};

class InterfaceManager {
  const DomainParams& params;
  const Grid& grid;
  const PhaseMap& phaseMap;
  public:
  std::vector<Interface> verticalInterfaces;
  std::vector<Interface> horizontalInterfaces;

  InterfaceManager(const DomainParams& p, const Grid& g, const PhaseMap& ph)
    : params(p), grid(g), phaseMap(ph) {
      identifyInterfaces();
    }

  void identifyInterfaces() {
    // Identify vertical interfaces along x-edges between cells:
    // These exist at poisitons i where a change in phase or special
    // interface occurs.
    for (int i = 1; i < params.Nx; ++i){ //interface between cell i-1 and i
      for (int j = 0; j <params.Ny; ++j) {
        CellPhase leftPhase = phaseMap.getPhase(i-1, j);
        CellPhase rightPhase = phaseMap.getPhase(i,j);
        if (leftPhase != rightPhase) {
          // What to put here?
          InterfaceType type = InterfaceType::None;

          // Transmission condition: mitsosis end flux doubling (Omega2 to
          // Omega 1 of next cycle)
          if (leftPhase == CellPhase::Omega2 && rightPhase == CellPhase::Omega1) {
            // Transition end of cell cycle (mitosis)
            type = InterfaceType::FluxDoubling;
          } else if ((leftPhase == CellPhase::Omega1 && rightPhase == CellPhase::Omega2)
              || (leftPhase == CellPhase::Omega1 && rightPhase == CellPhase::Omega2)) {
            // Transition start of cell cycle phase change
            type = InterfaceType::FluxContinuity;
          }
          // What other conditions can be added?
          if (type != InterfaceType::None) {
            verticalInterfaces.emplace_back(i, j, type, Interface::Vertical);
          }
        }
      }
    }

    // Identify horizontal interfaces where y = ys (maturity threshhold)
    // This will be at j = j_s (integer) where grid line matches Ys
    int j_s = static_cast<int>(std::round(params.Ys / params.dy));
    for (int i = 0; i < params.Nx; ++i) {
      // Interface between cells j_s-1 and j_s (horizontal interface)
      // This interface separates proliferate and differential areas
      CellPhase bottomPhase = phaseMap.getPhase(i, j_s-1);
      CellPhase topPhase = phaseMap.getPhase(i, j_s);
      if (bottomPhase != topPhase) {
        horizontalInterfaces.emplace_back(i, j_s, InterfaceType::Dirichlet,
            Interface::Horizontal);
      }
    }
  }
};


// Biological constants (example values, adjust as per paper/table 5.1)
struct BioParams {
  double gamma1 = 2.0;
  double gamma2 = 2.0;
  double tau_h = 0.7;   // maturation velocity scale
  double c1 = 0.68;
  double c2 = 0.08;
  double u_bar = 0.02;  // basal FSH scale
  double Lambda_bar = 0.1;
  double gamma_bar = 0.01; // source term spatial scale
  double Ys = 0.3;          // maturity threshold
};


// Initialize phi with a Gaussian bump centered near the start of the domain in Omega1 region
void initializeFieldPhi(const DomainParams& params, const Grid& grid, const PhaseMap& phaseMap,
    std::vector<GridField2D>& phi, int follicleIndex) {

  double sigma = 0.05; //TODO:
  double cx = 0.25 * params.Dc;  // center x near start of cell cycle 0
  double cy = 0.5 * params.Ys;   // center y in proliferative zone

  for (int i = 0; i < params.Nx; ++i) {
    double x = grid.xCenters[i];
    int p = static_cast<int>(x / params.Dc); // get the cycle index
                                             //
    for (int j = 0; j < params.Ny; ++j) {
      double y = grid.yCenters[j];
      CellPhase phase = phaseMap.getPhase(i,j);

      // Only initialize Omega1 cell in the first cell cycle (p == 0)
      if (phase == CellPhase::Omega1 && p == 0) {
        double dx = x - cx;
        double dy = y - cy;
        phi[follicleIndex][i][j] = std::exp(-(dx*dx + dy*dy) / (2 * sigma * sigma));
      } else {
        phi[follicleIndex][i][j] = 0.0;
      }
    }
  }
}


// Initialize control hormones u_f and U to physiological baseline values
void initializeHormones(double& u_f, double& U) {
  u_f = 0.05;  // initial local FSH level (example)
  U = 0.5;     // initial plasma FSH level (example)
}


// Compute velocity fields g, h, and source term Lambda based on current phi, u_f, U
void computeVelocitiesAndSource(const DomainParams& params, const Grid& grid, const PhaseMap& phaseMap,
    double u_f, double U,
    std::vector<GridField2D>& g,
    std::vector<GridField2D>& h,
    std::vector<GridField2D>& Lambda,
    const BioParams& bioParams, int follicleIndex) {

  for (int i = 0; i < params.Nx; ++i) {
    for (int j = 0; j < params.Ny; ++j) {
      double y = grid.yCenters[j];
      CellPhase phase = phaseMap.getPhase(i, j);
      switch (phase) {
        case CellPhase::Omega1:
          g[follicleIndex][i][j] = bioParams.gamma1 * u_f + bioParams.gamma2;
          h[follicleIndex][i][j] = bioParams.tau_h * ( - y * y + (bioParams.c1 * y + bioParams.c2) ) * 
            (1.0 - std::exp(-u_f / bioParams.u_bar));
          Lambda[follicleIndex][i][j] = bioParams.Lambda_bar *
            std::exp(-pow(y - bioParams.Ys, 2) / bioParams.gamma_bar) * (1 - U);
          break;

        case CellPhase::Omega2:
          g[follicleIndex][i][j] = 1.0;
          h[follicleIndex][i][j] = 0.0;
          Lambda[follicleIndex][i][j] = 0.0;
          break;

        case CellPhase::Omega3:
          g[follicleIndex][i][j] = 1.0;
          h[follicleIndex][i][j] = bioParams.tau_h * (- y * y + (bioParams.c1 * y + bioParams.c2)) *
            (1.0 - std::exp(-u_f / bioParams.u_bar));
          Lambda[follicleIndex][i][j] = bioParams.Lambda_bar *
            std::exp(-pow(y - bioParams.Ys, 2) / bioParams.gamma_bar) * (1 - U);
          break;
      }
    }
  }
}

// Structure holding the biological parameters related to hormone feedback
struct HormoneParams {
  double U_min = 0.5;
  double c = 2.0;
  double M_ref = 4.5;
  double b1 = 0.08;
  double b2 = 2.25;
  double b3 = 1450.0;
};


// Compute maturity moment m_f(t) = \int y * phi dx dy
double computeMaturityMoment(const DomainParams& params,
    const Grid& grid,
    const std::vector<std::vector<double>>& phi) {
  double m = 0.0;
  for (int i = 0; i < params.Nx; ++i) {
    for (int j = 0; j < params.Ny; ++j) {
      // Approximate integral by midpoint rule
      double cell_area = params.dx * params.dy;
      m += grid.yCenters[j] * phi[i][j] * cell_area;
    }
  }
  return m;
}

// Compute global maturity M(t) = sum_f m_f(t)
double computeGlobalMaturity(const std::vector<double>& m_f_vector) {
  double M = 0.0;
  for(double m_f : m_f_vector)
    M += m_f;
  return M;
}

// Compute plasma FSH level U(t)
double computePlasmaFSH(double M, const HormoneParams& params) {
  return params.U_min + (1.0 - params.U_min) / (1.0 + std::exp(params.c * (M - params.M_ref)));
}

// Compute bioavailable local FSH u_f(t)
double computeLocalFSH(double m_f, double U, const HormoneParams& params) {
  double val = (params.b1 + std::exp(params.b2 * m_f)) / params.b3;
  if(val > 1.0) val = 1.0;
  return val * U;
}

int main() {
  DomainParams params(1, 0.1, 0.5, 1, 4);
  Grid grid(params);
  PhaseMap phaseMap(params, grid);
  InterfaceManager interfaceManager(params, grid, phaseMap);

  // Initialise phi, g, h, Lambda for F follicles
  int numberOfFollicles = 1;
  std::vector<GridField2D> phi(numberOfFollicles);
  std::vector<GridField2D> g(numberOfFollicles);
  std::vector<GridField2D> h(numberOfFollicles);
  std::vector<GridField2D> Lambda(numberOfFollicles);
  std::vector<double> u_f(numberOfFollicles, 0.0);
  std::vector<double> m_f(numberOfFollicles, 0.0);
  double U = 0.0;
  double M = 0.0;

  //Define biological parameters:
  BioParams bioParams;
  bioParams.Ys = params.Ys;
  HormoneParams hParams;

  // Resize fields once
   for (int f = 0; f < numberOfFollicles; ++f) {
    phi[f].resize(params.Nx, std::vector<double>(params.Ny, 0.0));
    g[f].resize(params.Nx, std::vector<double>(params.Ny, 0.0));
    h[f].resize(params.Nx, std::vector<double>(params.Ny, 0.0));
    Lambda[f].resize(params.Nx, std::vector<double>(params.Ny, 0.0));
  }

  //Step 1: For each follicle initialize phi and compute maturity moment
  for (int f = 0; f < numberOfFollicles; ++f) {
    initializeFieldPhi(params, grid, phaseMap, phi, f); 
    m_f[f] = computeMaturityMoment(params, grid, phi[f]);
  }

  // Step 2: Compute global Maturity and plasma FSH
  M = computeGlobalMaturity(m_f);
  U = computePlasmaFSH(M, hParams);

  // Step 3: Compute local bioavailable FSH for each follicle
   for (int f = 0; f < numberOfFollicles; ++f) {
    u_f[f] = computeLocalFSH(m_f[f], U, hParams);
  }

   // Step 4: Compute velocities and source terms for each follicle
  for (int f = 0; f < numberOfFollicles; ++f) {
    computeVelocitiesAndSource(params, grid, phaseMap, u_f[f], U,
                              g, h, Lambda, bioParams, f);
  }

  std::cout << "Domain Parameters:\n";
std::cout << "Nc = " << params.Nc << ", Dc = " << params.Dc << ", Ys = " << params.Ys << "\n";
std::cout << "Nx = " << params.Nx << ", Ny = " << params.Ny << "\n";
std::cout << "Lx = " << params.Lx << ", Ly = " << params.Ly << "\n";
std::cout << "dx = " << params.dx << ", dy = " << params.dy << "\n\n";

std::cout << "Grid centers x (first 10): ";
for (int i = 0; i < std::min(10, params.Nx); ++i)
    std::cout << grid.xCenters[i] << " ";
std::cout << "\nGrid centers y (all): ";
for (int j = 0; j < params.Ny; ++j)
    std::cout << grid.yCenters[j] << " ";
std::cout << "\n\n";

std::cout << "Phase map (y descending):\n";
for (int j = params.Ny - 1; j >= 0; --j) {
    for (int i = 0; i < params.Nx; ++i) {
        CellPhase phase = phaseMap.getPhase(i, j);
        char c = (phase == CellPhase::Omega1) ? '1' :
                 (phase == CellPhase::Omega2) ? '2' :
                 (phase == CellPhase::Omega3) ? '3' : '?';
        std::cout << c;
    }
    std::cout << "\n";
}
std::cout << "\n";
std::cout << "Vertical interfaces:\n";
for (const auto& iface : interfaceManager.verticalInterfaces) {
    double x_interface = grid.xCoords[iface.i];
    double y_center = grid.yCenters[iface.j];
    std::cout << "At x-edge between cells " << (iface.i - 1) << " and " << iface.i
              << " at x=" << x_interface << ", y=" << y_center
              << ", type=" << iface.type << "\n";
}
std::cout << "\n";

std::cout << "Horizontal interfaces:\n";
for (const auto& iface : interfaceManager.horizontalInterfaces) {
    double y_interface = grid.yCoords[iface.j];
    double x_center = grid.xCenters[iface.i];
    std::cout << "At y-edge between cells " << (iface.j - 1) << " and " << iface.j
              << " at y=" << y_interface << ", x=" << x_center
              << ", type=" << iface.type << "\n";
}
std::cout << "\n";

int f = 0;  // follicle index to print
std::cout << "Sample values for follicle " << f << ":\n";

// Print phi at some representative cells
std::cout << "phi:\n";
for (int j = params.Ny - 1; j >= 0; --j) {
    for (int i = 0; i < params.Nx; ++i) {
        std::cout << phi[f][i][j] << " ";
    }
    std::cout << "\n";
}
std::cout << "\n";

// Print velocities g,h and source Lambda
std::cout << "g:\n";
for (int j = params.Ny - 1; j >= 0; --j) {
    for (int i = 0; i < params.Nx; ++i) {
        std::cout << g[f][i][j] << " ";
    }
    std::cout << "\n";
}
std::cout << "\n";

std::cout << "h:\n";
for (int j = params.Ny - 1; j >= 0; --j) {
    for (int i = 0; i < params.Nx; ++i) {
        std::cout << h[f][i][j] << " ";
    }
    std::cout << "\n";
}
std::cout << "\n";

std::cout << "Lambda:\n";
for (int j = params.Ny - 1; j >= 0; --j) {
    for (int i = 0; i < params.Nx; ++i) {
        std::cout << Lambda[f][i][j] << " ";
    }
    std::cout << "\n";
}
std::cout << "\n";
for (int f = 0; f < numberOfFollicles; ++f) {
    std::cout << "Follicle " << f << " local FSH (u_f): " << u_f[f] << "\n";
}
std::cout << "Plasma FSH (U): " << U << "\n\n";

  return 0;
}
