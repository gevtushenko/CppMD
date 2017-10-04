//
// Created by egi on 10/2/17.
//

#include <iostream>
#include <valarray>
#include <fstream>
#include <sstream>

#include "solvers/md_solver.h"
#include "core/atom.h"

MDSolver::MDSolver() {
    int idum = 0;
    m_random.set_seed(idum);

    m_dt = 0.005;
}

void MDSolver::read(std::string filename) {
    std::ifstream infile(filename);

    m_temperature   = 0.722;
    m_dt            = 0.005;
    m_friction      = 0.0;
    m_force_cutoff  = 2.5;
    m_list_cutoff   = 3.0;
    m_max_neighbors = 1000;
    m_idum          = 0;
    m_wrap_atoms    = false;

    std::string line;

    double cell_x, cell_y, cell_z;

    if(std::getline(infile, line)) {
        std::istringstream iss(line);

        if(!(iss >> cell_x >> cell_y >> cell_z)) {
            throw std::runtime_error("Error! Can't read file: '" + filename + "'!");
        }
    }

    m_cell = {cell_x, cell_y, cell_z};

    while(std::getline(infile, line)) {
        std::istringstream iss(line);

        double x, y, z;
        std::string name;

        if(!(iss >> name >> x >> y >> z)) {
            break;
        }
        m_atoms.push_back(Atom::create({x, y, z}));
    }

    std::cout << " Read " << m_atoms.size() << " atoms\n";

    randomize_velocities();
    renew_list();
    calculate_forces();
}

void MDSolver::write(std::string filename) {
    std::ofstream vtk(filename);

    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "vtk output\n";
    vtk << "ASCII\n";
    vtk << "DATASET UNSTRUCTURED_GRID\n";

    vtk << "POINTS " << m_atoms.size() << " double\n";

    for(auto& atom: m_atoms) {
        auto position = atom->position();

        vtk << position[0] << " " << position[1] << " " << position[2] << "\n";
    }

    vtk.close();
}

std::valarray<double> MDSolver::periodic_boundary_condition(std::valarray<double> cell, std::valarray<double> vec_in) {
    std::valarray<double> vec_out(std::size(vec_in));

    for(std::size_t i = 0; i < std::size(cell); ++i) {
        vec_out[i] = vec_in[i] - std::floor(vec_in[i] / cell[i] + 0.5) * cell[i];
    }

    return vec_out;
}

// Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
// it is a linear combination of old velocities and new, randomly chosen, velocity,
// with proper coefficients
void MDSolver::thermostat() {
    double c1 = std::exp(-m_friction * m_dt);

    for(auto& atom: m_atoms) {
        double c2 = std::sqrt((1.0 - c1*c1) * m_temperature / atom->mass());

        m_engint += atom->kinetic_energy().sum();
        atom->update_velocity(c1, c2, m_random);
        m_engint -= atom->kinetic_energy().sum();
    }
}

void MDSolver::randomize_velocities() {
    for(std::size_t aid = 0; aid < m_atoms.size(); ++aid) {
        m_atoms[aid]->randomize_velocity(m_temperature, m_random);
    }
}

void MDSolver::calculate_velocities() {
    for(std::size_t aid = 0; aid < m_atoms.size(); ++aid) {
        m_atoms[aid]->calculate_velocity(m_dt);
    }
}

void MDSolver::calculate_positions() {
    for(std::size_t aid = 0; aid < m_atoms.size(); ++aid) {
        m_atoms[aid]->calculate_position(m_dt);
    }
}

void MDSolver::calculate_engkin() {
    m_engkin = 0.0;

    for(auto& atom: m_atoms) {
        m_engkin += atom->kinetic_energy().sum();
    }
}

void MDSolver::calculate_forces() {
    double engcorrection; // Energy necessary shift the potential avoiding discontinuities
    double force_cutoff_2 = m_force_cutoff * m_force_cutoff;

    for(std::size_t aid = 0; aid < m_atoms.size(); ++aid) {
        m_atoms[aid]->clear_forces();
    }

    engcorrection = 4.0 * (1.0 / std::pow(force_cutoff_2, 6.0) - 1.0 / std::pow(force_cutoff_2, 3));

    // double engconf = 0.0;

#pragma omp parallel for
    for(std::size_t aid = 0; aid < m_atoms.size() - 1; ++aid) {
        auto atom = m_atoms[aid];
        auto neighbors = atom->neighbors();

        for(std::size_t nid = 0; nid < neighbors.size(); ++nid) {
            auto neighbor_atom = neighbors[nid].lock();

            std::valarray<double> distance = atom->position() - neighbor_atom->position();
            std::valarray<double> distance_pbc = periodic_boundary_condition(m_cell, distance);

            double distance_pbc_2 = (distance_pbc * distance_pbc).sum();

            if(distance_pbc_2 > force_cutoff_2) {
                continue;
            }

            double distance_pbc_6  = std::pow(distance_pbc_2, 3);
            double distance_pbc_8  = distance_pbc_6  * distance_pbc_2;
            double distance_pbc_12 = distance_pbc_6  * distance_pbc_6;
            double distance_pbc_14 = distance_pbc_12 * distance_pbc_2;

            // engconf += 4.0 * (1.0 / distance_pbc_12 - 1.0/distance_pbc_6) - engcorrection;

            std::valarray<double> f = 2.0 * distance_pbc * 4.0 * (6.0/distance_pbc_14 - 3.0/distance_pbc_8);

            atom->plus_force(f);
            neighbor_atom->minus_force(f);
        }
    }
}

void MDSolver::check_need_new_lists() {
    std::valarray<double> displacement;
    double delta_2 = std::pow(0.5 * (m_list_cutoff - m_force_cutoff), 2);
    m_renew_list = false;

    for(auto& atom: m_atoms) {
        displacement = atom->position() - atom->base_position();

        if((displacement * displacement).sum() > delta_2) {
            m_renew_list = true; break;
        }
    }
}

void MDSolver::renew_list() {
    double list_cutoff_2 = m_list_cutoff * m_list_cutoff;

    for(std::size_t aid = 0; aid < m_atoms.size(); ++aid) {
        auto atom = m_atoms[aid];

        atom->clear_neighbors_list();

        for(std::size_t naid = aid + 1; naid < m_atoms.size(); ++naid) {
            auto neighbor_atom = m_atoms[naid];

            std::valarray<double> distance = atom->position() - neighbor_atom->position();
            auto distance_pbc = periodic_boundary_condition(m_cell, distance);

            double d_2 = (distance_pbc * distance_pbc).sum();

            if(d_2 > list_cutoff_2) {
                continue;
            }

            atom->add_neighbor(neighbor_atom);
        }

        atom->set_base_position(atom->position());
    }
}

void MDSolver::solve() {
    thermostat();

    calculate_velocities();
    calculate_positions();

    check_need_new_lists();
    if(m_renew_list) {
        renew_list();
    }

    calculate_forces();
    calculate_velocities();

    thermostat();

    calculate_engkin();
}

double MDSolver::dt() const noexcept {
    return m_dt;
}
