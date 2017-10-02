//
// Created by egi on 10/2/17.
//

#ifndef CPPMD_MD_SOLVER_H
#define CPPMD_MD_SOLVER_H

#include <valarray>
#include <vector>
#include <memory>

#include "random/random.h"

class Atom;

class MDSolver {
public:
    using Ptr = std::shared_ptr<MDSolver>;

    static inline Ptr create() {
        return std::make_shared<MDSolver>();
    }

    MDSolver();

    void read(std::string filename);
    void write(std::string filename);

    std::valarray<double> periodic_boundary_condition(std::valarray<double> cell, std::valarray<double> vec_in);

    // Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
    // it is a linear combination of old velocities and new, randomly chosen, velocity,
    // with proper coefficients
    void thermostat();

    void randomize_velocities();
    void calculate_velocities();
    void calculate_positions();
    void calculate_engkin();
    void calculate_forces();
    void check_need_new_lists();
    void renew_list();
    void solve();

private:
    std::vector< std::shared_ptr<Atom> > m_atoms;

    Random m_random;

    double m_temperature;
    double m_friction;
    double m_force_cutoff;
    double m_list_cutoff;
    double m_engint;
    double m_engkin;
    double m_dt;

    bool m_renew_list;

    std::size_t m_max_neighbors;
    int m_idum;
    bool m_wrap_atoms;

    std::valarray<double> m_cell;
};

#endif //CPPMD_MD_SOLVER_H
