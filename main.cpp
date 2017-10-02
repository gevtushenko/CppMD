#include <iostream>
#include <valarray>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <cmath>

class Random {
public:
    static const int m_ia   = 16807,
                     m_im   = 2147483647,
                     m_iq   = 127773,
                     m_ir   = 2836,
                     m_ntab = 32;
    static const int m_ndiv = (1 + (m_im - 1) / m_ntab);

    const double m_fact = 5.9604644775390625e-8;     /* 1 / 2^24  */
    const double m_eps  = 3.0e-16;
    const double m_am   = 1.0/m_im;
    const double m_rnmx = (1.0 - m_eps);

    Random() {
        m_iv[0] = 0;
        set_seed(0);
    };

    void set_seed(int idum) {
        m_inc_prec = false;
        m_idum     = idum;
    }

    double RandU01 () {
        if (m_inc_prec) {
            return U01d();
        }
        else {
            return U01();
        }
    }

    double U01d () {
        double u = U01();

        u += U01() * m_fact;

        return (u < 1.0) ? u : (u - 1.0);
    }

    double U01() {
        int j,k;
        double temp;

        if (m_idum <= 0 || !m_iy) {
            if (-m_idum < 1) {
                m_idum=1;
            }
            else {
                m_idum = -m_idum;
            }

            for(j=m_ntab + 7; j >= 0; j--) {
                k      = m_idum/m_iq;
                m_idum = m_ia * (m_idum - k * m_iq) - m_ir * k;
                if (m_idum < 0) {
                    m_idum += m_im;
                }
                if (j < m_ntab) {
                    m_iv[j] = m_idum;
                }
            }
            m_iy = m_iv[0];
        }

        k      = m_idum / m_iq;
        m_idum = m_ia * (m_idum - k * m_iq) - m_ir * k;

        if (m_idum < 0) {
            m_idum += m_im;
        }

        j       = m_iy / m_ndiv;
        m_iy    = m_iv[j];
        m_iv[j] = m_idum;

        if ((temp = m_am * m_iy) > m_rnmx) {
            return m_rnmx;
        }

        else {
            return temp;
        }
    }

    double gaussian() {
        double v1, v2, rsq;

        if(m_switch_gaussian) {
            m_switch_gaussian = false;
            return m_save_gaussian;
        }

        while(true) {
            v1 = 2.0 * RandU01() - 1.0;
            v2 = 2.0 * RandU01() - 1.0;

            rsq = v1 * v1 + v2 * v2;

            if(rsq < 1.0 && rsq > 0.0) {
                break;
            }
        }

        double fac = std::sqrt(-2.0 * std::log(rsq) / rsq);

        m_save_gaussian = v1 * fac;
        m_switch_gaussian = true;

        return v2 * fac;
    }

private:
    double m_save_gaussian;
    bool   m_switch_gaussian, m_inc_prec;

    int m_iv[m_ntab];
    int m_idum, m_iy;
};

class Atom {
public:
    using Ptr = std::shared_ptr<Atom>;

    static inline Ptr create(std::valarray<double> position) {
        return std::make_shared<Atom>(position);
    }

    Atom(std::valarray<double> position)
        : m_mass(1.0)
    {
        m_position = position;
        m_velocity.resize(3, 0.0);
        m_forces.resize(3, 0.0);
    }

    double mass() const noexcept {
        return m_mass;
    }

    std::valarray<double> kinetic_energy() {
        return 0.5 * m_mass * m_velocity * m_velocity;
    }

    std::valarray<double> velocity() {
        return m_velocity;
    }

    std::valarray<double> position() {
        return m_position;
    }

    std::valarray<double> base_position() {
        return m_base_position;
    }

    void set_base_position(std::valarray<double> new_position) {
        m_base_position = new_position;
    }

    double update_velocity(double& c1, double& c2, Random& random) {
        m_velocity = c1 * m_velocity + c2 * random.gaussian();

        if(std::isnan(m_velocity[0])) {
            std::cerr << "Error! Nan velocity!" << std::endl;
        }
    }

    void calculate_velocity(double dt) {
        m_velocity += m_forces * 0.5 * dt / m_mass;

        if(std::isnan(m_velocity[0])) {
            std::cerr << "Error! Nan velocity! (forces: " << m_forces[0] << ", " << m_forces[1] << ", " << m_forces[2] << ")" << std::endl;
        }
    }

    void calculate_position(double dt) {
        m_position += m_velocity * dt;

        if(std::isnan(m_position[0])) {
            std::cerr << "Error! Nan position!" << std::endl;
        }
    }

    void clear_forces() {
        m_forces = 0.0;
    }

    double distance_to(Atom::Ptr atom) {
        auto dif    = atom->position() - position();
        auto dif_2  = dif * dif;
        return std::sqrt(dif_2.sum());
    }

    void clear_neighbors_list() {
        m_neighbors.clear();
    }

    void add_neighbor(Atom::Ptr neighbor) {
        m_neighbors.push_back(neighbor);
    }

    std::vector<Atom::Ptr>& neighbors() {
        return m_neighbors;
    }

    void randomize_velocity(double& temperature, Random& random) {
        m_velocity = std::sqrt(temperature / m_mass) * random.gaussian();
    }

    void plus_force(const std::valarray<double>& f) {
        m_forces += f;
    }

    void minus_force(const std::valarray<double>& f) {
        m_forces -= f;
    }

private:
    std::valarray<double> m_position;
    std::valarray<double> m_base_position;

    std::valarray<double> m_velocity;
    std::valarray<double> m_forces;

    std::vector<Ptr> m_neighbors;

    double m_mass;
};

class MDSolver {
public:
    using Ptr = std::shared_ptr<MDSolver>;

    static inline Ptr create() {
        return std::make_shared<MDSolver>();
    }

    MDSolver() {
        int idum = 0;
        m_random.set_seed(idum);

        m_dt = 0.005;
    }

    void read(std::string filename) {
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

    void write(std::string filename) {
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

    std::valarray<double> periodic_boundary_condition(std::valarray<double> cell, std::valarray<double> vec_in) {
        std::valarray<double> vec_out(std::size(vec_in));

        for(std::size_t i = 0; i < std::size(cell); ++i) {
            vec_out[i] = vec_in[i] - std::floor(vec_in[i] / cell[i] + 0.5) * cell[i];
        }

        return vec_out;
    }

    // Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
    // it is a linear combination of old velocities and new, randomly chosen, velocity,
    // with proper coefficients
    void thermostat() {
        double c1 = std::exp(-m_friction * m_dt);

        for(auto& atom: m_atoms) {
            double c2 = std::sqrt((1.0 - c1*c1) * m_temperature / atom->mass());

            m_engint += atom->kinetic_energy().sum();
                        atom->update_velocity(c1, c2, m_random);
            m_engint -= atom->kinetic_energy().sum();
        }
    }

    void randomize_velocities() {
        for(auto& atom: m_atoms) {
            atom->randomize_velocity(m_temperature, m_random);
        }
    }

    void calculate_velocities() {
        for(auto& atom: m_atoms) {
            atom->calculate_velocity(m_dt);
        }
    }

    void calculate_positions() {
        for(auto& atom: m_atoms) {
            atom->calculate_position(m_dt);
        }
    }

    void calculate_engkin() {
        m_engkin = 0.0;

        for(auto& atom: m_atoms) {
            m_engkin += atom->kinetic_energy().sum();
        }
    }

    void calculate_forces() {
        double engcorrection; // Energy necessary shift the potential avoiding discontinuities
        double force_cutoff_2 = m_force_cutoff * m_force_cutoff;

        for(auto& atom: m_atoms) {
            atom->clear_forces();
        }

        engcorrection = 4.0 * (1.0 / std::pow(force_cutoff_2, 6.0) - 1.0 / std::pow(force_cutoff_2, 3));

        double engconf = 0.0;

        for(auto& atom: m_atoms) {
            for(auto& neighbor_atom: atom->neighbors()) {
                auto distance = atom->position() - neighbor_atom->position();
                auto distance_pbc = periodic_boundary_condition(m_cell, distance);

                double distance_pbc_2 = (distance_pbc * distance_pbc).sum();

                if(distance_pbc_2 > force_cutoff_2) {
                    continue;
                }

                double distance_pbc_6  = std::pow(distance_pbc_2, 3);
                double distance_pbc_8  = distance_pbc_6  * distance_pbc_2;
                double distance_pbc_12 = distance_pbc_6  * distance_pbc_6;
                double distance_pbc_14 = distance_pbc_12 * distance_pbc_2;

                engconf += 4.0 * (1.0 / distance_pbc_12 - 1.0/distance_pbc_6) - engcorrection;

                auto f = 2.0 * distance_pbc * 4.0 * (6.0/distance_pbc_14 - 3.0/distance_pbc_8);

                if(distance_pbc_12 == 0.0) {
                    std::cerr << "Error! Zero distance on atom: " << atom->position()[0] << ", " << atom->position()[1] << ", " << atom->position()[2] << " and atom " << ", " << neighbor_atom->position()[0] << ", " << neighbor_atom->position()[1] << ", " << neighbor_atom->position()[2] << std::endl;
                    return;
                }

                atom->plus_force(f);
                neighbor_atom->minus_force(f);
            }
        }
    }

    void check_need_new_lists() {
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

    void renew_list() {
        double list_cutoff_2 = m_list_cutoff * m_list_cutoff;

        for(std::size_t aid = 0; aid < m_atoms.size(); ++aid) {
            auto atom = m_atoms[aid];

            atom->clear_neighbors_list();

            for(std::size_t naid = aid + 1; naid < m_atoms.size(); ++naid) {
                auto neighbor_atom = m_atoms[naid];

                auto distance = atom->position() - neighbor_atom->position();
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

    void solve() {
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

private:
    std::vector<Atom::Ptr> m_atoms;

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

int main() {
    auto max_step = 1000;
    auto solver   = MDSolver::create();
         solver->read("../examples/liquid.xyz");

    for(std::size_t step = 0; step < max_step; ++step) {
        std::cout << " Step: " << step << std::endl;
        solver->solve();
        solver->write("step_" + std::to_string(step) + ".vtk");
    }

    return 0;
}