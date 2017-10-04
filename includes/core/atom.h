//
// Created by egi on 10/2/17.
//

#ifndef CPPMD_ATOM_H
#define CPPMD_ATOM_H

#include <valarray>
#include <vector>
#include <memory>

class Random;

class Atom {
public:
    using Ptr = std::shared_ptr<Atom>;

    static inline Ptr create(std::valarray<double> position) {
        return std::make_shared<Atom>(position);
    }

    Atom(std::valarray<double> position);

    double mass() const noexcept;
    std::valarray<double> kinetic_energy();
    std::valarray<double> velocity();
    std::valarray<double> position();
    std::valarray<double> base_position();
    void set_base_position(std::valarray<double> new_position);
    double update_velocity(double& c1, double& c2, Random& random);
    void calculate_velocity(double dt);
    void calculate_position(double dt);
    void clear_forces();
    double distance_to(Atom::Ptr atom);
    void clear_neighbors_list();
    void add_neighbor(Atom::Ptr neighbor);
    std::vector< std::weak_ptr<Atom> >& neighbors();
    void randomize_velocity(double& temperature, Random& random);
    void plus_force(const std::valarray<double>& f);
    void minus_force(const std::valarray<double>& f);

private:
    std::valarray<double> m_position;
    std::valarray<double> m_base_position;

    std::valarray<double> m_velocity;
    std::valarray<double> m_forces;

    std::vector< std::weak_ptr<Atom> > m_neighbors;

    double m_mass;
};

#endif //CPPMD_ATOM_H
