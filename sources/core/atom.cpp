//
// Created by egi on 10/2/17.
//

#include <iostream>

#include "random/random.h"
#include "core/atom.h"

Atom::Atom(std::valarray<double> position)
        : m_mass(1.0)
{
    m_position = position;
    m_velocity.resize(3, 0.0);
    m_forces.resize(3, 0.0);
}

double Atom::mass() const noexcept {
    return m_mass;
}

std::valarray<double> Atom::kinetic_energy() {
    return 0.5 * m_mass * m_velocity * m_velocity;
}

std::valarray<double> Atom::velocity() {
    return m_velocity;
}

std::valarray<double> Atom::position() {
    return m_position;
}

std::valarray<double> Atom::base_position() {
    return m_base_position;
}

void Atom::set_base_position(std::valarray<double> new_position) {
    m_base_position = new_position;
}

double Atom::update_velocity(double& c1, double& c2, Random& random) {
    m_velocity = c1 * m_velocity + c2 * random.gaussian();

    if(std::isnan(m_velocity[0])) {
        std::cerr << "Error! Nan velocity!" << std::endl;
    }
}

void Atom::calculate_velocity(double dt) {
    m_velocity += m_forces * 0.5 * dt / m_mass;

    if(std::isnan(m_velocity[0])) {
        std::cerr << "Error! Nan velocity! (forces: " << m_forces[0] << ", " << m_forces[1] << ", " << m_forces[2] << ")" << std::endl;
    }
}

void Atom::calculate_position(double dt) {
    m_position += m_velocity * dt;

    if(std::isnan(m_position[0])) {
        std::cerr << "Error! Nan position!" << std::endl;
    }
}

void Atom::clear_forces() {
    m_forces = 0.0;
}

double Atom::distance_to(Atom::Ptr atom) {
    auto dif    = atom->position() - m_position;
    auto dif_2  = dif * dif;
    return std::sqrt(dif_2.sum());
}

void Atom::clear_neighbors_list() {
    m_neighbors.clear();
}

void Atom::add_neighbor(Atom::Ptr neighbor) {
    m_neighbors.push_back(neighbor);
}

std::vector<Atom::Ptr>& Atom::neighbors() {
    return m_neighbors;
}

void Atom::randomize_velocity(double& temperature, Random& random) {
    m_velocity = std::sqrt(temperature / m_mass) * random.gaussian();
}

void Atom::plus_force(const std::valarray<double>& f) {
    m_forces += f;
}

void Atom::minus_force(const std::valarray<double>& f) {
    m_forces -= f;
}
