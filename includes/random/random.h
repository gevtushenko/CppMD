//
// Created by egi on 10/2/17.
//

#ifndef CPPMD_RANDOM_H
#define CPPMD_RANDOM_H

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

    Random();

    void set_seed(int idum);

    double rand_u01();
    double u01d();
    double u01();

    double gaussian();

private:
    double m_save_gaussian;
    bool   m_switch_gaussian, m_inc_prec;

    int m_iv[m_ntab];
    int m_idum, m_iy;
};

#endif //CPPMD_RANDOM_H
