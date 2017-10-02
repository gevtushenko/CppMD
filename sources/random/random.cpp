//
// Created by egi on 10/2/17.
//

#include "random/random.h"

#include <cmath>

Random::Random() {
    m_iv[0] = 0;
    set_seed(0);
};

void Random::set_seed(int idum) {
    m_inc_prec = false;
    m_idum     = idum;
}

double Random::rand_u01() {
    if (m_inc_prec) {
        return u01d();
    }
    else {
        return u01();
    }
}

double Random::u01d () {
    double u = u01();

    u += u01() * m_fact;

    return (u < 1.0) ? u : (u - 1.0);
}

double Random::u01() {
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

double Random::gaussian() {
    double v1, v2, rsq;

    if(m_switch_gaussian) {
        m_switch_gaussian = false;
        return m_save_gaussian;
    }

    while(true) {
        v1 = 2.0 * rand_u01() - 1.0;
        v2 = 2.0 * rand_u01() - 1.0;

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

