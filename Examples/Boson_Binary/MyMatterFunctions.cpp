/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"
#include <fstream>
#include "SphericalHarmonics.hpp"
#include <cmath>
#include <complex>

Real ScalarField::my_potential_function(const Real &phi_here) const
{
    return 0.5 * pow(m_matter_params.scalar_mass * phi_here, 2.0);
}


Real ScalarField::my_phi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + (loc[1] - m_matter_params.offset_scalar) * (loc[1] - m_matter_params.offset_scalar) + loc[2] * loc[2]);

    const double spacing = 0.01; // in r for the values


    
    // // Interpolate data from read in values
    const int indxL = static_cast<int>(floor(rr / spacing));
    const int indxH = static_cast<int>(ceil(rr / spacing));

    // the field values
    double phi =
    m_matter_params.dPhi[indxL] +
        (rr / spacing - indxL) * (m_matter_params.dPhi[indxH] - m_matter_params.dPhi[indxL]);

    int ell_SH = 1;
    int m_SH = 1;
    int s_SH = 0;

    SphericalHarmonics::Y_lm_t<double> spi_Y_2Plusl_m = SphericalHarmonics::spin_Y_lm(loc[0], (loc[1] - m_matter_params.offset_scalar), loc[2], s_SH, 2.0 + ell_SH, m_SH);
    SphericalHarmonics::Y_lm_t<double> spi_Y_l_m = SphericalHarmonics::spin_Y_lm(loc[0], (loc[1] - m_matter_params.offset_scalar), loc[2], s_SH, ell_SH, m_SH);

    const complex<double> i(0.0,1.0);
    complex<double> spin_Y_2Plusl_m = spi_Y_2Plusl_m.Real + i * spi_Y_2Plusl_m.Im;
    complex<double> spin_Y_l_m = spi_Y_l_m.Real + i * spi_Y_l_m.Im;
    complex<double> nom_gamm_sq = i * 2.0 * sqrt(6) * spin_Y_2Plusl_m;
    complex<double> denom_gamm_sq = i * 50.0 * sqrt(21);
    
    complex<double> phi_ang = spin_Y_l_m+ m_matter_params.spheroidicity_param * nom_gamm_sq / denom_gamm_sq;

    double phi_ang_real = real(phi_ang);

    double phi_tot = phi* phi_ang_real;
    
    // double phi_tot = exp(-rr);// exp(-rr);

    // std::cout << "Interpolated phi: " << phi_tot << std::endl;

    return phi_tot;
}


Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    // Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    return 0.0;
}
