/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include <fstream>

namespace MatterParams
{

struct params_t
{
    Real scalar_mass;
    Real spheroidicity_param;
    Real offset_scalar;
    Real *dPhi;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("scalar_mass", matter_params.scalar_mass);
    pp.get("spheroidicity_param", matter_params.spheroidicity_param);
    pp.get("offset_scalar", matter_params.offset_scalar);
    
    #define NUMRAD 2000
    static Real inputdPhi[NUMRAD];
    
    std::array<double, 1> tmp = {0.0};

    std::string phi_file("radial_profile_alpha02_short.csv");

    ifstream ifs0(phi_file);

    for (int i = 0; i < NUMRAD; ++i)
    {
        ifs0 >> tmp[0];

        inputdPhi[i]= tmp[0];
    }
    matter_params.dPhi = inputdPhi;

}

}; // namespace MatterParams

#endif
