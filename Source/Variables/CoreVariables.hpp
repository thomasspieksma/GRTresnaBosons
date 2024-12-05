#ifndef MULTIGRIDUSERVARIABLES_HPP
#define MULTIGRIDUSERVARIABLES_HPP

#include "ArrayTools.hpp"

// assign an enum to each variable

// Core vars
enum
{
    c_psi_reg,

    c_V1_0,
    c_V2_0,
    c_V3_0,
    c_U_0,

    c_K_0,

    c_A11_0,
    c_A12_0,
    c_A13_0,
    c_A22_0,
    c_A23_0,
    c_A33_0,

    NUM_METRIC_VARS
};

namespace MetricVariables
{
static const std::array<std::string, NUM_METRIC_VARS> variable_names = {
    "psi_reg",

    "V1_0",    "V2_0",  "V3_0",  "U_0",

    "K_0",

    "A11_0",   "A12_0", "A13_0", "A22_0", "A23_0", "A33_0"};

static constexpr std::array<int, NUM_METRIC_VARS> const vars_parity = {
    0, 1, 2, 3, 0, 0, 0, 4, 6, 0, 5, 0};

} // namespace MetricVariables

// assign an enum to each constraint variable
enum
{
    c_psi,

    c_V1,
    c_V2,
    c_V3,
    c_U,

    NUM_CONSTRAINT_VARS
};

namespace ConstraintUserVariables
{
static const std::array<std::string, NUM_CONSTRAINT_VARS> variable_names = {
    "psi", "V1", "V2", "V3", "U"};

static constexpr std::array<int, NUM_CONSTRAINT_VARS> const vars_parity = {
    0, 1, 2, 3, 0};

} // namespace ConstraintUserVariables

#endif /* MULTIGRIDUSERVARIABLES_HPP */