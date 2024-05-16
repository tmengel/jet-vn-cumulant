#ifndef BKGDFUNCTIONS_H
#define BKGDFUNCTIONS_H

#include <vector>

namespace toymodel  
{

const double pi_vn_params[3][7] = {
    {-0.0220529, 0.172125, -0.0353618, -0.003559, 0.00113968, 4.0, 0.5},
    {-0.0066372, 0.0262161, 0.0372216, -0.0187145, 0.00228567, 4.0, 0.5},
    {0.000152642, 0.00135534, 0.0523496, -0.0225954, 0.0025451, 4.0, 0.5}
};

const double k_vn_params[3][7] = {
    {-0.0424241, 0.152629, -0.00506494, -0.0151633, 0.00254353, 4.0, 0.5},
    {-0.0149325, 0.0253627, 0.0329371, -0.0153877, 0.00170996, 4.0, 0.5},
    {-0.0171898, 0.0261749, 0.032913, -0.0180592, 0.00240376, 4.0, 0.5}
};

const double pro_vn_params[3][7] = { 
    {0.0128407, -0.0812974, 0.196424, -0.0729275, 0.0081403, 4.0, 0.5},
    {0.0216277, -0.0905268, 0.125852, -0.0410326, 0.00433817, 4.0, 0.5},
    {0.0296393, -0.113592, 0.137947, -0.0424535, 0.00422479, 4.0, 0.5}
};

const double blastwave_params[5][6] = { 
    {1.39570e-01, 1.39570e-01, 4.93680e-01, 4.93680e-01, 9.38270e-01, 9.38270e-01},
    {7.73497e-01, 7.72956e-01, 7.56607e-01, 7.54914e-01, 6.92553e-01, 6.61272e-01},
    {1.74370e-01, 1.73122e-01, 1.68712e-01, 1.64637e-01, 2.49546e-01, 2.40833e-01},
    {2.90054e+01, 2.84080e+01, 1.31964e+01, 1.24268e+01, 5.48703e+01, 3.47224e+01},
    {3.82884e+03, 4.03029e+03, 2.41347e+03, 2.66737e+03, 3.60574e+02, 3.41159e+02}
};

const std::vector<int> particle_ids = {211, -211, 321, -321, 2212, -2212};
const std::vector<double> particle_masses = {0.139570, 0.139570, 0.493677, 0.493677, 0.938272, 0.938272};
const std::vector<int> particle_yeilds = {145, 145,  22,  22,  16,  13};
const unsigned int n_species = 6;

double BWintegrand(double *x, double *p);
double MyStaticBGdNdPtTimesPt(double *x, double *p);

double BGBW(double *x, double *p);

double VnFunction(double *x, double *p);

double TruthRefVn(int harmonic, double min_particle_pt);

} // namespace toymodel

#endif