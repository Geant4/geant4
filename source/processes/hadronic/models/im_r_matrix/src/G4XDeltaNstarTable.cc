//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4XDeltaNstarTable
//
//      Author:      Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)  
// 
//      Creation date: 4 June 1999
//
//      Modifications: 
//
// Hadron Kinetic Model
// p p -> Delta N* cross section tables
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4XDeltaNstarTable.hh"
#include "G4PhysicsFreeVector.hh"


const G4int G4XDeltaNstarTable::sizeDeltaNstar = 121;


// Energies (GeV) corresponding to the cross section table
// Units are assigned while filling the PhysicsVector

const G4double G4XDeltaNstarTable::energyTable[121] =
{
  0.0, 
  2.014, 2.014, 2.016, 2.018, 2.022, 2.026, 2.031, 2.037, 2.044, 2.052,  
  2.061, 2.071, 2.082, 2.094, 2.107, 2.121, 2.135, 2.151, 2.168, 2.185,  
  2.204, 2.223, 2.244, 2.265, 2.287, 2.311, 2.335, 2.360, 2.386, 2.413,  
  2.441, 2.470, 2.500, 2.531, 2.562, 2.595, 2.629, 2.664, 2.699, 2.736,
  2.773, 2.812, 2.851, 2.891, 2.933, 2.975, 3.018, 3.062, 3.107, 3.153,  
  3.200, 3.248, 3.297, 3.347, 3.397, 3.449, 3.502, 3.555, 3.610, 3.666,  
  3.722, 3.779, 3.838, 3.897, 3.957, 4.018, 4.081, 4.144, 4.208, 4.273,  
  4.339, 4.406, 4.473, 4.542, 4.612, 4.683, 4.754, 4.827, 4.900, 4.975,
  5.000, 6.134, 7.269, 8.403, 9.538, 10.672, 11.807, 12.941, 14.076, 15.210, 
 16.345, 17.479, 18.613, 19.748, 20.882, 22.017, 23.151, 24.286, 25.420, 26.555, 
 27.689, 28.824, 29.958, 31.092, 32.227, 33.361, 34.496, 35.630, 36.765, 37.899, 
 39.034, 40.168, 41.303, 42.437, 43.571, 44.706, 45.840, 46.975, 48.109, 49.244
};

// Cross-sections in mb, from S.A. Bass et al., Prog.Part.Nucl.Phys.41:225-370,1998 
// Units are assigned while filling the PhysicsVector

const G4double G4XDeltaNstarTable::sigmaDN1440[121] = 
{
  0.0, 
  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
  0.000, 0.000, 0.000, 0.001, 0.003, 0.006, 0.011, 0.021,
  0.041, 0.000, 0.000, 0.002, 0.011, 0.047, 0.131, 0.257,
  0.408, 0.568, 0.729, 0.886,  1.036,  1.178,  1.309,  1.431,
  1.542,  1.644,  1.736,  1.819,  1.894,  1.960,  2.017,  2.068,
  2.111,  2.148,  2.178,  2.203,  2.222,  2.237,  2.247,  2.253,
  2.254,  2.252,  2.247,  2.239,  2.228,  2.215,  2.199,  2.181,
  2.161,  2.140,  2.117,  2.093,  2.068,  2.042,  2.014,  1.986,
  1.977,  1.558,  1.224, 0.976, 0.797, 0.654, 0.548, 0.465,
  0.399, 0.346, 0.303, 0.268, 0.238, 0.213,  0.192,  0.173,
  0.158,  0.144,  0.132,  0.121,  0.112,  0.103,  0.096,  0.089,
  0.083,  0.078,  0.073,  0.069,  0.065,  0.061,  0.057,  0.054,
  0.051,  0.049,  0.046,  0.044,  0.042,  0.040,  0.038,  0.036 
};

const G4double G4XDeltaNstarTable::sigmaDN1520[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.001,  0.001,
  0.003,  0.006,  0.013,  0.027,  0.000,  0.001,  0.010,  0.048,
  0.152,  0.318,  0.513,  0.713,  0.908,  1.091,  1.262,  1.417,
  1.559,  1.686,  1.801,  1.903,  1.992,  2.071,  2.140,  2.198,
  2.248,  2.289,  2.323,  2.349,  2.369,  2.383,  2.392,  2.395,
  2.394,  2.389,  2.380,  2.368,  2.352,  2.334,  2.313,  2.290,
  2.266,  2.239,  2.211,  2.182,  2.151,  2.119,  2.087,  2.054,
  2.043,  1.571,  1.211,  0.953,  0.769,  0.627,  0.521,  0.440,
  0.376,  0.325,  0.284,  0.250,  0.222,  0.198,  0.178,  0.160,
  0.146,  0.133,  0.121,  0.112,  0.103,  0.095,  0.088,  0.082,
  0.076,  0.071,  0.067,  0.063,  0.059,  0.056,  0.052,  0.050,
  0.047,  0.045,  0.042,  0.040,  0.038,  0.036,  0.035,  0.033 
};

const G4double G4XDeltaNstarTable::sigmaDN1535[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.001,  0.001,  0.002,  0.004,  0.006,  0.010,
  0.014,  0.021,  0.030,  0.044,  0.000,  0.000,  0.002,  0.010,
  0.039,  0.097,  0.174,  0.257,  0.339,  0.417,  0.490,  0.556,
  0.615,  0.667,  0.714,  0.754,  0.790,  0.820,  0.846,  0.867,
  0.885,  0.899,  0.910,  0.918,  0.924,  0.927,  0.928,  0.927,
  0.924,  0.920,  0.914,  0.907,  0.899,  0.890,  0.880,  0.870,
  0.858,  0.847,  0.834,  0.822,  0.809,  0.795,  0.782,  0.768,
  0.764,  0.576,  0.439,  0.343,  0.275,  0.223,  0.185,  0.156,
  0.133,  0.115,  0.100,  0.088,  0.078,  0.070,  0.062,  0.056,
  0.051,  0.047,  0.043,  0.039,  0.036,  0.033,  0.031,  0.029,
  0.027,  0.025,  0.023,  0.022,  0.021,  0.019,  0.018,  0.017,
  0.016,  0.016,  0.015,  0.014,  0.013,  0.013,  0.012,  0.012 
};

const G4double G4XDeltaNstarTable::sigmaDN1650[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.001,  0.001,  0.002,  0.003,
  0.004,  0.006,  0.008,  0.011,  0.014,  0.019,  0.026,  0.000,
  0.000,  0.001,  0.004,  0.018,  0.047,  0.083,  0.122,  0.159,
  0.194,  0.227,  0.256,  0.282,  0.305,  0.326,  0.343,  0.359,
  0.372,  0.383,  0.393,  0.400,  0.406,  0.411,  0.415,  0.417,
  0.418,  0.418,  0.418,  0.417,  0.415,  0.412,  0.409,  0.405,
  0.401,  0.397,  0.392,  0.387,  0.382,  0.376,  0.371,  0.365,
  0.363,  0.278,  0.214,  0.168,  0.135,  0.110,  0.091,  0.077,
  0.066,  0.057,  0.049,  0.043,  0.038,  0.034,  0.031,  0.028,
  0.025,  0.023,  0.021,  0.019,  0.018,  0.016,  0.015,  0.014,
  0.013,  0.012,  0.012,  0.011,  0.010,  0.010,  0.009,  0.009,
  0.008,  0.008,  0.007,  0.007,  0.007,  0.006,  0.006,  0.006 
};

const G4double G4XDeltaNstarTable::sigmaDN1675[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.001,  0.001,  0.002,  0.004,  0.007,  0.013,
  0.000,  0.000,  0.003,  0.021,  0.072,  0.155,  0.251,  0.350,
  0.446,  0.538,  0.622,  0.699,  0.771,  0.835,  0.893,  0.945,
  0.992,  1.033,  1.068,  1.099,  1.126,  1.148,  1.167,  1.182,
  1.194,  1.203,  1.209,  1.212,  1.214,  1.213,  1.210,  1.205,
  1.199,  1.191,  1.182,  1.172,  1.161,  1.148,  1.135,  1.122,
  1.117,  0.891,  0.700,  0.556,  0.454,  0.371,  0.310,  0.262,
  0.225,  0.195,  0.170,  0.150,  0.133,  0.119,  0.107,  0.097,
  0.088,  0.080,  0.073,  0.067,  0.062,  0.057,  0.053,  0.050,
  0.046,  0.043,  0.040,  0.038,  0.036,  0.034,  0.032,  0.030,
  0.028,  0.027,  0.026,  0.024,  0.023,  0.022,  0.021,  0.020 
};

const G4double G4XDeltaNstarTable::sigmaDN1680[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.001,  0.001,  0.002,  0.004,  0.009,
  0.000,  0.000,  0.003,  0.019,  0.070,  0.157,  0.258,  0.361,
  0.461,  0.555,  0.641,  0.720,  0.792,  0.857,  0.915,  0.966,
  1.012,  1.051,  1.086,  1.116,  1.141,  1.162,  1.179,  1.192,
  1.202,  1.209,  1.214,  1.216,  1.215,  1.212,  1.208,  1.202,
  1.194,  1.185,  1.175,  1.163,  1.151,  1.137,  1.123,  1.108,
  1.103,  0.870,  0.678,  0.536,  0.437,  0.355,  0.296,  0.250,
  0.214,  0.185,  0.162,  0.142,  0.126,  0.113,  0.101,  0.092,
  0.083,  0.076,  0.069,  0.064,  0.059,  0.054,  0.050,  0.047,
  0.044,  0.041,  0.038,  0.036,  0.034,  0.032,  0.030,  0.028,
  0.027,  0.025,  0.024,  0.023,  0.022,  0.021,  0.020,  0.019
};

const G4double G4XDeltaNstarTable::sigmaDN1700[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.001,  0.001,  0.002,  0.004,  0.007,
  0.013,  0.000,  0.001,  0.006,  0.028,  0.078,  0.142,  0.209,
  0.273,  0.332,  0.385,  0.432,  0.475,  0.512,  0.545,  0.573,
  0.598,  0.620,  0.638,  0.653,  0.665,  0.675,  0.683,  0.688,
  0.692,  0.694,  0.695,  0.694,  0.692,  0.689,  0.685,  0.680,
  0.674,  0.668,  0.660,  0.653,  0.645,  0.636,  0.627,  0.618,
  0.615,  0.477,  0.368,  0.290,  0.234,  0.190,  0.158,  0.134,
  0.114,  0.099,  0.086,  0.076,  0.067,  0.060,  0.054,  0.049,
  0.044,  0.040,  0.037,  0.034,  0.031,  0.029,  0.027,  0.025,
  0.023,  0.022,  0.020,  0.019,  0.018,  0.017,  0.016,  0.015,
  0.014,  0.013,  0.013,  0.012,  0.012,  0.011,  0.010,  0.010 
};

const G4double G4XDeltaNstarTable::sigmaDN1710[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.001,  0.001,  0.001,  0.003,
  0.005,  0.000,  0.000,  0.002,  0.009,  0.028,  0.056,  0.086,
  0.116,  0.144,  0.170,  0.194,  0.215,  0.234,  0.251,  0.266,
  0.279,  0.290,  0.300,  0.308,  0.315,  0.321,  0.326,  0.329,
  0.332,  0.334,  0.335,  0.336,  0.335,  0.334,  0.333,  0.331,
  0.329,  0.327,  0.324,  0.320,  0.317,  0.313,  0.309,  0.305,
  0.304,  0.239,  0.187,  0.161,  0.121,  0.098,  0.082,  0.069,
  0.059,  0.051,  0.045,  0.040,  0.035,  0.031,  0.028,  0.025,
  0.023,  0.021,  0.019,  0.018,  0.016,  0.015,  0.014,  0.013,
  0.012,  0.011,  0.011,  0.010,  0.009,  0.009,  0.008,  0.008,
  0.007,  0.007,  0.007,  0.006,  0.006,  0.006,  0.006,  0.005 
};

const G4double G4XDeltaNstarTable::sigmaDN1720[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.001,  0.001,  0.002,  0.003,  0.005,  0.008,
  0.014,  0.000,  0.000,  0.001,  0.009,  0.033,  0.075,  0.124,
  0.175,  0.225,  0.271,  0.315,  0.355,  0.391,  0.423,  0.452,
  0.478,  0.501,  0.521,  0.538,  0.553,  0.566,  0.576,  0.585,
  0.592,  0.597,  0.601,  0.603,  0.604,  0.604,  0.603,  0.601,
  0.598,  0.595,  0.590,  0.585,  0.580,  0.574,  0.568,  0.561,
  0.559,  0.446,  0.351,  0.279,  0.228,  0.202,  0.156,  0.132,
  0.113,  0.098,  0.086,  0.076,  0.067,  0.060,  0.054,  0.049,
  0.044,  0.040,  0.037,  0.034,  0.031,  0.029,  0.027,  0.025,
  0.023,  0.022,  0.020,  0.019,  0.018,  0.017,  0.016,  0.015,
  0.014,  0.014,  0.013,  0.012,  0.012,  0.011,  0.011,  0.010 
};

const G4double G4XDeltaNstarTable::sigmaDN1900[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.001,  0.001,  0.002,  0.002,  0.003,
  0.005,  0.007,  0.009,  0.014,  0.021,  0.000,  0.000,  0.001,
  0.003,  0.009,  0.017,  0.028,  0.038,  0.049,  0.060,  0.070,
  0.080,  0.090,  0.099,  0.107,  0.115,  0.122,  0.128,  0.134,
  0.140,  0.145,  0.149,  0.153,  0.156,  0.159,  0.162,  0.164,
  0.166,  0.168,  0.169,  0.170,  0.170,  0.171,  0.171,  0.171,
  0.171,  0.155,  0.131,  0.110,  0.095,  0.079,  0.067,  0.058,
  0.051,  0.044,  0.039,  0.035,  0.031,  0.028,  0.026,  0.023,
  0.021,  0.019,  0.018,  0.016,  0.016,  0.015,  0.013,  0.012,
  0.011,  0.011,  0.010,  0.009,  0.009,  0.008,  0.008,  0.008,
  0.007,  0.007,  0.006,  0.006,  0.006,  0.006,  0.005,  0.005 
};

const G4double G4XDeltaNstarTable::sigmaDN1990[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.001,  0.001,  0.002,  0.004,  0.008,  0.016,
  0.000,  0.000,  0.001,  0.005,  0.013,  0.024,  0.035,  0.046,
  0.057,  0.068,  0.078,  0.088,  0.097,  0.105,  0.113,  0.120,
  0.126,  0.132,  0.138,  0.143,  0.147,  0.151,  0.155,  0.158,
  0.161,  0.164,  0.166,  0.168,  0.169,  0.171,  0.172,  0.172,
  0.173,  0.164,  0.143,  0.122,  0.105,  0.089,  0.076,  0.066,
  0.058,  0.053,  0.045,  0.040,  0.036,  0.032,  0.029,  0.027,
  0.024,  0.022,  0.021,  0.019,  0.018,  0.016,  0.015,  0.014,
  0.013,  0.012,  0.012,  0.011,  0.010,  0.010,  0.009,  0.009,
  0.008,  0.008,  0.007,  0.007,  0.007,  0.006,  0.006,  0.006 
};

const G4double G4XDeltaNstarTable::sigmaDN2090[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.001,  0.001,  0.002,  0.003,  0.004,  0.006,
  0.010,  0.000,  0.000,  0.001,  0.005,  0.013,  0.022,  0.033,
  0.043,  0.053,  0.063,  0.071,  0.079,  0.086,  0.093,  0.099,
  0.104,  0.108,  0.112,  0.116,  0.119,  0.121,  0.123,  0.125,
  0.126,  0.127,  0.128,  0.129,  0.129,  0.129,  0.129,  0.128,
  0.128,  0.110,  0.090,  0.073,  0.061,  0.049,  0.041,  0.035,
  0.030,  0.026,  0.023,  0.020,  0.018,  0.016,  0.015,  0.013,
  0.012,  0.011,  0.010,  0.009,  0.008,  0.008,  0.007,  0.007,
  0.006,  0.006,  0.006,  0.005,  0.005,  0.005,  0.004,  0.004,
  0.004,  0.004,  0.004,  0.003,  0.003,  0.003,  0.003,  0.003 
};

const G4double G4XDeltaNstarTable::sigmaDN2190[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.001,  0.001,
  0.002,  0.003,  0.006,  0.010,  0.000,  0.000,  0.000,  0.003,
  0.007,  0.013,  0.019,  0.026,  0.032,  0.038,  0.044,  0.050,
  0.055,  0.060,  0.065,  0.069,  0.073,  0.077,  0.080,  0.083,
  0.086,  0.088,  0.090,  0.092,  0.094,  0.095,  0.096,  0.097,
  0.098,  0.097,  0.086,  0.073,  0.062,  0.053,  0.046,  0.039,
  0.034,  0.030,  0.027,  0.024,  0.021,  0.019,  0.017,  0.016,
  0.015,  0.013,  0.012,  0.011,  0.010,  0.010,  0.009,  0.008,
  0.008,  0.007,  0.007,  0.006,  0.006,  0.006,  0.005,  0.005,
  0.005,  0.005,  0.004,  0.004,  0.004,  0.004,  0.004,  0.003 
};

const G4double G4XDeltaNstarTable::sigmaDN2220[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.001,  0.001,  0.003,  0.006,  0.011,  0.000,  0.000,  0.001,
  0.005,  0.010,  0.017,  0.024,  0.030,  0.037,  0.044,  0.050,
  0.055,  0.061,  0.066,  0.070,  0.074,  0.078,  0.082,  0.085,
  0.088,  0.091,  0.093,  0.095,  0.097,  0.099,  0.100,  0.101,
  0.102,  0.103,  0.092,  0.079,  0.067,  0.058,  0.050,  0.043,
  0.038,  0.033,  0.029,  0.027,  0.023,  0.021,  0.019,  0.017,
  0.016,  0.014,  0.013,  0.012,  0.011,  0.011,  0.010,  0.009,
  0.009,  0.008,  0.007,  0.007,  0.007,  0.006,  0.006,  0.006,
  0.005,  0.005,  0.005,  0.005,  0.004,  0.004,  0.004,  0.004, 
};

const G4double G4XDeltaNstarTable::sigmaDN2250[121] = 
{
  0.0,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
  0.001,  0.001,  0.002,  0.003,  0.006,  0.000,  0.000,  0.000,
  0.003,  0.007,  0.014,  0.021,  0.028,  0.035,  0.042,  0.049,
  0.055,  0.060,  0.066,  0.071,  0.075,  0.079,  0.083,  0.087,
  0.090,  0.093,  0.095,  0.098,  0.100,  0.101,  0.103,  0.104,
  0.105,  0.105,  0.093,  0.079,  0.067,  0.057,  0.049,  0.042,
  0.037,  0.032,  0.029,  0.025,  0.023,  0.020,  0.018,  0.017,
  0.015,  0.014,  0.013,  0.012,  0.011,  0.010,  0.009,  0.009,
  0.008,  0.008,  0.007,  0.007,  0.006,  0.006,  0.006,  0.005,
  0.005,  0.005,  0.005,  0.004,  0.004,  0.004,  0.004,  0.004 
};


G4XDeltaNstarTable::G4XDeltaNstarTable() 
{
  xMap["N(1440)0"] = (G4double*) sigmaDN1440;
  xMap["N(1440)+"] = (G4double*) sigmaDN1440;
  
  xMap["N(1520)0"] = (G4double*) sigmaDN1520;
  xMap["N(1520)+"] = (G4double*) sigmaDN1520;
  
  xMap["N(1535)0"] = (G4double*) sigmaDN1535;
  xMap["N(1535)+"] = (G4double*) sigmaDN1535;
  
  xMap["N(1650)0"] = (G4double*) sigmaDN1650;
  xMap["N(1650)+"] = (G4double*) sigmaDN1650;
  
  xMap["N(1675)0"] = (G4double*) sigmaDN1675;
  xMap["N(1675)+"] = (G4double*) sigmaDN1675;
  
  xMap["N(1680)0"] = (G4double*) sigmaDN1680;
  xMap["N(1680)+"] = (G4double*) sigmaDN1680;
  
  xMap["N(1700)0"] = (G4double*) sigmaDN1700;
  xMap["N(1700)+"] = (G4double*) sigmaDN1700;
  
  xMap["N(1710)0"] = (G4double*) sigmaDN1710;
  xMap["N(1710)+"] = (G4double*) sigmaDN1710;
  
  xMap["N(1720)0"] = (G4double*) sigmaDN1720;
  xMap["N(1720)+"] = (G4double*) sigmaDN1720;
  
  xMap["N(1900)0"] = (G4double*) sigmaDN1900;
  xMap["N(1900)+"] = (G4double*) sigmaDN1900;
    
  xMap["N(1990)0"] = (G4double*) sigmaDN1990;
  xMap["N(1990)+"] = (G4double*) sigmaDN1990;
  
  xMap["N(2090)0"] = (G4double*) sigmaDN2090;
  xMap["N(2090)+"] = (G4double*) sigmaDN2090;
  
  xMap["N(2190)0"] = (G4double*) sigmaDN2190;
  xMap["N(2190)+"] = (G4double*) sigmaDN2190;
  
  xMap["N(2220)0"] = (G4double*) sigmaDN2220;
  xMap["N(2220)+"] = (G4double*) sigmaDN2220;
  
  xMap["N(2250)0"] = (G4double*) sigmaDN2250;
  xMap["N(2250)+"] = (G4double*) sigmaDN2250; 
}


G4XDeltaNstarTable::~G4XDeltaNstarTable()
{ }


const G4PhysicsVector* G4XDeltaNstarTable::CrossSectionTable(const G4String& particleName) const
{
   // NOTE: the returned pointer is owned by the client

  if (xMap.find(particleName) != xMap.end())
    {
      // Cross section table for the requested particle available in the Map
      G4PhysicsFreeVector* sigmaVector = new G4PhysicsFreeVector(sizeDeltaNstar);
      std::map <G4String, G4double*, std::less<G4String> >::const_iterator iter;
      G4double* sigmaPointer = 0;
      for (iter = xMap.begin(); iter != xMap.end(); ++iter)
	{
	  G4String str = (*iter).first;
          if (str == particleName)
	    {
	      sigmaPointer = (*iter).second; 
	    }
	}

      G4int i;
      for (i=0; i<sizeDeltaNstar; i++)
	{
	  G4double value = *(sigmaPointer + i) * 0.5* millibarn;
	  G4double energy = energyTable[i] * GeV;
	  sigmaVector->PutValue(i,energy,value);
	}	    
      return sigmaVector;
    }
  else
    // No cross section table for the requested particle is available in the Map
    return 0;
}












