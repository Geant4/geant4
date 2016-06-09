//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4RDGenerator2BN
//
// Author:        Andreia Trindade (andreia@lip.pt)
//                Pedro Rodrigues  (psilva@lip.pt)
//                Luis Peralta     (luis@lip.pt)
//
// Creation date: 21 June 2003
//
// Modifications: 
// 21 Jun 2003                                 First implementation acording with new design
// 05 Nov 2003  MGP                            Fixed compilation warning
// 14 Mar 2004                                 Code optimization
//
// Class Description: 
//
// Concrete base class for Bremsstrahlung Angular Distribution Generation - 2BN Distribution
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//    

#include "G4RDGenerator2BN.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
//    

G4double G4RDGenerator2BN::Atab[320] =
         { 0.0264767, 0.0260006, 0.0257112, 0.0252475, 0.024792, 0.0243443, 0.0240726, 0.0236367,
           0.023209,  0.0227892, 0.0225362, 0.0221284, 0.0217282,0.0214894, 0.0211005, 0.0207189,
           0.0204935, 0.0201227, 0.0197588, 0.019546,  0.0191923,0.0188454, 0.0186445, 0.0183072,
           0.0179763, 0.0177866, 0.0174649, 0.0172828, 0.0169702,0.0166634, 0.0164915, 0.0161933,
           0.0160283, 0.0157384, 0.0155801, 0.0152981, 0.0151463,0.0148721, 0.0147263, 0.0144598,
           0.01432,   0.0140607, 0.0139267, 0.0136744, 0.0135459,0.0133005, 0.0131773, 0.0129385,
           0.0128205, 0.0125881, 0.012475,  0.0122488, 0.0121406,0.0119204, 0.0118167, 0.0117158,
           0.0115032, 0.0114067, 0.0111995, 0.0111072, 0.0110175,0.0108173, 0.0107316, 0.0105365,
           0.0104547, 0.0102646, 0.0101865, 0.010111,  0.00992684,0.0098548,0.00967532,0.00960671,
           0.00943171,0.00936643,0.00930328,0.0091337, 0.00907372,0.00890831,0.00885141,0.00869003,
           0.00863611,0.00858428,0.00842757,0.00837854,0.0082256,0.00817931,0.00803,0.00798639,
           0.00784058,0.00779958,0.00776046,0.00761866,0.00758201,0.00744346,0.00740928,0.00727384,
           0.00724201,0.00710969,0.00708004,0.00695074,0.00692333,0.00679688,0.00677166,0.00664801,
           0.00662484,0.00650396,0.00648286,0.00636458,0.00634545,0.00622977,0.00621258,0.00609936,
           0.00608412,0.00597331,0.00595991,0.00585143,0.00583988,0.0057337,0.0057239,0.00561991,
           0.0056119, 0.00551005,0.00550377,0.00540399,0.00539938,0.00530162,0.00529872,0.00520292,
           0.0051091, 0.00510777,0.00501582,0.00501608,0.00492594,0.00492781,0.00483942,0.0048429,
           0.00475622,0.00476127,0.00467625,0.00468287,0.00459947,0.00451785,0.00452581,0.00444573,
           0.00445522,0.00437664,0.00438768,0.00431057,0.00432316,0.00424745,0.0042616,0.00418726,
           0.004203,  0.00413,   0.00405869,0.00407563,0.00400561,0.00402414,0.00395536,0.00397553,
           0.00390795,0.00392975,0.00386339,0.00379862,0.00382167,0.00375805,0.00378276,0.00372031,
           0.00374678,0.00368538,0.00371363,0.00365335,0.00359463,0.0036242, 0.00356653,0.003598,
           0.00354139,0.00357481,0.00351921,0.00355464,0.00350005,0.0034471,0.00348403,0.00343208,
           0.0034712, 0.00342026,0.00346165,0.00341172,0.00345548,0.00340657,0.00335944,0.00340491,
           0.00335885,0.00340692,0.00336191,0.00341273,0.00336879,0.00342249,0.00337962,0.00333889,
           0.00339463,0.00335506,0.00341401,0.00337558,0.00343797,0.00340067,0.00336584,0.00343059,
           0.0033969, 0.00346557,0.00343302,0.00350594,0.00347448,0.00344563,0.00352163,0.00349383,
           0.00357485,0.00354807,0.00352395,0.00360885,0.00358571,0.00367661,0.00365446,0.00375194,
           0.00373078,0.00371234,0.00381532,0.00379787,0.00390882,0.00389241,0.00387881,0.00399675,
           0.00398425,0.00411183,0.00410042,0.00409197,0.00422843,0.00422123,0.00436974,0.0043637,
           0.00436082,0.00452075,0.00451934,0.00452125,0.00469406,0.00469756,0.00488741,0.00489221,
           0.00490102,0.00510782,0.00511801,0.00513271,0.0053589, 0.00537524,0.00562643,0.00564452,
           0.0056677, 0.00594482,0.00596999,0.0059999, 0.00630758,0.00634014,0.00637849,0.00672136,
           0.00676236,0.00680914,0.00719407,0.0072439, 0.00730063,0.0077349, 0.00779494,0.00786293,
           0.00835577,0.0084276, 0.00850759,0.00907162,0.00915592,0.00924925,0.00935226,0.00999779,
           0.0101059, 0.0102249, 0.0109763, 0.0111003, 0.0112367, 0.0113862, 0.0122637, 0.0124196,
           0.0125898, 0.0136311, 0.0138081, 0.0140011, 0.0142112, 0.0154536, 0.0156723, 0.0159099,
           0.016168,  0.0176664, 0.0179339, 0.0182246, 0.0185396, 0.020381,  0.0207026, 0.0210558,
           0.0214374, 0.0237377, 0.0241275, 0.0245528, 0.0250106, 0.0255038, 0.0284158, 0.0289213,
           0.0294621, 0.0300526, 0.0338619, 0.0344537, 0.0351108, 0.0358099, 0.036554,  0.0416399
         };

G4double G4RDGenerator2BN::ctab[320] =
    { 0.482253, 0.482253, 0.489021, 0.489021, 0.489021, 0.489021, 0.495933,
      0.495933, 0.495933, 0.495933, 0.502993, 0.502993, 0.502993, 0.510204,
      0.510204, 0.510204, 0.517572, 0.517572, 0.517572, 0.5251,   0.5251,
      0.5251,   0.532793, 0.532793, 0.532793, 0.540657, 0.540657, 0.548697,
      0.548697, 0.548697, 0.556917, 0.556917, 0.565323, 0.565323, 0.573921,
      0.573921, 0.582717, 0.582717, 0.591716, 0.591716, 0.600925, 0.600925,
      0.610352, 0.610352, 0.620001, 0.620001, 0.629882, 0.629882, 0.64,
      0.64,     0.650364, 0.650364, 0.660982, 0.660982, 0.671862, 0.683013,
      0.683013, 0.694444, 0.694444, 0.706165, 0.718184, 0.718184, 0.730514,
      0.730514, 0.743163, 0.743163, 0.756144, 0.769468, 0.769468, 0.783147,
      0.783147, 0.797194, 0.797194, 0.811622, 0.826446, 0.826446, 0.84168,
      0.84168,  0.857339, 0.857339, 0.873439, 0.889996, 0.889996, 0.907029,
      0.907029, 0.924556, 0.924556, 0.942596, 0.942596, 0.961169, 0.980296,
      0.980296, 1,        1,        1.0203,   1.0203,   1.04123,  1.04123,
      1.06281,  1.06281,  1.08507,  1.08507,  1.10803,  1.10803,  1.13173,
      1.13173,  1.1562,   1.1562,   1.18147,  1.18147,  1.20758,  1.20758,
      1.23457,  1.23457,  1.26247,  1.26247,  1.29132,  1.29132,  1.32118,
      1.32118,  1.35208,  1.35208,  1.38408,  1.38408,  1.41723,  1.41723,
      1.45159,  1.45159,  1.45159,  1.48721,  1.48721,  1.52416,  1.52416,
      1.5625,   1.5625,   1.60231,  1.60231,  1.64366,  1.64366,  1.68663,
      1.68663,  1.68663,  1.7313,   1.7313,   1.77778,  1.77778,  1.82615,
      1.82615,  1.87652,  1.87652,  1.92901,  1.92901,  1.98373,  1.98373,
      1.98373,  2.04082,  2.04082,  2.1004,   2.1004,   2.16263,  2.16263,
      2.22767,  2.22767,  2.22767,  2.29568,  2.29568,  2.36686,  2.36686,
      2.44141,  2.44141,  2.51953,  2.51953,  2.51953,  2.60146,  2.60146,
      2.68745,  2.68745,  2.77778,  2.77778,  2.87274,  2.87274,  2.87274,
      2.97265,  2.97265,  3.07787,  3.07787,  3.18878,  3.18878,  3.30579,
      3.30579,  3.30579,  3.42936,  3.42936,  3.55999,  3.55999,  3.69822,
      3.69822,  3.84468,  3.84468,  3.84468,  4,        4,        4.16493,
      4.16493,  4.34028,  4.34028,  4.34028,  4.52694,  4.52694,  4.7259,
      4.7259,   4.93827,  4.93827,  4.93827,  5.16529,  5.16529,  5.40833,
      5.40833,  5.40833,  5.66893,  5.66893,  5.94884,  5.94884,  6.25,
      6.25,     6.25,     6.57462,  6.57462,  6.92521,  6.92521,  6.92521,
      7.3046,   7.3046,   7.71605,  7.71605,  7.71605,  8.16327,  8.16327,
      8.65052,  8.65052,  8.65052,  9.18274,  9.18274,  9.18274,  9.76562,
      9.76562,  10.4058,  10.4058,  10.4058,  11.1111,  11.1111,  11.1111,
      11.8906,  11.8906,  12.7551,  12.7551,  12.7551,  13.7174,  13.7174,
      13.7174,  14.7929,  14.7929,  14.7929,  16,       16,       16,
      17.3611,  17.3611,  17.3611,  18.9036,  18.9036,  18.9036,  20.6612,
      20.6612,  20.6612,  22.6757,  22.6757,  22.6757,  22.6757,  25,
      25,       25,       27.7008,  27.7008,  27.7008,  27.7008,  30.8642,
      30.8642,  30.8642,  34.6021,  34.6021,  34.6021,  34.6021,  39.0625,
      39.0625,  39.0625,  39.0625,  44.4444,  44.4444,  44.4444,  44.4444,
      51.0204,  51.0204,  51.0204,  51.0204,  59.1716,  59.1716,  59.1716,
      59.1716,  59.1716,  69.4444,  69.4444,  69.4444,  69.4444,  82.6446,
      82.6446,  82.6446,  82.6446,  82.6446,  100
    };


G4RDGenerator2BN::G4RDGenerator2BN(const G4String& name):G4RDVBremAngularDistribution(name)
{
  b = 1.2;
  index_min = -300;
  index_max = 20;

  // Set parameters minimum limits Ekelectron = 250 eV and kphoton = 100 eV
  kmin = 100*eV;
  Ekmin = 250*eV;
  kcut = 100*eV;

  // Increment Theta value for surface interpolation
  dtheta = 0.1*rad;

  // Construct Majorant Surface to 2BN cross-section
  // (we dont need this. Pre-calculated values are used instead due to performance issues
  // ConstructMajorantSurface();
}

//    

G4RDGenerator2BN::~G4RDGenerator2BN() 
{;}

//

G4double G4RDGenerator2BN::PolarAngle(const G4double initial_energy,
				    const G4double final_energy,
				    const G4int) // Z
{

  G4double theta = 0;

  G4double k = initial_energy - final_energy;

  theta = Generate2BN(initial_energy, k);

  return theta;
}
//

G4double G4RDGenerator2BN::CalculateFkt(G4double k, G4double theta, G4double A, G4double c) const
{
  G4double Fkt_value = 0;
  Fkt_value = A*std::pow(k,-b)*theta/(1+c*theta*theta);
  return Fkt_value;
}

G4double G4RDGenerator2BN::Calculatedsdkdt(G4double kout, G4double theta, G4double Eel) const
{

  G4double dsdkdt_value = 0.;
  G4double Z = 1;
  // classic radius (in cm)
  G4double r0 = 2.82E-13;
  //squared classic radius (in barn)
  G4double r02 = r0*r0*1.0E+24;

  // Photon energy cannot be greater than electron kinetic energy
  if(kout > (Eel-electron_mass_c2)){
    dsdkdt_value = 0;
    return dsdkdt_value;
  }

     G4double E0 = Eel/electron_mass_c2;
     G4double k =  kout/electron_mass_c2;
     G4double E =  E0 - k;

     if(E <= 1*MeV ){                                 // Kinematic limit at 1 MeV !!
       dsdkdt_value = 0;
       return dsdkdt_value;
     }


     G4double p0 = std::sqrt(E0*E0-1);
     G4double p = std::sqrt(E*E-1);
     G4double L = std::log((E*E0-1+p*p0)/(E*E0-1-p*p0));
     G4double delta0 = E0 - p0*std::cos(theta);
     G4double epsilon = std::log((E+p)/(E-p));
     G4double Z2 = Z*Z;
     G4double sintheta2 = std::sin(theta)*std::sin(theta);
     G4double E02 = E0*E0;
     G4double E2 = E*E;
     G4double p02 = E0*E0-1;
     G4double k2 = k*k;
     G4double delta02 = delta0*delta0;
     G4double delta04 =  delta02* delta02;
     G4double Q = std::sqrt(p02+k2-2*k*p0*std::cos(theta));
     G4double Q2 = Q*Q;
     G4double epsilonQ = std::log((Q+p)/(Q-p));


     dsdkdt_value = Z2 * (r02/(8*pi*137)) * (1/k) * (p/p0) *
       ( (8 * (sintheta2*(2*E02+1))/(p02*delta04)) -
         ((2*(5*E02+2*E*E0+3))/(p02 * delta02)) -
         ((2*(p02-k2))/((Q2*delta02))) +
         ((4*E)/(p02*delta0)) +
         (L/(p*p0))*(
                 ((4*E0*sintheta2*(3*k-p02*E))/(p02*delta04)) +
                 ((4*E02*(E02+E2))/(p02*delta02)) +
                 ((2-2*(7*E02-3*E*E0+E2))/(p02*delta02)) +
                 (2*k*(E02+E*E0-1))/((p02*delta0))
                 ) -
         ((4*epsilon)/(p*delta0)) +
         ((epsilonQ)/(p*Q))*
         (4/delta02-(6*k/delta0)-(2*k*(p02-k2))/(Q2*delta0))
         );



     dsdkdt_value = dsdkdt_value*std::sin(theta);
     return dsdkdt_value;
}

void G4RDGenerator2BN::ConstructMajorantSurface()
{

  G4double Eel;
  G4int vmax;
  G4double Ek;
  G4double k, theta, k0, theta0, thetamax;
  G4double ds, df, dsmax, dk, dt;
  G4double ratmin;
  G4double ratio = 0;
  G4double Vds, Vdf;
  G4double A, c;

  G4cout << "**** Constructing Majorant Surface for 2BN Distribution ****" << G4endl;

  if(kcut > kmin) kmin = kcut;

  G4int i = 0;
  // Loop on energy index
  for(G4int index = index_min; index < index_max; index++){

  G4double fraction = index/100.;
  Ek = std::pow(10.,fraction);
  Eel = Ek + electron_mass_c2;

  // find x-section maximum at k=kmin
  dsmax = 0.;
  thetamax = 0.;

  for(theta = 0.; theta < pi; theta = theta + dtheta){

    ds = Calculatedsdkdt(kmin, theta, Eel);

    if(ds > dsmax){
      dsmax = ds;
      thetamax = theta;
    }
  }

  // Compute surface paremeters at kmin
  if(Ek < kmin || thetamax == 0){
    c = 0;
    A = 0;
  }else{
    c = 1/(thetamax*thetamax);
    A = 2*std::sqrt(c)*dsmax/(std::pow(kmin,-b));
  }

  // look for correction factor to normalization at kmin 
  ratmin = 1.;

  // Volume under surfaces
  Vds = 0.;
  Vdf = 0.;
  k0 = 0.;
  theta0 = 0.;

  vmax = G4int(100.*std::log10(Ek/kmin));

  for(G4int v = 0; v < vmax; v++){
    G4double fraction = (v/100.);
    k = std::pow(10.,fraction)*kmin;

    for(theta = 0.; theta < pi; theta = theta + dtheta){
      dk = k - k0;
      dt = theta - theta0;
      ds = Calculatedsdkdt(k,theta, Eel);
      Vds = Vds + ds*dk*dt;
      df = CalculateFkt(k,theta, A, c);
      Vdf = Vdf + df*dk*dt;

      if(ds != 0.){
	if(df != 0.) ratio = df/ds;
      }

      if(ratio < ratmin && ratio != 0.){
	ratmin = ratio;
      }
    }
  }


  // sampling efficiency
  Atab[i] = A/ratmin * 1.04;
  ctab[i] = c;

//  G4cout << Ek << " " << i << " " << index << " " << Atab[i] << " " << ctab[i] << " " << G4endl;
  i++;
  }

}

G4double G4RDGenerator2BN::Generate2BN(G4double Ek, G4double k) const 
{

  G4double Eel;
  G4double t;
  G4double cte2;
  G4double y, u;
  G4double fk, ft;
  G4double ds;
  G4double A2;
  G4double A, c;

  G4int trials = 0;
  G4int index;

  // find table index
  index = G4int(std::log10(Ek)*100) - index_min;
  Eel = Ek + electron_mass_c2;

  c = ctab[index];
  A = Atab[index];
  if(index < index_max){
    A2 = Atab[index+1];
    if(A2 > A) A = A2;
  }

  do{
  // generate k accordimg to std::pow(k,-b)
  trials++;

  // normalization constant 
//  cte1 = (1-b)/(std::pow(kmax,(1-b))-std::pow(kmin2,(1-b)));
//  y = G4UniformRand();
//  k = std::pow(((1-b)/cte1*y+std::pow(kmin2,(1-b))),(1/(1-b)));

  // generate theta accordimg to theta/(1+c*std::pow(theta,2)
  // Normalization constant
  cte2 = 2*c/std::log(1+c*pi2);

  y = G4UniformRand();
  t = std::sqrt((std::exp(2*c*y/cte2)-1)/c);
  u = G4UniformRand();

  // point acceptance algorithm
  fk = std::pow(k,-b);
  ft = t/(1+c*t*t);
  ds = Calculatedsdkdt(k,t,Eel);

  // violation point
  if(ds > (A*fk*ft)) G4cout << "WARNING IN G4RDGenerator2BN !!!" << Ek << " " <<  (ds-A*fk*ft)/ds << G4endl;

  }while(u*A*fk*ft > ds);

  return t;

}

void G4RDGenerator2BN::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Bremsstrahlung Angular Generator is 2BN Generator from 2BN Koch & Motz distribution (Rev Mod Phys 31(4), 920 (1959))" << G4endl;
  G4cout << "\n" << G4endl;
} 

