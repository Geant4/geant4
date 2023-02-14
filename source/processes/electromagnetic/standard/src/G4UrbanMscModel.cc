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
// -------------------------------------------------------------------
//   
// GEANT4 Class file
//    
//
// File name:   G4UrbanMscModel
//
// Author:      Laszlo Urban
//
// Creation date: 19.02.2013
//
// Created from G4UrbanMscModel96
//
// New parametrization for theta0
// Correction for very small step length
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526 and others

// -------------------------------------------------------------------
// In its present form the model can be  used for simulation 
//   of the e-/e+ multiple scattering


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UrbanMscModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Positron.hh"
#include "G4EmParameters.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4ProductionCutsTable.hh"

#include "G4Poisson.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4AutoLock.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4UrbanMscModel::mscData*> G4UrbanMscModel::msc;

namespace
{
  G4Mutex theUrbanMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UrbanMscModel::G4UrbanMscModel(const G4String& nam)
  : G4VMscModel(nam)
{
  masslimite    = 0.6*CLHEP::MeV;
  fr            = 0.02;
  taubig        = 8.0;
  tausmall      = 1.e-16;
  taulim        = 1.e-6;
  currentTau    = taulim;
  tlimitminfix  = 0.01*CLHEP::nm;             
  tlimitminfix2 = 1.*CLHEP::nm;
  stepmin       = tlimitminfix;
  smallstep     = 1.e10;
  currentRange  = 0. ;
  rangeinit     = 0.;
  tlimit        = 1.e10*CLHEP::mm;
  tlimitmin     = 10.*tlimitminfix;            
  tgeom         = 1.e50*CLHEP::mm;
  geombig       = tgeom;
  geommin       = 1.e-3*CLHEP::mm;
  geomlimit     = geombig;
  presafety     = 0.*CLHEP::mm;

  particle      = nullptr;

  positron      = G4Positron::Positron();
  rndmEngineMod = G4Random::getTheEngine();

  firstStep     = true; 
  insideskin    = false;
  latDisplasmentbackup = false;
  dispAlg96 = true;

  drr      = 0.35;
  finalr   = 10.*CLHEP::um;

  tlow = 5.*CLHEP::keV;
  invmev = 1.0/CLHEP::MeV;

  skindepth = skin*stepmin;

  mass = CLHEP::proton_mass_c2;
  charge = chargeSquare = 1.0;
  currentKinEnergy = currentRadLength = lambda0 = lambdaeff = tPathLength 
    = zPathLength = par1 = par2 = par3 = rndmarray[0] = rndmarray[1] = 0;
  currentLogKinEnergy = LOG_EKIN_MIN;

  idx = 0;
  fParticleChange = nullptr;
  couple = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UrbanMscModel::~G4UrbanMscModel()
{
  if(isFirstInstance) {
    for(auto & ptr : msc) { delete ptr; }
    msc.clear();
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UrbanMscModel::Initialise(const G4ParticleDefinition* p,
                                 const G4DataVector&)
{
  // set values of some data members
  SetParticle(p);
  fParticleChange = GetParticleChangeForMSC(p);
  InitialiseParameters(p);

  latDisplasmentbackup = latDisplasment;
  dispAlg96 = G4EmParameters::Instance()->LateralDisplacementAlg96();
  fPosiCorrection = G4EmParameters::Instance()->MscPositronCorrection();

  // initialise cache only once
  if(0 == msc.size()) {
    G4AutoLock l(&theUrbanMutex);
    if(0 == msc.size()) {
      isFirstInstance = true;
      msc.resize(1, nullptr);
    }
    l.unlock();
  }
  // initialise cache for each new run
  if(isFirstInstance) { InitialiseModelCache(); }

  /*
  G4cout << "### G4UrbanMscModel::Initialise done for " 
 	 << p->GetParticleName() << " type= " << steppingAlgorithm << G4endl;
  G4cout << "    RangeFact= " << facrange << " GeomFact= " << facgeom
	 << " SafetyFact= " << facsafety << " LambdaLim= " << lambdalimit
	 << G4endl;
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UrbanMscModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* part,
                                   G4double kinEnergy,
                                   G4double atomicNumber,G4double,
                                   G4double, G4double)
{
  static const G4double epsmin = 1.e-4 , epsmax = 1.e10;

  static const G4double Zdat[15] = { 4.,  6., 13., 20., 26., 29., 32., 38.,47.,
                                     50., 56., 64., 74., 79., 82. };

  // corr. factors for e-/e+ lambda for T <= Tlim
  static const G4double celectron[15][22] =
          {{1.125,1.072,1.051,1.047,1.047,1.050,1.052,1.054,
            1.054,1.057,1.062,1.069,1.075,1.090,1.105,1.111,
            1.112,1.108,1.100,1.093,1.089,1.087            },
           {1.408,1.246,1.143,1.096,1.077,1.059,1.053,1.051,
            1.052,1.053,1.058,1.065,1.072,1.087,1.101,1.108,
            1.109,1.105,1.097,1.090,1.086,1.082            },
           {2.833,2.268,1.861,1.612,1.486,1.309,1.204,1.156,
            1.136,1.114,1.106,1.106,1.109,1.119,1.129,1.132,
            1.131,1.124,1.113,1.104,1.099,1.098            },
           {3.879,3.016,2.380,2.007,1.818,1.535,1.340,1.236,
            1.190,1.133,1.107,1.099,1.098,1.103,1.110,1.113,
            1.112,1.105,1.096,1.089,1.085,1.098            },
           {6.937,4.330,2.886,2.256,1.987,1.628,1.395,1.265,
            1.203,1.122,1.080,1.065,1.061,1.063,1.070,1.073,
            1.073,1.070,1.064,1.059,1.056,1.056            },
           {9.616,5.708,3.424,2.551,2.204,1.762,1.485,1.330,
            1.256,1.155,1.099,1.077,1.070,1.068,1.072,1.074,
            1.074,1.070,1.063,1.059,1.056,1.052            },
           {11.72,6.364,3.811,2.806,2.401,1.884,1.564,1.386,
            1.300,1.180,1.112,1.082,1.073,1.066,1.068,1.069,
            1.068,1.064,1.059,1.054,1.051,1.050            },
           {18.08,8.601,4.569,3.183,2.662,2.025,1.646,1.439,
            1.339,1.195,1.108,1.068,1.053,1.040,1.039,1.039,
            1.039,1.037,1.034,1.031,1.030,1.036            },
           {18.22,10.48,5.333,3.713,3.115,2.367,1.898,1.631,
            1.498,1.301,1.171,1.105,1.077,1.048,1.036,1.033,
            1.031,1.028,1.024,1.022,1.021,1.024            },
           {14.14,10.65,5.710,3.929,3.266,2.453,1.951,1.669,
            1.528,1.319,1.178,1.106,1.075,1.040,1.027,1.022,
            1.020,1.017,1.015,1.013,1.013,1.020            },
           {14.11,11.73,6.312,4.240,3.478,2.566,2.022,1.720,
            1.569,1.342,1.186,1.102,1.065,1.022,1.003,0.997,
            0.995,0.993,0.993,0.993,0.993,1.011            },
           {22.76,20.01,8.835,5.287,4.144,2.901,2.219,1.855,
            1.677,1.410,1.224,1.121,1.073,1.014,0.986,0.976,
            0.974,0.972,0.973,0.974,0.975,0.987            },
           {50.77,40.85,14.13,7.184,5.284,3.435,2.520,2.059,
            1.837,1.512,1.283,1.153,1.091,1.010,0.969,0.954,
            0.950,0.947,0.949,0.952,0.954,0.963            },
           {65.87,59.06,15.87,7.570,5.567,3.650,2.682,2.182,
            1.939,1.579,1.325,1.178,1.108,1.014,0.965,0.947,
            0.941,0.938,0.940,0.944,0.946,0.954            },
           {55.60,47.34,15.92,7.810,5.755,3.767,2.760,2.239,
            1.985,1.609,1.343,1.188,1.113,1.013,0.960,0.939,
            0.933,0.930,0.933,0.936,0.939,0.949            }};
            
  static const G4double cpositron[15][22] = {
           {2.589,2.044,1.658,1.446,1.347,1.217,1.144,1.110,
            1.097,1.083,1.080,1.086,1.092,1.108,1.123,1.131,
            1.131,1.126,1.117,1.108,1.103,1.100            },
           {3.904,2.794,2.079,1.710,1.543,1.325,1.202,1.145,
            1.122,1.096,1.089,1.092,1.098,1.114,1.130,1.137,
            1.138,1.132,1.122,1.113,1.108,1.102            },
           {7.970,6.080,4.442,3.398,2.872,2.127,1.672,1.451,
            1.357,1.246,1.194,1.179,1.178,1.188,1.201,1.205,
            1.203,1.190,1.173,1.159,1.151,1.145            },
           {9.714,7.607,5.747,4.493,3.815,2.777,2.079,1.715,
            1.553,1.353,1.253,1.219,1.211,1.214,1.225,1.228,
            1.225,1.210,1.191,1.175,1.166,1.174            },
           {17.97,12.95,8.628,6.065,4.849,3.222,2.275,1.820,
            1.624,1.382,1.259,1.214,1.202,1.202,1.214,1.219,
            1.217,1.203,1.184,1.169,1.160,1.151            },
           {24.83,17.06,10.84,7.355,5.767,3.707,2.546,1.996,
            1.759,1.465,1.311,1.252,1.234,1.228,1.238,1.241,
            1.237,1.222,1.201,1.184,1.174,1.159            },
           {23.26,17.15,11.52,8.049,6.375,4.114,2.792,2.155,
            1.880,1.535,1.353,1.281,1.258,1.247,1.254,1.256,
            1.252,1.234,1.212,1.194,1.183,1.170            },
           {22.33,18.01,12.86,9.212,7.336,4.702,3.117,2.348,
            2.015,1.602,1.385,1.297,1.268,1.251,1.256,1.258,
            1.254,1.237,1.214,1.195,1.185,1.179            },
           {33.91,24.13,15.71,10.80,8.507,5.467,3.692,2.808,
            2.407,1.873,1.564,1.425,1.374,1.330,1.324,1.320,
            1.312,1.288,1.258,1.235,1.221,1.205            },
           {32.14,24.11,16.30,11.40,9.015,5.782,3.868,2.917,
            2.490,1.925,1.596,1.447,1.391,1.342,1.332,1.327,
            1.320,1.294,1.264,1.240,1.226,1.214            },
           {29.51,24.07,17.19,12.28,9.766,6.238,4.112,3.066,
            2.602,1.995,1.641,1.477,1.414,1.356,1.342,1.336,
            1.328,1.302,1.270,1.245,1.231,1.233            },
           {38.19,30.85,21.76,15.35,12.07,7.521,4.812,3.498,
            2.926,2.188,1.763,1.563,1.484,1.405,1.382,1.371,
            1.361,1.330,1.294,1.267,1.251,1.239            },
           {49.71,39.80,27.96,19.63,15.36,9.407,5.863,4.155,
            3.417,2.478,1.944,1.692,1.589,1.480,1.441,1.423,
            1.409,1.372,1.330,1.298,1.280,1.258            },
           {59.25,45.08,30.36,20.83,16.15,9.834,6.166,4.407,
            3.641,2.648,2.064,1.779,1.661,1.531,1.482,1.459,
            1.442,1.400,1.354,1.319,1.299,1.272            },
           {56.38,44.29,30.50,21.18,16.51,10.11,6.354,4.542,
            3.752,2.724,2.116,1.817,1.692,1.554,1.499,1.474,
            1.456,1.412,1.364,1.328,1.307,1.282            }};

  //data/corrections for T > Tlim  
                                             
  static const G4double hecorr[15] = {
    120.70, 117.50, 105.00, 92.92, 79.23,  74.510,  68.29,
    57.39,  41.97,  36.14, 24.53, 10.21,  -7.855, -16.84,
    -22.30};

  G4double sigma;
  SetParticle(part);

  G4double Z23 = G4Pow::GetInstance()->Z23(G4lrint(atomicNumber));

  // correction if particle .ne. e-/e+
  // compute equivalent kinetic energy
  // lambda depends on p*beta ....

  G4double eKineticEnergy = kinEnergy;

  if(mass > CLHEP::electron_mass_c2)
  {
     G4double TAU = kinEnergy/mass ;
     G4double c = mass*TAU*(TAU+2.)/(CLHEP::electron_mass_c2*(TAU+1.)) ;
     G4double w = c-2.;
     G4double tau = 0.5*(w+std::sqrt(w*w+4.*c)) ;
     eKineticEnergy = CLHEP::electron_mass_c2*tau ;
  }

  G4double eTotalEnergy = eKineticEnergy + CLHEP::electron_mass_c2 ;
  G4double beta2 = eKineticEnergy*(eTotalEnergy+CLHEP::electron_mass_c2)
                                 /(eTotalEnergy*eTotalEnergy);
  G4double bg2   = eKineticEnergy*(eTotalEnergy+CLHEP::electron_mass_c2)
                                 /(CLHEP::electron_mass_c2*CLHEP::electron_mass_c2);

  static const G4double epsfactor = 2.*CLHEP::electron_mass_c2*
    CLHEP::electron_mass_c2*CLHEP::Bohr_radius*CLHEP::Bohr_radius
    /(CLHEP::hbarc*CLHEP::hbarc);
  G4double eps = epsfactor*bg2/Z23;

  if     (eps<epsmin)  sigma = 2.*eps*eps;
  else if(eps<epsmax)  sigma = G4Log(1.+2.*eps)-2.*eps/(1.+2.*eps);
  else                 sigma = G4Log(2.*eps)-1.+1./eps;

  sigma *= chargeSquare*atomicNumber*atomicNumber/(beta2*bg2);

  // interpolate in AtomicNumber and beta2 
  G4double c1,c2,cc1;

  // get bin number in Z
  G4int iZ = 14;
  // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  while ((iZ>=0)&&(Zdat[iZ]>=atomicNumber)) { --iZ; }

  iZ = std::min(std::max(iZ, 0), 13);

  G4double ZZ1 = Zdat[iZ];
  G4double ZZ2 = Zdat[iZ+1];
  G4double ratZ = (atomicNumber-ZZ1)*(atomicNumber+ZZ1)/
                  ((ZZ2-ZZ1)*(ZZ2+ZZ1));

  static const G4double Tlim = 10.*CLHEP::MeV;
  static const G4double sigmafactor =
    CLHEP::twopi*CLHEP::classic_electr_radius*CLHEP::classic_electr_radius;
  static const G4double beta2lim = Tlim*(Tlim+2.*CLHEP::electron_mass_c2)/
    ((Tlim+CLHEP::electron_mass_c2)*(Tlim+CLHEP::electron_mass_c2));
  static const G4double bg2lim   = Tlim*(Tlim+2.*CLHEP::electron_mass_c2)/
    (CLHEP::electron_mass_c2*CLHEP::electron_mass_c2);

  static const G4double sig0[15] = {
    0.2672*CLHEP::barn,  0.5922*CLHEP::barn,  2.653*CLHEP::barn, 6.235*CLHEP::barn,
    11.69*CLHEP::barn  , 13.24*CLHEP::barn  , 16.12*CLHEP::barn, 23.00*CLHEP::barn,
    35.13*CLHEP::barn  , 39.95*CLHEP::barn  , 50.85*CLHEP::barn, 67.19*CLHEP::barn,
    91.15*CLHEP::barn  , 104.4*CLHEP::barn  , 113.1*CLHEP::barn};

  static const G4double Tdat[22] = { 
    100*CLHEP::eV,  200*CLHEP::eV,  400*CLHEP::eV,  700*CLHEP::eV,
    1*CLHEP::keV,   2*CLHEP::keV,   4*CLHEP::keV,   7*CLHEP::keV,
    10*CLHEP::keV,  20*CLHEP::keV,  40*CLHEP::keV,  70*CLHEP::keV,
    100*CLHEP::keV, 200*CLHEP::keV, 400*CLHEP::keV, 700*CLHEP::keV,
    1*CLHEP::MeV,   2*CLHEP::MeV,   4*CLHEP::MeV,   7*CLHEP::MeV,
    10*CLHEP::MeV,  20*CLHEP::MeV};

  if(eKineticEnergy <= Tlim) 
  {
    // get bin number in T (beta2)
    G4int iT = 21;
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    while ((iT>=0)&&(Tdat[iT]>=eKineticEnergy)) iT -= 1;

    iT = std::min(std::max(iT, 0), 20);

    //  calculate betasquare values
    G4double T = Tdat[iT];
    G4double E = T + CLHEP::electron_mass_c2;
    G4double b2small = T*(E+CLHEP::electron_mass_c2)/(E*E);

    T = Tdat[iT+1];
    E = T + CLHEP::electron_mass_c2;
    G4double b2big = T*(E+CLHEP::electron_mass_c2)/(E*E);
    G4double ratb2 = (beta2-b2small)/(b2big-b2small);

    if (charge < 0.)
    {
       c1 = celectron[iZ][iT];
       c2 = celectron[iZ+1][iT];
       cc1 = c1+ratZ*(c2-c1);

       c1 = celectron[iZ][iT+1];
       c2 = celectron[iZ+1][iT+1];
    }
    else              
    {
       c1 = cpositron[iZ][iT];
       c2 = cpositron[iZ+1][iT];
       cc1 = c1+ratZ*(c2-c1);

       c1 = cpositron[iZ][iT+1];
       c2 = cpositron[iZ+1][iT+1];
    }
    G4double cc2 = c1+ratZ*(c2-c1);
    sigma *= sigmafactor/(cc1+ratb2*(cc2-cc1));
  }
  else
  {
    c1 = bg2lim*sig0[iZ]*(1.+hecorr[iZ]*(beta2-beta2lim))/bg2;
    c2 = bg2lim*sig0[iZ+1]*(1.+hecorr[iZ+1]*(beta2-beta2lim))/bg2;
    if((atomicNumber >= ZZ1) && (atomicNumber <= ZZ2))
      sigma = c1+ratZ*(c2-c1) ;
    else if(atomicNumber < ZZ1)
      sigma = atomicNumber*atomicNumber*c1/(ZZ1*ZZ1);
    else if(atomicNumber > ZZ2)
      sigma = atomicNumber*atomicNumber*c2/(ZZ2*ZZ2);
  }
  // low energy correction based on theory 
  sigma *= (1.+0.30/(1.+std::sqrt(1000.*eKineticEnergy)));  

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UrbanMscModel::StartTracking(G4Track* track)
{
  SetParticle(track->GetDynamicParticle()->GetDefinition());
  firstStep = true; 
  insideskin = false;
  fr = facrange;
  tlimit = tgeom = rangeinit = geombig;
  smallstep     = 1.e10;
  stepmin       = tlimitminfix;
  tlimitmin     = 10.*tlimitminfix;            
  rndmEngineMod = G4Random::getTheEngine();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UrbanMscModel::ComputeTruePathLengthLimit(
                             const G4Track& track,
                             G4double& currentMinimalStep)
{
  tPathLength = currentMinimalStep;
  const G4DynamicParticle* dp = track.GetDynamicParticle();
  
  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  G4StepStatus stepStatus = sp->GetStepStatus();
  couple = track.GetMaterialCutsCouple();
  SetCurrentCouple(couple); 
  idx = couple->GetIndex();
  currentKinEnergy = dp->GetKineticEnergy();
  currentLogKinEnergy = dp->GetLogKineticEnergy();
  currentRange = GetRange(particle,currentKinEnergy,couple,currentLogKinEnergy);
  lambda0 = GetTransportMeanFreePath(particle,currentKinEnergy,
                                              currentLogKinEnergy);
  tPathLength = std::min(tPathLength,currentRange);
  /*  
  G4cout << "G4Urban::StepLimit tPathLength= " << tPathLength 
  << " range= " <<currentRange<< " lambda= "<<lambda0
            <<G4endl;
  */
  
  // stop here if small step
  if(tPathLength < tlimitminfix) { 
    latDisplasment = false;   
    return ConvertTrueToGeom(tPathLength, currentMinimalStep); 
  }

  // upper limit for the straight line distance the particle can travel
  // for electrons and positrons
  G4double distance = (mass < masslimite) 
    ? currentRange*msc[idx]->doverra
    // for muons, hadrons
    : currentRange*msc[idx]->doverrb;

  presafety = (stepStatus == fGeomBoundary) ? sp->GetSafety()
              : ComputeSafety(sp->GetPosition(),tPathLength);
  /*  
  G4cout << "G4Urban::StepLimit tPathLength= " 
            <<tPathLength<<" safety= " << presafety
          << " range= " <<currentRange<< " lambda= "<<lambda0
            << " Alg: " << steppingAlgorithm <<G4endl;
  */
  // far from geometry boundary
  if(distance < presafety)
    {
      latDisplasment = false;   
      return ConvertTrueToGeom(tPathLength, currentMinimalStep);  
    }

  latDisplasment = latDisplasmentbackup;
  // ----------------------------------------------------------------
  // distance to boundary 
  if (steppingAlgorithm == fUseDistanceToBoundary)
    {
      //compute geomlimit and presafety 
      geomlimit = ComputeGeomLimit(track, presafety, currentRange);
      /*
        G4cout << "G4Urban::Distance to boundary geomlimit= "
            <<geomlimit<<" safety= " << presafety<<G4endl;
      */

      smallstep += 1.;
      insideskin = false;

      // initialisation at firs step and at the boundary
      if(firstStep || (stepStatus == fGeomBoundary))
        {
          rangeinit = currentRange;
          if(!firstStep) { smallstep = 1.; }

          //stepmin ~ lambda_elastic
          stepmin = ComputeStepmin();
          skindepth = skin*stepmin;
          tlimitmin = ComputeTlimitmin();
        /* 
          G4cout << "rangeinit= " << rangeinit << " stepmin= " << stepmin
                 << " tlimitmin= " << tlimitmin << " geomlimit= " 
                 << geomlimit <<G4endl;
        */
          // constraint from the geometry

          if((geomlimit < geombig) && (geomlimit > geommin))
            {
              // geomlimit is a geometrical step length
              // transform it to true path length (estimation)
              if(lambda0 > geomlimit) {
                geomlimit = -lambda0*G4Log(1.-geomlimit/lambda0)+tlimitmin;
	      }
              tgeom = (stepStatus == fGeomBoundary)
                ? geomlimit/facgeom : 2.*geomlimit/facgeom;
            }
          else
	    {
	      tgeom = geombig;
	    }
        }

      //step limit 
      tlimit = (currentRange > presafety) ?
        std::max(facrange*rangeinit, facsafety*presafety) : currentRange;

      //lower limit for tlimit
      tlimit = std::min(std::max(tlimit,tlimitmin), tgeom);
      /*
      G4cout << "tgeom= " << tgeom << " geomlimit= " << geomlimit  
            << " tlimit= " << tlimit << " presafety= " << presafety << G4endl;
      */
      // shortcut
      if((tPathLength < tlimit) && (tPathLength < presafety) &&
         (smallstep > skin) && (tPathLength < geomlimit-0.999*skindepth))
      {
        return ConvertTrueToGeom(tPathLength, currentMinimalStep);   
      }

      // step reduction near to boundary
      if(smallstep <= skin)
        {
          tlimit = stepmin;
          insideskin = true;
        }
      else if(geomlimit < geombig)
        {
          if(geomlimit > skindepth)
            {
              tlimit = std::min(tlimit, geomlimit-0.999*skindepth);
            }
          else
            {
              insideskin = true;
              tlimit = std::min(tlimit, stepmin);
            }
        }

      tlimit = std::max(tlimit, stepmin); 

      // randomise if not 'small' step and step determined by msc
      tPathLength = (tlimit < tPathLength && smallstep > skin && !insideskin) 
        ? std::min(tPathLength, Randomizetlimit()) 
	: std::min(tPathLength, tlimit);
    }
  // ----------------------------------------------------------------
  // for simulation with or without magnetic field 
  // there no small step/single scattering at boundaries
  else if(steppingAlgorithm == fUseSafety)
    {
      if(firstStep || (stepStatus == fGeomBoundary)) {
        rangeinit = currentRange;
        fr = facrange;
        // stepping for e+/e- only (not for muons,hadrons)
        if(mass < masslimite) 
          {
            rangeinit = std::max(rangeinit, lambda0);
            if(lambda0 > lambdalimit) {
              fr *= (0.75+0.25*lambda0/lambdalimit);
            }
          }
        //lower limit for tlimit
        stepmin = ComputeStepmin();
        tlimitmin = ComputeTlimitmin();
      }

      //step limit
      tlimit = (currentRange > presafety) ?
        std::max(fr*rangeinit, facsafety*presafety) : currentRange;
  
      //lower limit for tlimit
      tlimit = std::max(tlimit, tlimitmin); 
     
      // randomise if step determined by msc
      tPathLength = (tlimit < tPathLength) ?
        std::min(tPathLength, Randomizetlimit()) : tPathLength; 
    }
  // ----------------------------------------------------------------
  // for simulation with or without magnetic field 
  // there is small step/single scattering at boundaries
  else if(steppingAlgorithm == fUseSafetyPlus)
    {
      if(firstStep || (stepStatus == fGeomBoundary)) {
        rangeinit = currentRange;
        fr = facrange;
        if(mass < masslimite)
          {
            if(lambda0 > lambdalimit) {
              fr *= (0.84+0.16*lambda0/lambdalimit);
            }
          }
        //lower limit for tlimit
        stepmin = ComputeStepmin();
        tlimitmin = ComputeTlimitmin();
      }
      //step limit
      tlimit = (currentRange > presafety) ?
	std::max(fr*rangeinit, facsafety*presafety) : currentRange;

      //lower limit for tlimit
      tlimit = std::max(tlimit, tlimitmin);

      // condition for tPathLength from drr and finalr
      if(currentRange > finalr) {
        G4double tmax = drr*currentRange+
                        finalr*(1.-drr)*(2.-finalr/currentRange);
        tPathLength = std::min(tPathLength,tmax); 
      }

      // randomise if step determined by msc
      tPathLength = (tlimit < tPathLength) ?
        std::min(tPathLength, Randomizetlimit()) : tPathLength; 
    }

  // ----------------------------------------------------------------
  // simple step limitation
  else
    {
      if (stepStatus == fGeomBoundary)
        {
          tlimit = (currentRange > lambda0) 
	    ? facrange*currentRange : facrange*lambda0;
          tlimit = std::max(tlimit, tlimitmin);
        }
      // randomise if step determined by msc
      tPathLength = (tlimit < tPathLength) ?
        std::min(tPathLength, Randomizetlimit()) : tPathLength; 
    }

  // ----------------------------------------------------------------
  firstStep = false; 
  return ConvertTrueToGeom(tPathLength, currentMinimalStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UrbanMscModel::ComputeGeomPathLength(G4double)
{
  lambdaeff = lambda0;
  par1 = -1. ;  
  par2 = par3 = 0. ;  

  // this correction needed to run MSC with eIoni and eBrem inactivated
  // and makes no harm for a normal run
  tPathLength = std::min(tPathLength,currentRange); 

  //  do the true -> geom transformation
  zPathLength = tPathLength;

  // z = t for very small tPathLength
  if(tPathLength < tlimitminfix2) return zPathLength;

  /*
  G4cout << "ComputeGeomPathLength: tpl= " <<  tPathLength
         << " R= " << currentRange << " L0= " << lambda0
         << " E= " << currentKinEnergy << "  " 
         << particle->GetParticleName() << G4endl;
  */
  G4double tau = tPathLength/lambda0 ;

  if ((tau <= tausmall) || insideskin) {
    zPathLength = std::min(tPathLength, lambda0); 

  } else if (tPathLength < currentRange*dtrl) {
    zPathLength = (tau < taulim) ? tPathLength*(1.-0.5*tau)
      : lambda0*(1.-G4Exp(-tau));

  } else if(currentKinEnergy < mass || tPathLength == currentRange)  {
    par1 = 1./currentRange;
    par2 = currentRange/lambda0;
    par3 = 1.+par2;
    if(tPathLength < currentRange) {
      zPathLength = 
        (1.-G4Exp(par3*G4Log(1.-tPathLength/currentRange)))/(par1*par3);
    } else {
      zPathLength = 1./(par1*par3);
    }

  } else {
    G4double rfin = std::max(currentRange-tPathLength, 0.01*currentRange);
    G4double T1 = GetEnergy(particle,rfin,couple);
    G4double lambda1 = GetTransportMeanFreePath(particle,T1);

    par1 = (lambda0-lambda1)/(lambda0*tPathLength);
    //G4cout << "par1= " << par1 << " L1= " << lambda1 << G4endl;
    par2 = 1./(par1*lambda0);
    par3 = 1.+par2;
    zPathLength = (1.-G4Exp(par3*G4Log(lambda1/lambda0)))/(par1*par3);
  }

  zPathLength = std::min(zPathLength, lambda0);
  //G4cout<< "zPathLength= "<< zPathLength<< " L0= " << lambda0 << G4endl;
  return zPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UrbanMscModel::ComputeTrueStepLength(G4double geomStepLength)
{
  // step defined other than transportation 
  if(geomStepLength == zPathLength) { 
    //G4cout << "Urban::ComputeTrueLength: tPathLength= " << tPathLength 
    //           << " step= " << geomStepLength << " *** " << G4endl;
    return tPathLength; 
  }

  zPathLength = geomStepLength;

  // t = z for very small step
  if(geomStepLength < tlimitminfix2) { 
    tPathLength = geomStepLength; 
  
  // recalculation
  } else {

    G4double tlength = geomStepLength;
    if((geomStepLength > lambda0*tausmall) && !insideskin) {

      if(par1 <  0.) {
        tlength = -lambda0*G4Log(1.-geomStepLength/lambda0) ;
      } else {
        const G4double par4 = par1*par3;
        if(par4*geomStepLength < 1.) {
          tlength = (1.-G4Exp(G4Log(1.-par4*geomStepLength)/par3))/par1;
        } else {
          tlength = currentRange;
        }
      }

      if(tlength < geomStepLength)   { tlength = geomStepLength; }
      else if(tlength > tPathLength) { tlength = tPathLength; }
    }  
    tPathLength = tlength; 
  }
  //G4cout << "Urban::ComputeTrueLength: tPathLength= " << tPathLength 
  //         << " step= " << geomStepLength << " &&& " << G4endl;

  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector& 
G4UrbanMscModel::SampleScattering(const G4ThreeVector& oldDirection,
                                  G4double /*safety*/)
{
  fDisplacement.set(0.0,0.0,0.0);
  G4double kinEnergy = currentKinEnergy;
  if (tPathLength > currentRange*dtrl) {
    kinEnergy = GetEnergy(particle,currentRange-tPathLength,couple);
  } else if(tPathLength > currentRange*0.01) {
    kinEnergy -= tPathLength*GetDEDX(particle,currentKinEnergy,couple,
                                     currentLogKinEnergy);
  }

  if((tPathLength <= tlimitminfix) || (tPathLength < tausmall*lambda0) || 
     (kinEnergy <= CLHEP::eV)) { return fDisplacement; }

  G4double cth = SampleCosineTheta(tPathLength,kinEnergy);

  // protection against 'bad' cth values
  if(std::abs(cth) >= 1.0) { return fDisplacement; } 

  G4double sth = std::sqrt((1.0 - cth)*(1.0 + cth));
  G4double phi = CLHEP::twopi*rndmEngineMod->flat();
  G4ThreeVector newDirection(sth*std::cos(phi),sth*std::sin(phi),cth);
  newDirection.rotateUz(oldDirection);

  fParticleChange->ProposeMomentumDirection(newDirection);
  /*
  G4cout << "G4UrbanMscModel::SampleSecondaries: e(MeV)= " << kineticEnergy
         << " sinTheta= " << sth << " safety(mm)= " << safety
         << " trueStep(mm)= " << tPathLength
         << " geomStep(mm)= " << zPathLength
         << G4endl;
  */

  if (latDisplasment && currentTau >= tausmall) {
    if(dispAlg96) { SampleDisplacement(sth, phi); }
    else          { SampleDisplacementNew(cth, phi); }
    fDisplacement.rotateUz(oldDirection);
  }
  return fDisplacement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UrbanMscModel::SampleCosineTheta(G4double trueStepLength,
                                            G4double kinEnergy)
{
  G4double cth = 1.0;
  G4double tau = trueStepLength/lambda0;

  // mean tau value
  if(currentKinEnergy != kinEnergy) {
    G4double lambda1 = GetTransportMeanFreePath(particle, kinEnergy);
    if(std::abs(lambda1 - lambda0) > lambda0*0.01 && lambda1 > 0.) {
      tau = trueStepLength*G4Log(lambda0/lambda1)/(lambda0-lambda1);
    }
  }

  currentTau = tau;
  lambdaeff = trueStepLength/currentTau;
  currentRadLength = couple->GetMaterial()->GetRadlen();

  if (tau >= taubig) { cth = -1.+2.*rndmEngineMod->flat(); }
  else if (tau >= tausmall) {
    static const G4double numlim = 0.01;
    static const G4double onethird = 1./3.;
    if(tau < numlim) {
      xmeanth = 1.0 - tau*(1.0 - 0.5*tau);
      x2meanth= 1.0 - tau*(5.0 - 6.25*tau)*onethird;
    } else {
      xmeanth = G4Exp(-tau);
      x2meanth = (1.+2.*G4Exp(-2.5*tau))*onethird;
    }

    // too large step of low-energy particle
    G4double relloss = 1. - kinEnergy/currentKinEnergy;
    static const G4double rellossmax= 0.50;
    if(relloss > rellossmax) {
      return SimpleScattering();
    }
    // is step extreme small ?
    G4bool extremesmallstep = false;
    G4double tsmall = std::min(tlimitmin,lambdalimit);

    G4double theta0;
    if(trueStepLength > tsmall) {
      theta0 = ComputeTheta0(trueStepLength,kinEnergy);
    } else {
      theta0 = std::sqrt(trueStepLength/tsmall)
	*ComputeTheta0(tsmall,kinEnergy);
      extremesmallstep = true;
    }

    static const G4double onesixth = 1./6.;
    static const G4double one12th = 1./12.;
    static const G4double theta0max = CLHEP::pi*onesixth;

    // protection for very small angles
    G4double theta2 = theta0*theta0;

    if(theta2 < tausmall) { return cth; }
    if(theta0 > theta0max) { return SimpleScattering(); }

    G4double x = theta2*(1.0 - theta2*one12th);
    if(theta2 > numlim) {
      G4double sth = 2*std::sin(0.5*theta0);
      x = sth*sth;
    }

    // parameter for tail
    G4double ltau = G4Log(tau);
    G4double u = !extremesmallstep ? G4Exp(ltau*onesixth) 
      : G4Exp(G4Log(tsmall/lambda0)*onesixth); 

    G4double xx  = G4Log(lambdaeff/currentRadLength);
    G4double xsi = msc[idx]->coeffc1 + 
      u*(msc[idx]->coeffc2+msc[idx]->coeffc3*u)+msc[idx]->coeffc4*xx;

    // tail should not be too big
    xsi = std::max(xsi, 1.9); 
      /*
      if(KineticEnergy > 20*MeV && xsi < 1.6) {
        G4cout << "G4UrbanMscModel::SampleCosineTheta: E(GeV)= " 
               << KineticEnergy/GeV 
               << " !!** c= " << xsi
               << " **!! length(mm)= " << trueStepLength << " Zeff= " << Zeff 
               << " " << couple->GetMaterial()->GetName()
               << " tau= " << tau << G4endl;
      }
      */

    G4double c = xsi;

    if(std::abs(c-3.) < 0.001)      { c = 3.001; }
    else if(std::abs(c-2.) < 0.001) { c = 2.001; }

    G4double c1 = c-1.;
    G4double ea = G4Exp(-xsi);
    G4double eaa = 1.-ea ;
    G4double xmean1 = 1.-(1.-(1.+xsi)*ea)*x/eaa;
    G4double x0 = 1. - xsi*x;

    // G4cout << " xmean1= " << xmean1 << "  xmeanth= " << xmeanth << G4endl;

    if(xmean1 <= 0.999*xmeanth) { return SimpleScattering(); }

    //from continuity of derivatives
    G4double b = 1.+(c-xsi)*x;

    G4double b1 = b+1.;
    G4double bx = c*x;

    G4double eb1 = G4Exp(G4Log(b1)*c1);
    G4double ebx = G4Exp(G4Log(bx)*c1);
    G4double d = ebx/eb1;

    G4double xmean2 = (x0 + d - (bx - b1*d)/(c-2.))/(1. - d);

    G4double f1x0 = ea/eaa;
    G4double f2x0 = c1/(c*(1. - d));
    G4double prob = f2x0/(f1x0+f2x0);

    G4double qprob = xmeanth/(prob*xmean1+(1.-prob)*xmean2);

    // sampling of costheta
    //G4cout << "c= " << c << " qprob= " << qprob << " eb1= " << eb1
    // << " c1= " << c1 << " b1= " << b1 << " bx= " << bx << " eb1= " << eb1
    //             << G4endl;
    rndmEngineMod->flatArray(2, rndmarray);
    if(rndmarray[0] < qprob)
    {
      G4double var = 0;
      if(rndmarray[1] < prob) {
        cth = 1.+G4Log(ea+rndmEngineMod->flat()*eaa)*x;
      } else {
        var = (1.0 - d)*rndmEngineMod->flat();
        if(var < numlim*d) {
          var /= (d*c1); 
          cth = -1.0 + var*(1.0 - 0.5*var*c)*(2. + (c - xsi)*x);
        } else {
          cth = 1. + x*(c - xsi - c*G4Exp(-G4Log(var + d)/c1));
        }
      } 
    } else {
      cth = -1.+2.*rndmarray[1];
    }
  }
  return cth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UrbanMscModel::ComputeTheta0(G4double trueStepLength,
                                        G4double kinEnergy)
{
  // for all particles take the width of the central part
  //  from a  parametrization similar to the Highland formula
  // ( Highland formula: Particle Physics Booklet, July 2002, eq. 26.10)
  G4double invbetacp = (kinEnergy+mass)/(kinEnergy*(kinEnergy+2.*mass));
  if(currentKinEnergy != kinEnergy) {
    invbetacp = std::sqrt(invbetacp*(currentKinEnergy+mass)/
			  (currentKinEnergy*(currentKinEnergy+2.*mass)));
  }
  G4double y = trueStepLength/currentRadLength;

  if(fPosiCorrection && particle == positron)
  {
    static const G4double xl= 0.6;
    static const G4double xh= 0.9;
    static const G4double e = 113.0;
    G4double corr;

    G4double tau = std::sqrt(currentKinEnergy*kinEnergy)/mass;
    G4double x = std::sqrt(tau*(tau+2.)/((tau+1.)*(tau+1.)));
    G4double a = msc[idx]->posa; 
    G4double b = msc[idx]->posb;
    G4double c = msc[idx]->posc; 
    G4double d = msc[idx]->posd;
    if(x < xl) {
      corr = a*(1.-G4Exp(-b*x));  
    } else if(x > xh) {
      corr = c+d*G4Exp(e*(x-1.)); 
    } else {
      G4double yl = a*(1.-G4Exp(-b*xl));
      G4double yh = c+d*G4Exp(e*(xh-1.));
      G4double y0 = (yh-yl)/(xh-xl);
      G4double y1 = yl-y0*xl;
      corr = y0*x+y1;
    }
    //==================================================================
    y *= corr*msc[idx]->pose;
  }

  static const G4double c_highland = 13.6*CLHEP::MeV;
  G4double theta0 = c_highland*std::abs(charge)*std::sqrt(y)*invbetacp;
 
  // correction factor from e- scattering data
  theta0 *= (msc[idx]->coeffth1+msc[idx]->coeffth2*G4Log(y));
  return theta0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UrbanMscModel::SampleDisplacement(G4double, G4double phi)
{ 
  // simple and fast sampling
  // based on single scattering results
  // u = r/rmax : mean value

  G4double rmax = std::sqrt((tPathLength-zPathLength)*(tPathLength+zPathLength));
  if(rmax > 0.)
  {
    G4double r = 0.73*rmax;

    // simple distribution for v=Phi-phi=psi ~exp(-beta*v)
    // beta determined from the requirement that distribution should give
    // the same mean value than that obtained from the ss simulation

    static const G4double cbeta  = 2.160;
    static const G4double cbeta1 = 1. - G4Exp(-cbeta*CLHEP::pi);
    rndmEngineMod->flatArray(2, rndmarray);
    G4double psi = -G4Log(1. - rndmarray[0]*cbeta1)/cbeta;
    G4double Phi = (rndmarray[1] < 0.5) ? phi+psi : phi-psi;
    fDisplacement.set(r*std::cos(Phi),r*std::sin(Phi),0.0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
   
void G4UrbanMscModel::SampleDisplacementNew(G4double, G4double phi)
{
  // best sampling based on single scattering results
  G4double rmax = 
    std::sqrt((tPathLength-zPathLength)*(tPathLength+zPathLength));
  G4double r(0.0);
  G4double u(0.0);
  static const G4double reps = 5.e-3;

  if(rmax > 0.)
  {
    static const G4double umax = 0.855;
    static const G4double wlow = 0.750;

    static const G4double  ralpha = 6.83e+0;          
    static const G4double  ra1 =-4.16179e+1;       
    static const G4double  ra2 = 1.12548e+2;       
    static const G4double  ra3 =-8.66665e+1;        
    static const G4double  ralpha1 = 0.751*ralpha;        
    static const G4double  ralpha2 =ralpha-ralpha1;        
    static const G4double  rwa1 = G4Exp(ralpha1*reps);
    static const G4double  rwa2 = G4Exp(ralpha1*umax)-rwa1;
    static const G4double  rejamax = 1.16456;

    static const G4double  rbeta = 2.18e+1; 
    static const G4double  rb0 = 4.81382e+2;
    static const G4double  rb1 =-1.12842e+4;
    static const G4double  rb2 = 4.57745e+4;
    static const G4double  rbeta1 = 0.732*rbeta;
    static const G4double  rbeta2 = rbeta-rbeta1;
    static const G4double  rwb1 = G4Exp(-rbeta1*umax);
    static const G4double  rwb2 = rwb1-G4Exp(-rbeta1*(1.-reps));
    static const G4double  rejbmax = 1.62651;

    G4int count = 0;
    G4double uc,rej;

    if(rndmEngineMod->flat() < wlow)
    {
      do {
	   rndmEngineMod->flatArray(2, rndmarray);
           u = G4Log(rwa1+rwa2*rndmarray[0])/ralpha1;
           uc = umax-u;
           rej = G4Exp(-ralpha2*uc)*
                (1.+ralpha*uc+ra1*uc*uc+ra2*uc*uc*uc+ra3*uc*uc*uc*uc);
         } while (rejamax*rndmarray[1] > rej && ++count < 1000);
    }
    else
    {
      do {
	   rndmEngineMod->flatArray(2, rndmarray);
           u = -G4Log(rwb1-rwb2*rndmarray[0])/rbeta1;
           uc = u-umax;
           rej = G4Exp(-rbeta2*uc)*
                (1.+rbeta*uc+rb0*uc*uc+rb1*uc*uc*uc+rb2*uc*uc*uc*uc);
         } while (rejbmax*rndmarray[1] > rej && ++count < 1000);
    }
    r = rmax*u;
  }
                               
  if(r > 0.)
  {
    // sample Phi using lateral correlation
    // and r/rmax - (Phi-phi) correlation
    // v = Phi-phi = acos(latcorr/(r*sth))
    // from SS simulation f(v)*g(v)               
    // f(v) ~ exp(-a1*v) normalized distribution
    // g(v) rejection function (0 < g(v) <= 1)
    G4double v, rej;
                          
    static const G4double peps = 1.e-4;
    static const G4double palpha[10] = {2.300e+0,2.490e+0,2.610e+0,2.820e+0,2.710e+0,
                                        2.750e+0,2.910e+0,3.400e+0,4.150e+0,5.400e+0};
    static const G4double palpha1[10]= {4.600e-2,1.245e-1,2.610e-1,2.820e-1,2.710e-1,
                                        6.875e-1,1.019e+0,1.360e+0,1.660e+0,2.430e+0};
    static const G4double pejmax[10] = {3.513,1.968,1.479,1.239,1.116,
                                        1.081,1.064,1.073,1.103,1.158};

    static const G4double pa1[10] = { 3.218e+0, 2.412e+0, 2.715e+0, 2.787e+0, 2.541e+0,
                                      2.508e+0, 2.600e+0, 3.231e+0, 4.588e+0, 6.584e+0};
    static const G4double pa2[10] = {-5.528e-1, 2.523e+0, 1.738e+0, 2.082e+0, 1.423e+0,
                                      4.682e-1,-6.883e-1,-2.147e+0,-5.127e+0,-1.054e+1};
    static const G4double pa3[10] = { 3.618e+0, 2.032e+0, 2.341e+0, 2.172e+0, 7.205e-1,
                                      4.655e-1, 6.318e-1, 1.255e+0, 2.425e+0, 4.938e+0};
    static const G4double pa4[10] = { 2.437e+0, 9.450e-1, 4.349e-1, 2.221e-1, 1.130e-1,
                                      5.405e-2, 2.245e-2, 7.370e-3, 1.456e-3, 1.508e-4};
    static const G4double pw1[10] = {G4Exp(-palpha1[0]*peps),G4Exp(-palpha1[1]*peps),
                                     G4Exp(-palpha1[2]*peps),G4Exp(-palpha1[3]*peps),
                                     G4Exp(-palpha1[4]*peps),G4Exp(-palpha1[5]*peps),
                                     G4Exp(-palpha1[6]*peps),G4Exp(-palpha1[7]*peps),
                                     G4Exp(-palpha1[8]*peps),G4Exp(-palpha1[9]*peps)};
    static const G4double pw2[10] = {pw1[0]-G4Exp(-palpha1[0]*(CLHEP::pi-peps)),
                                     pw1[1]-G4Exp(-palpha1[1]*(CLHEP::pi-peps)),
                                     pw1[2]-G4Exp(-palpha1[2]*(CLHEP::pi-peps)),
                                     pw1[3]-G4Exp(-palpha1[3]*(CLHEP::pi-peps)),
                                     pw1[4]-G4Exp(-palpha1[4]*(CLHEP::pi-peps)),
                                     pw1[5]-G4Exp(-palpha1[5]*(CLHEP::pi-peps)),
                                     pw1[6]-G4Exp(-palpha1[6]*(CLHEP::pi-peps)),
                                     pw1[7]-G4Exp(-palpha1[7]*(CLHEP::pi-peps)),
                                     pw1[8]-G4Exp(-palpha1[8]*(CLHEP::pi-peps)),
                                     pw1[9]-G4Exp(-palpha1[9]*(CLHEP::pi-peps))};

    G4int iphi = (G4int)(u*10.);
    if(iphi < 0)      { iphi = 0; }
    else if(iphi > 9) { iphi = 9; }
    G4int count = 0;

    do {
      rndmEngineMod->flatArray(2, rndmarray);
      v = -G4Log(pw1[iphi]-pw2[iphi]*rndmarray[0])/palpha1[iphi];
      rej = (G4Exp(-palpha[iphi]*v)*
            (1+pa1[iphi]*v+pa2[iphi]*v*v+pa3[iphi]*v*v*v)+pa4[iphi])/
            G4Exp(-pw1[iphi]*v);
    }
    // Loop checking, 5-March-2018, Vladimir Ivanchenko
    while (pejmax[iphi]*rndmarray[1] > rej && ++count < 1000);

    G4double Phi = (rndmEngineMod->flat() < 0.5) ? phi+v : phi-v;
    fDisplacement.set(r*std::cos(Phi),r*std::sin(Phi),0.0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UrbanMscModel::InitialiseModelCache()                                   
{
  // it is assumed, that for the second run only addition
  // of a new G4MaterialCutsCouple is possible
  auto theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();
  if(numOfCouples != msc.size()) { msc.resize(numOfCouples, nullptr); }
  
  for(G4int j=0; j<(G4int)numOfCouples; ++j) {
    auto aCouple = theCoupleTable->GetMaterialCutsCouple(j);

    // new couple
    msc[j] = new mscData();
    G4double Zeff = aCouple->GetMaterial()->GetIonisation()->GetZeffective();
    msc[j]->sqrtZ = std::sqrt(Zeff);
    G4double lnZ = G4Log(Zeff);
    // correction in theta0 formula
    G4double w = G4Exp(lnZ/6.);
    G4double facz = 0.990395+w*(-0.168386+w*0.093286);
    msc[j]->coeffth1 = facz*(1. - 8.7780e-2/Zeff);
    msc[j]->coeffth2 = facz*(4.0780e-2 + 1.7315e-4*Zeff);

    // tail parameters
    G4double Z13 = w*w;
    msc[j]->coeffc1 = 2.3785    - Z13*(4.1981e-1 - Z13*6.3100e-2);
    msc[j]->coeffc2 = 4.7526e-1 + Z13*(1.7694    - Z13*3.3885e-1);
    msc[j]->coeffc3 = 2.3683e-1 - Z13*(1.8111    - Z13*3.2774e-1);
    msc[j]->coeffc4 = 1.7888e-2 + Z13*(1.9659e-2 - Z13*2.6664e-3);

    msc[j]->Z23 = Z13*Z13;       

    msc[j]->stepmina = 27.725/(1.+0.203*Zeff);
    msc[j]->stepminb =  6.152/(1.+0.111*Zeff);

    // 21.07.2020
    msc[j]->doverra = 9.6280e-1 - 8.4848e-2*msc[j]->sqrtZ + 4.3769e-3*Zeff;

    // 06.10.2020 
    // msc[j]->doverra = 7.7024e-1 - 6.7878e-2*msc[j]->sqrtZ + 3.5015e-3*Zeff;
    msc[j]->doverrb = 1.15 - 9.76e-4*Zeff;

    // corrections for e+
    msc[j]->posa = 0.994-4.08e-3*Zeff;
    msc[j]->posb = 7.16+(52.6+365./Zeff)/Zeff;
    msc[j]->posc = 1.000-4.47e-3*Zeff;
    msc[j]->posd = 1.21e-3*Zeff;
    msc[j]->pose = 1.+Zeff*(1.84035e-4*Zeff-1.86427e-2)+0.41125; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
