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
// $Id:   $
// GEANT4 tag $Name:  $
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
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UrbanMscModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"
#include "G4ParticleChangeForMSC.hh"

#include "G4Poisson.hh"
#include "G4Pow.hh"
#include "globals.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UrbanMscModel::G4UrbanMscModel(const G4String& nam)
  : G4VMscModel(nam)
{
  masslimite    = 0.6*MeV;
  lambdalimit   = 1.*mm;
  fr            = 0.02;
  taubig        = 8.0;
  tausmall      = 1.e-16;
  taulim        = 1.e-6;
  currentTau    = taulim;
  tlimitminfix  = 0.01*nm;             
  tlimitminfix2 =   1.*nm;             
  stepmin       = tlimitminfix;
  smallstep     = 1.e10;
  currentRange  = 0. ;
  rangeinit     = 0.;
  tlimit        = 1.e10*mm;
  tlimitmin     = 10.*tlimitminfix;            
  tgeom         = 1.e50*mm;
  geombig       = 1.e50*mm;
  geommin       = 1.e-3*mm;
  geomlimit     = geombig;
  presafety     = 0.*mm;

  facsafety     = 0.6;

  Zold          = 0.;
  Zeff          = 1.;
  Z2            = 1.;                
  Z23           = 1.;                    
  lnZ           = 0.;
  coeffth1      = 0.;
  coeffth2      = 0.;
  coeffc1       = 0.;
  coeffc2       = 0.;
  coeffc3       = 0.;
  coeffc4       = 0.;
  particle      = 0;

  positron      = G4Positron::Positron();
  theManager    = G4LossTableManager::Instance(); 
  rndmEngineMod = G4Random::getTheEngine();

  firstStep     = true; 
  insideskin    = false;
  latDisplasmentbackup = false;
  dispAlg96 = true;

  rangecut = geombig;
  drr      = 0.35 ;
  finalr   = 10.*um ;

  skindepth = skin*stepmin;

  mass = proton_mass_c2;
  charge = ChargeSquare = 1.0;
  currentKinEnergy = currentRadLength = lambda0 = lambdaeff = tPathLength 
    = zPathLength = par1 = par2 = par3 = 0;

  currentMaterialIndex = -1;
  fParticleChange = nullptr;
  couple = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UrbanMscModel::~G4UrbanMscModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UrbanMscModel::Initialise(const G4ParticleDefinition* p,
                                 const G4DataVector&)
{
  // set values of some data members
  SetParticle(p);
  fParticleChange = GetParticleChangeForMSC(p);

  latDisplasmentbackup = latDisplasment;
  dispAlg96 = (G4EmParameters::Instance()->LateralDisplacementAlg96());

  //G4cout << "### G4UrbanMscModel::Initialise done for " 
  //	 << p->GetParticleName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UrbanMscModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* part,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber,G4double,
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

  Z23 = G4Pow::GetInstance()->Z23(G4lrint(AtomicNumber));

  // correction if particle .ne. e-/e+
  // compute equivalent kinetic energy
  // lambda depends on p*beta ....

  G4double eKineticEnergy = KineticEnergy;

  if(mass > electron_mass_c2)
  {
     G4double TAU = KineticEnergy/mass ;
     G4double c = mass*TAU*(TAU+2.)/(electron_mass_c2*(TAU+1.)) ;
     G4double w = c-2. ;
     G4double tau = 0.5*(w+sqrt(w*w+4.*c)) ;
     eKineticEnergy = electron_mass_c2*tau ;
  }

  G4double eTotalEnergy = eKineticEnergy + electron_mass_c2 ;
  G4double beta2 = eKineticEnergy*(eTotalEnergy+electron_mass_c2)
                                 /(eTotalEnergy*eTotalEnergy);
  G4double bg2   = eKineticEnergy*(eTotalEnergy+electron_mass_c2)
                                 /(electron_mass_c2*electron_mass_c2);

  static const G4double epsfactor = 2.*CLHEP::electron_mass_c2*
    CLHEP::electron_mass_c2*CLHEP::Bohr_radius*CLHEP::Bohr_radius
    /(CLHEP::hbarc*CLHEP::hbarc);
  G4double eps = epsfactor*bg2/Z23;

  if     (eps<epsmin)  sigma = 2.*eps*eps;
  else if(eps<epsmax)  sigma = G4Log(1.+2.*eps)-2.*eps/(1.+2.*eps);
  else                 sigma = G4Log(2.*eps)-1.+1./eps;

  sigma *= ChargeSquare*AtomicNumber*AtomicNumber/(beta2*bg2);

  // interpolate in AtomicNumber and beta2 
  G4double c1,c2,cc1,cc2,corr;

  // get bin number in Z
  G4int iZ = 14;
  // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  while ((iZ>=0)&&(Zdat[iZ]>=AtomicNumber)) iZ -= 1;
  if (iZ==14)                               iZ = 13;
  if (iZ==-1)                               iZ = 0 ;

  G4double ZZ1 = Zdat[iZ];
  G4double ZZ2 = Zdat[iZ+1];
  G4double ratZ = (AtomicNumber-ZZ1)*(AtomicNumber+ZZ1)/
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
    if(iT==21)                                  iT = 20;
    if(iT==-1)                                  iT = 0 ;

    //  calculate betasquare values
    G4double T = Tdat[iT],   E = T + electron_mass_c2;
    G4double b2small = T*(E+electron_mass_c2)/(E*E);

    T = Tdat[iT+1]; E = T + electron_mass_c2;
    G4double b2big = T*(E+electron_mass_c2)/(E*E);
    G4double ratb2 = (beta2-b2small)/(b2big-b2small);

    if (charge < 0.)
    {
       c1 = celectron[iZ][iT];
       c2 = celectron[iZ+1][iT];
       cc1 = c1+ratZ*(c2-c1);

       c1 = celectron[iZ][iT+1];
       c2 = celectron[iZ+1][iT+1];
       cc2 = c1+ratZ*(c2-c1);

       corr = cc1+ratb2*(cc2-cc1);

       sigma *= sigmafactor/corr;
    }
    else              
    {
       c1 = cpositron[iZ][iT];
       c2 = cpositron[iZ+1][iT];
       cc1 = c1+ratZ*(c2-c1);

       c1 = cpositron[iZ][iT+1];
       c2 = cpositron[iZ+1][iT+1];
       cc2 = c1+ratZ*(c2-c1);

       corr = cc1+ratb2*(cc2-cc1);

       sigma *= sigmafactor/corr;
    }
  }
  else
  {
    c1 = bg2lim*sig0[iZ]*(1.+hecorr[iZ]*(beta2-beta2lim))/bg2;
    c2 = bg2lim*sig0[iZ+1]*(1.+hecorr[iZ+1]*(beta2-beta2lim))/bg2;
    if((AtomicNumber >= ZZ1) && (AtomicNumber <= ZZ2))
      sigma = c1+ratZ*(c2-c1) ;
    else if(AtomicNumber < ZZ1)
      sigma = AtomicNumber*AtomicNumber*c1/(ZZ1*ZZ1);
    else if(AtomicNumber > ZZ2)
      sigma = AtomicNumber*AtomicNumber*c2/(ZZ2*ZZ2);
  }
  return sigma;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UrbanMscModel::StartTracking(G4Track* track)
{
  SetParticle(track->GetDynamicParticle()->GetDefinition());
  firstStep = true; 
  insideskin = false;
  fr = facrange;
  tlimit = tgeom = rangeinit = rangecut = geombig;
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
  currentMaterialIndex = couple->GetIndex();
  currentKinEnergy = dp->GetKineticEnergy();
  currentRange = GetRange(particle,currentKinEnergy,couple);
  lambda0 = GetTransportMeanFreePath(particle,currentKinEnergy);
  tPathLength = min(tPathLength,currentRange);
  /*
  G4cout << "G4Urban::StepLimit tPathLength= " << tPathLength 
  << " range= " <<currentRange<< " lambda= "<<lambda0
            <<G4endl;
  */
  // set flag to default values
  Zeff = couple->GetMaterial()->GetIonisation()->GetZeffective();
  //         couple->GetMaterial()->GetTotNbOfAtomsPerVolume();

  if(Zold != Zeff)
    UpdateCache();
  
  // stop here if small step
  if(tPathLength < tlimitminfix) { 
    latDisplasment = false;   
    return ConvertTrueToGeom(tPathLength, currentMinimalStep); 
  }

  // upper limit for the straight line distance the particle can travel
  // for electrons and positrons
  G4double distance = currentRange;
  // for muons, hadrons
  if(mass > masslimite) {
    distance *= (1.15-9.76e-4*Zeff);
  } else {
    distance *= (1.20-Zeff*(1.62e-2-9.22e-5*Zeff));
  }
  presafety = sp->GetSafety();
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
  static const G4double invmev = 1.0/CLHEP::MeV;
  // standard  version
  //
  if (steppingAlgorithm == fUseDistanceToBoundary)
    {
      //compute geomlimit and presafety 
      geomlimit = ComputeGeomLimit(track, presafety, currentRange);
      /*
        G4cout << "G4Urban::Distance to boundary geomlimit= "
            <<geomlimit<<" safety= " << presafety<<G4endl;
      */

      // is it far from boundary ?
      if(distance < presafety)
        {
          latDisplasment = false;   
          return ConvertTrueToGeom(tPathLength, currentMinimalStep);   
        }

      smallstep += 1.;
      insideskin = false;

      // initialisation at firs step and at the boundary
      if(firstStep || (stepStatus == fGeomBoundary))
        {
          rangeinit = currentRange;
          if(!firstStep) { smallstep = 1.; }

          //define stepmin here (it depends on lambda!)
          //rough estimation of lambda_elastic/lambda_transport
          G4double rat = currentKinEnergy*invmev;
          rat = 1.e-3/(rat*(10.+rat)) ;
          //stepmin ~ lambda_elastic
          stepmin = rat*lambda0;
          skindepth = skin*stepmin;
          tlimitmin = max(10*stepmin,tlimitminfix);
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
              if((1.-geomlimit/lambda0) > 0.)
                geomlimit = -lambda0*G4Log(1.-geomlimit/lambda0)+tlimitmin ;

              if(stepStatus == fGeomBoundary)
                tgeom = geomlimit/facgeom;
              else
                tgeom = 2.*geomlimit/facgeom;
            }
          else
            tgeom = geombig;
        }

      //step limit 
      tlimit = facrange*rangeinit;              

      //lower limit for tlimit
      tlimit = max(tlimit,tlimitmin);
      tlimit = min(tlimit,tgeom); 
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
              tlimit = min(tlimit, geomlimit-0.999*skindepth);
            }
          else
            {
              insideskin = true;
              tlimit = min(tlimit, stepmin);
            }
        }

      tlimit = max(tlimit, stepmin); 

      // randomise if not 'small' step and step determined by msc
      if((tlimit < tPathLength) && (smallstep > skin) && !insideskin) 
        { 
          tPathLength = min(tPathLength, Randomizetlimit());
        }
      else
        {  
          tPathLength = min(tPathLength, tlimit); 
        }

    }
    // for 'normal' simulation with or without magnetic field 
    //  there no small step/single scattering at boundaries
  else if(steppingAlgorithm == fUseSafety)
    {
      if(stepStatus != fGeomBoundary)  {
        presafety = ComputeSafety(sp->GetPosition(),tPathLength); 
      }
      /*
      G4cout << "presafety= " << presafety
             << " firstStep= " << firstStep
             << " stepStatus= " << stepStatus 
             << G4endl;
      */
      // is far from boundary
      if(distance < presafety)
        {
          latDisplasment = false;
          return ConvertTrueToGeom(tPathLength, currentMinimalStep);  
        }

      if(firstStep || (stepStatus == fGeomBoundary)) {
        rangeinit = currentRange;
        fr = facrange;
        // 9.1 like stepping for e+/e- only (not for muons,hadrons)
        if(mass < masslimite) 
          {
            rangeinit = max(rangeinit, lambda0);
            if(lambda0 > lambdalimit) {
              fr *= (0.75+0.25*lambda0/lambdalimit);
            }
          }
        //lower limit for tlimit
        G4double rat = currentKinEnergy*invmev;
        rat = 1.e-3/(rat*(10 + rat)) ;
        stepmin = lambda0*rat;
        tlimitmin = max(10*stepmin, tlimitminfix);
      }

      //step limit
      tlimit = max(fr*rangeinit, facsafety*presafety);
  
      //lower limit for tlimit
      tlimit = max(tlimit, tlimitmin); 
     
      // randomise if step determined by msc
      if(tlimit < tPathLength)
      {
        tPathLength = min(tPathLength, Randomizetlimit());
      }
      else { tPathLength = min(tPathLength, tlimit); }
    }
  // new stepping mode UseSafetyPlus
  else if(steppingAlgorithm == fUseSafetyPlus)
    {
      if(stepStatus != fGeomBoundary)  {
        presafety = ComputeSafety(sp->GetPosition(),tPathLength);
      }
      /*
      G4cout << "presafety= " << presafety
             << " firstStep= " << firstStep
             << " stepStatus= " << stepStatus
             << G4endl;
      */
      // is far from boundary
      if(distance < presafety)
        {
          latDisplasment = false;
          return ConvertTrueToGeom(tPathLength, currentMinimalStep);
        }

      if(firstStep || (stepStatus == fGeomBoundary)) {
        rangeinit = currentRange;
        fr = facrange;
        rangecut = geombig;
        if(mass < masslimite)
          {
            G4int index = 1;
            if(charge > 0.) index = 2;
            rangecut = couple->GetProductionCuts()->GetProductionCut(index);
            if(lambda0 > lambdalimit) {
              fr *= (0.84+0.16*lambda0/lambdalimit);
            }
          }
        //lower limit for tlimit
        G4double rat = currentKinEnergy*invmev;
        rat = 1.e-3/(rat*(10 + rat)) ;
        stepmin = lambda0*rat;
        tlimitmin = max(10*stepmin, tlimitminfix);
      }
      //step limit
      tlimit = max(fr*rangeinit, facsafety*presafety);

      //lower limit for tlimit
      tlimit = max(tlimit, tlimitmin);

      // condition for tPathLength from drr and finalr
      if(currentRange > finalr) {
        G4double tmax = drr*currentRange+
                        finalr*(1.-drr)*(2.-finalr/currentRange);
        tPathLength = min(tPathLength,tmax); 
      }

      // condition safety
      if(currentRange > rangecut) {
        if(firstStep) {
          tPathLength = min(tPathLength,facsafety*presafety);
        } else if(stepStatus != fGeomBoundary && presafety > stepmin) {
          tPathLength = min(tPathLength,presafety);
        } 
      }

      // randomise if step determined by msc
      if(tPathLength < tlimit)
      {
        tPathLength = min(tPathLength, Randomizetlimit());
      }
      else { tPathLength = min(tPathLength, tlimit); }
    }

  // version similar to 7.1 (needed for some experiments)
  else
    {
      if (stepStatus == fGeomBoundary)
        {
          if (currentRange > lambda0) { tlimit = facrange*currentRange; }
          else                        { tlimit = facrange*lambda0; }

          tlimit = max(tlimit, tlimitmin);
        }
      // randomise if step determined by msc
      if(tlimit < tPathLength)                      
      {
        tPathLength = min(tPathLength, Randomizetlimit());
      }
      else { tPathLength = min(tPathLength, tlimit); }
    }
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

  // VI: it is already checked
  // if(tPathLength > currentRange)
  //  tPathLength = currentRange ;
  /*
  G4cout << "ComputeGeomPathLength: tpl= " <<  tPathLength
         << " R= " << currentRange << " L0= " << lambda0
         << " E= " << currentKinEnergy << "  " 
         << particle->GetParticleName() << G4endl;
  */
  G4double tau = tPathLength/lambda0 ;

  if ((tau <= tausmall) || insideskin) {
    zPathLength = min(tPathLength, lambda0); 

  } else  if (tPathLength < currentRange*dtrl) {
    if(tau < taulim) zPathLength = tPathLength*(1.-0.5*tau) ;
    else             zPathLength = lambda0*(1.-G4Exp(-tau));

  } else if(currentKinEnergy < mass || tPathLength == currentRange)  {
    par1 = 1./currentRange ;
    par2 = 1./(par1*lambda0) ;
    par3 = 1.+par2 ;
    if(tPathLength < currentRange) {
      zPathLength = 
        (1.-G4Exp(par3*G4Log(1.-tPathLength/currentRange)))/(par1*par3);
    } else {
      zPathLength = 1./(par1*par3);
    }

  } else {
    G4double rfin = max(currentRange-tPathLength, 0.01*currentRange);
    G4double T1 = GetEnergy(particle,rfin,couple);
    G4double lambda1 = GetTransportMeanFreePath(particle,T1);

    par1 = (lambda0-lambda1)/(lambda0*tPathLength);
    //G4cout << "par1= " << par1 << " L1= " << lambda1 << G4endl;
    par2 = 1./(par1*lambda0);
    par3 = 1.+par2 ;
    zPathLength = (1.-G4Exp(par3*G4Log(lambda1/lambda0)))/(par1*par3);
  }

  zPathLength = min(zPathLength, lambda0);
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
        if(par1*par3*geomStepLength < 1.) {
          tlength = (1.-G4Exp(G4Log(1.-par1*par3*geomStepLength)/par3))/par1 ;
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
  G4double kineticEnergy = currentKinEnergy;
  if (tPathLength > currentRange*dtrl) {
    kineticEnergy = GetEnergy(particle,currentRange-tPathLength,couple);
  } else {
    kineticEnergy -= tPathLength*GetDEDX(particle,currentKinEnergy,couple);
  }

  if((kineticEnergy <= eV) || (tPathLength <= tlimitminfix) ||
     (tPathLength < tausmall*lambda0)) { return fDisplacement; }

  G4double cth = SampleCosineTheta(tPathLength,kineticEnergy);

  // protection against 'bad' cth values
  if(std::abs(cth) >= 1.0) { return fDisplacement; } 

  /*
  if(cth < 1.0 - 1000*tPathLength/lambda0 && cth < 0.5 &&
     kineticEnergy > 20*MeV) { 
    G4cout << "### G4UrbanMscModel::SampleScattering for "
           << particle->GetParticleName()
           << " E(MeV)= " << kineticEnergy/MeV
           << " Step(mm)= " << tPathLength/mm
           << " in " << CurrentCouple()->GetMaterial()->GetName()
           << " CosTheta= " << cth << G4endl;
  }
  */
  G4double sth  = sqrt((1.0 - cth)*(1.0 + cth));
  G4double phi  = twopi*rndmEngineMod->flat(); 
  G4double dirx = sth*cos(phi);
  G4double diry = sth*sin(phi);

  G4ThreeVector newDirection(dirx,diry,cth);
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
                                            G4double KineticEnergy)
{
  G4double cth = 1. ;
  G4double tau = trueStepLength/lambda0;
  currentTau   = tau;
  lambdaeff    = lambda0;

  G4double lambda1 = GetTransportMeanFreePath(particle,KineticEnergy);
  if(std::abs(lambda1 - lambda0) > lambda0*0.01 && lambda1 > 0.)
  {
    // mean tau value
    tau = trueStepLength*G4Log(lambda0/lambda1)/(lambda0-lambda1);
  }

  currentTau = tau ;
  lambdaeff = trueStepLength/currentTau;
  currentRadLength = couple->GetMaterial()->GetRadlen();

  if (tau >= taubig) { cth = -1.+2.*rndmEngineMod->flat(); }
  else if (tau >= tausmall) {
    static const G4double numlim = 0.01;
    G4double xmeanth, x2meanth;
    if(tau < numlim) {
      xmeanth = 1.0 - tau*(1.0 - 0.5*tau);
      x2meanth= 1.0 - tau*(5.0 - 6.25*tau)/3.;
    } else {
      xmeanth = G4Exp(-tau);
      x2meanth = (1.+2.*G4Exp(-2.5*tau))/3.;
    }

    // too large step of low-energy particle
    G4double relloss = 1. - KineticEnergy/currentKinEnergy;
    static const G4double rellossmax= 0.50;
    if(relloss > rellossmax) {
      return SimpleScattering(xmeanth,x2meanth);
    }
    // is step extreme small ?
    G4bool extremesmallstep = false ;
    G4double tsmall = std::min(tlimitmin,lambdalimit);
    G4double theta0 = 0.;
    if(trueStepLength > tsmall) {
      theta0 = ComputeTheta0(trueStepLength,KineticEnergy);
    } else {
      theta0 = sqrt(trueStepLength/tsmall)*ComputeTheta0(tsmall,KineticEnergy);
      extremesmallstep = true ;
    }

    static const G4double theta0max = CLHEP::pi/6.;
    //G4cout << "Theta0= " << theta0 << " theta0max= " << theta0max 
    //             << "  sqrt(tausmall)= " << sqrt(tausmall) << G4endl;

    // protection for very small angles
    G4double theta2 = theta0*theta0;

    if(theta2 < tausmall) { return cth; }
    
    if(theta0 > theta0max) {
      return SimpleScattering(xmeanth,x2meanth);
    }

    G4double x = theta2*(1.0 - theta2/12.);
    if(theta2 > numlim) {
      G4double sth = 2*sin(0.5*theta0);
      x = sth*sth;
    }

    // parameter for tail
    G4double ltau= G4Log(tau);
    G4double u   = G4Exp(ltau/6.);
    if(extremesmallstep) { u = G4Exp(G4Log(tsmall/lambda0)/6.); }
    G4double xx  = G4Log(lambdaeff/currentRadLength);
    G4double xsi = coeffc1+u*(coeffc2+coeffc3*u)+coeffc4*xx;

    // tail should not be too big
    if(xsi < 1.9) { 
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
      xsi = 1.9; 
    }

    G4double c = xsi;

    if(std::abs(c-3.) < 0.001)      { c = 3.001; }
    else if(std::abs(c-2.) < 0.001) { c = 2.001; }

    G4double c1 = c-1.;

    G4double ea = G4Exp(-xsi);
    G4double eaa = 1.-ea ;
    G4double xmean1 = 1.-(1.-(1.+xsi)*ea)*x/eaa;
    G4double x0 = 1. - xsi*x;

    // G4cout << " xmean1= " << xmean1 << "  xmeanth= " << xmeanth << G4endl;

    if(xmean1 <= 0.999*xmeanth) {
      return SimpleScattering(xmeanth,x2meanth);
    }
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
    if(rndmEngineMod->flat() < qprob)
    {
      G4double var = 0;
      if(rndmEngineMod->flat() < prob) {
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
      /*
      if(KineticEnergy > 5*GeV && cth < 0.9) {
        G4cout << "G4UrbanMscModel::SampleCosineTheta: E(GeV)= " 
               << KineticEnergy/GeV 
               << " 1-cosT= " << 1 - cth
               << " length(mm)= " << trueStepLength << " Zeff= " << Zeff 
               << " tau= " << tau
               << " prob= " << prob << " var= " << var << G4endl;
        G4cout << "  c= " << c << " qprob= " << qprob << " eb1= " << eb1
               << " ebx= " << ebx
               << " c1= " << c1 << " b= " << b << " b1= " << b1 
               << " bx= " << bx << " d= " << d
               << " ea= " << ea << " eaa= " << eaa << G4endl;
      }
      */
    }
    else {
      cth = -1.+2.*rndmEngineMod->flat();
      /*
      if(KineticEnergy > 5*GeV) {
        G4cout << "G4UrbanMscModel::SampleCosineTheta: E(GeV)= " 
               << KineticEnergy/GeV 
               << " length(mm)= " << trueStepLength << " Zeff= " << Zeff 
               << " qprob= " << qprob << G4endl;
      }
      */
    }
  }
  return cth ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UrbanMscModel::ComputeTheta0(G4double trueStepLength,
                                        G4double KineticEnergy)
{
  // for all particles take the width of the central part
  //  from a  parametrization similar to the Highland formula
  // ( Highland formula: Particle Physics Booklet, July 2002, eq. 26.10)
  G4double invbetacp = std::sqrt((currentKinEnergy+mass)*(KineticEnergy+mass)/
                                 (currentKinEnergy*(currentKinEnergy+2.*mass)*
                                  KineticEnergy*(KineticEnergy+2.*mass)));
  G4double y = trueStepLength/currentRadLength;

  if(particle == positron)
  {
    static const G4double xl= 0.6;
    static const G4double xh= 0.9;
    static const G4double e = 113.0;
    G4double corr;

    G4double tau = std::sqrt(currentKinEnergy*KineticEnergy)/mass;
    G4double x = std::sqrt(tau*(tau+2.)/((tau+1.)*(tau+1.)));
    G4double a = 0.994-4.08e-3*Zeff;
    G4double b = 7.16+(52.6+365./Zeff)/Zeff;
    G4double c = 1.000-4.47e-3*Zeff;
    G4double d = 1.21e-3*Zeff;
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
    y *= corr*(1.+Zeff*(1.84035e-4*Zeff-1.86427e-2)+0.41125);
  }

  static const G4double c_highland = 13.6*CLHEP::MeV;
  G4double theta0 = c_highland*std::abs(charge)*std::sqrt(y)*invbetacp;
 
  // correction factor from e- scattering data
  theta0 *= (coeffth1+coeffth2*G4Log(y));
  return theta0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UrbanMscModel::SampleDisplacement(G4double sth, G4double phi)
{
  G4double rmax = sqrt((tPathLength-zPathLength)*(tPathLength+zPathLength));

  static const G4double third  = 1./3.;
  G4double r = rmax*G4Exp(G4Log(rndmEngineMod->flat())*third);
  /*    
    G4cout << "G4UrbanMscModel::SampleSecondaries: e(MeV)= " << kineticEnergy
           << " sinTheta= " << sth << " r(mm)= " << r
           << " trueStep(mm)= " << tPathLength
           << " geomStep(mm)= " << zPathLength
           << G4endl;
  */

  if(r > 0.) {
    static const G4double kappa = 2.5;
    static const G4double kappami1 = 1.5;
  
    G4double latcorr = 0.;
    if((currentTau >= tausmall) && !insideskin) {
      if(currentTau < taulim) {
	latcorr = lambdaeff*kappa*currentTau*currentTau*
	  (1.-(kappa+1.)*currentTau*third)*third;

      } else {
	G4double etau = (currentTau < taubig) ? G4Exp(-currentTau) : 0.;
	latcorr = -kappa*currentTau;
	latcorr = G4Exp(latcorr)/kappami1;
	latcorr += 1.-kappa*etau/kappami1 ;
	latcorr *= 2.*lambdaeff*third;
      }
    }
    latcorr = std::min(latcorr, r);

    // sample direction of lateral displacement
    // compute it from the lateral correlation
    G4double Phi;
    if(std::abs(r*sth) < latcorr) {
      Phi  = twopi*rndmEngineMod->flat();

    } else {
      //G4cout << "latcorr= " << latcorr << "  r*sth= " << r*sth 
      //     << " ratio= " << latcorr/(r*sth) <<  G4endl;
      G4double psi = std::acos(latcorr/(r*sth));
      G4double rdm = rndmEngineMod->flat();
      Phi = (rdm < 0.5) ? phi+psi : phi-psi;
    }
    fDisplacement.set(r*std::cos(Phi),r*std::sin(Phi),0.0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
   
void G4UrbanMscModel::SampleDisplacementNew(G4double , G4double phi)
{
  //sample displacement r

  G4double rmax = sqrt((tPathLength-zPathLength)*(tPathLength+zPathLength));
  // u = (r/rmax)**2 , v=1-u
  // paramerization from ss simulation
  // f(u) = p0*exp(p1*log(v)-p2*v)+v*(p3+p4*v) 
  G4double u ,v , rej;
  G4int count = 0;

  static const G4double reps = 1.e-6;
  static const G4double  rp0 = 2.2747e+4;
  static const G4double  rp1 = 4.5980e+0;
  static const G4double  rp2 = 1.5580e+1;
  static const G4double  rp3 = 7.1287e-1;
  static const G4double  rp4 =-5.7069e-1;

  do {
    u = reps+(1.-2.*reps)*rndmEngineMod->flat();
    v = 1.-u ;
    rej = rp0*G4Exp(rp1*G4Log(v)-rp2*v) + v*(rp3+rp4*v);
  }
  // Loop checking, 15-Sept-2015, Vladimir Ivanchenko
  while (rndmEngineMod->flat() > rej && ++count < 1000);
  G4double r = rmax*sqrt(u);

  if(r > 0.)
  {
    // sample Phi using lateral correlation
    // v = Phi-phi = acos(latcorr/(r*sth))
    // v has a universal distribution which can be parametrized from ss
    // simulation as
    // f(v) = 1.49e-2*exp(-v**2/(2*0.320))+2.50e-2*exp(-31.0*log(1.+6.30e-2*v))+
    //        1.96e-5*exp(8.42e-1*log(1.+1.45e1*v))
    static const G4double probv1 = 0.305533;
    static const G4double probv2 = 0.955176;
    static const G4double vhigh  = 3.15;     
    static const G4double w2v = 1./G4Exp(30.*G4Log(1. + 6.30e-2*vhigh));
    static const G4double w3v = 1./G4Exp(-1.842*G4Log(1. + 1.45e1*vhigh));

    G4double Phi;
    G4double random = rndmEngineMod->flat();
    if(random < probv1) {
      do {
	v = G4RandGauss::shoot(rndmEngineMod,0.,0.320);
      }
      // Loop checking, 15-Sept-2015, Vladimir Ivanchenko
      while (std::abs(v) >= vhigh);
      Phi = phi + v;

    } else {

      G4double rnd = rndmEngineMod->flat();
      v = (random < probv2) 
        ? (-1.+1./G4Exp(G4Log(1.-rnd*(1.-w2v))/30.))/6.30e-2
        : (-1.+1./G4Exp(G4Log(1.-rnd*(1.-w3v))/-1.842))/1.45e1;

      rnd = rndmEngineMod->flat();
      Phi = (rnd < 0.5) ? phi+v : phi-v;
    }
    fDisplacement.set(r*std::cos(Phi),r*std::sin(Phi),0.0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
