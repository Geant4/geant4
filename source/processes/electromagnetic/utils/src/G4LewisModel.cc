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
//
// GEANT4 Class file
//
//
// File name:   G4LewisModel
//
// Author:      Laszlo Urban
//
// Creation date: 03.03.2001
//
// Modifications:
//
// 27-03-03 Move model part from G4MultipleScattering (V.Ivanchenko)
//

// Class Description:
//
// Implementation of Lewis model of multiple scattering
// H.W.Lewis Phys Rev 78 (1950) 526

// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LewisModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
#include "G4PhysicsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LewisModel::G4LewisModel(G4double& m_dtrl, G4double& m_NuclCorrPar, 
                           G4double& m_FactPar, G4double& m_facxsi, 
			   G4bool& m_samplez, const G4String& nam)
  : G4VEmModel(nam),
  highKinEnergy(10.0*TeV),
  lowKinEnergy(0.1*keV),
  taubig(8.0),
  tausmall(1.e-14),
  taulim(1.e-5),
  dtrl(m_dtrl),
  NuclCorrPar (m_NuclCorrPar),
  FactPar(m_FactPar),
  facxsi(m_facxsi),
  samplez(m_samplez),
  stepmin(1.e-6*mm)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LewisModel::~G4LewisModel()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4LewisModel::IsInCharge(const G4ParticleDefinition* p)
{
  return (p->GetPDGCharge() != 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LewisModel::Initialise(const G4ParticleDefinition* p,
                              const G4DataVector&)
{
  // set values of some data members
  if((p->GetParticleName() == "e-") || (p->GetParticleName() == "e+"))
  {
    // parameters for e+/e-
    alfa1 = 1.45 ;
    alfa2 = 0.60 ;
    alfa3 = 0.30 ;
    b = 1. ;
    xsi = facxsi*2.22 ;
    c0 = 2.30 ;
  }
  else
  {
    // parameters for heavy particles
    alfa1 = 1.10 ;
    alfa2 = 0.14 ;
    alfa3 = 0.07 ;
    b = 1. ;
    xsi = facxsi*2.70 ;
    c0 = 1.40 ;
  }
  sigmafactor = twopi*classic_electr_radius*classic_electr_radius;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LewisModel::CrossSection(const G4Material* material,
                                    const G4ParticleDefinition* p,
                                          G4double kineticEnergy,
                                          G4double,
                                          G4double)
{
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* NbOfAtomsPerVolume = material->GetVecNbOfAtomsPerVolume();
  G4int NumberOfElements = material->GetNumberOfElements();

  // loop for element in the material
  G4double sigma = 0.0;

  for (G4int iel=0; iel<NumberOfElements; iel++)
  {
    G4double atomicNumber = (*theElementVector)[iel]->GetZ();
    G4double atomicWeight = (*theElementVector)[iel]->GetA();
    sigma += NbOfAtomsPerVolume[iel]*ComputeTransportCrossSection(p,
             kineticEnergy,atomicNumber,atomicWeight);
  }
  sigma *= sigmafactor;

  // Calculate lambda
  if ( sigma > 0.0) sigma = 1.0/sigma;
  else              sigma = DBL_MAX;

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LewisModel::ComputeTransportCrossSection(
                             const G4ParticleDefinition* particle,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber,
                                   G4double AtomicWeight)
{
  const G4double epsfactor = 2.*electron_mass_c2*electron_mass_c2*
                             Bohr_radius*Bohr_radius/(hbarc*hbarc);
  const G4double epsmin = 1.e-4 , epsmax = 1.e10;

  const G4double Zdat[15] = { 4., 6.,13.,20.,26.,29.,32.,38.,47.,
                             50.,56.,64.,74.,79.,82. };

  const G4double Tdat[23] = {0.0001*MeV,0.0002*MeV,0.0004*MeV,0.0007*MeV,
                             0.001*MeV,0.002*MeV,0.004*MeV,0.007*MeV,
                             0.01*MeV,0.02*MeV,0.04*MeV,0.07*MeV,
                             0.1*MeV,0.2*MeV,0.4*MeV,0.7*MeV,
                             1.*MeV,2.*MeV,4.*MeV,7.*MeV,10.*MeV,20.*MeV,
                             10000.0*MeV};

 // corr. factors for e-/e+ lambda

          G4double celectron[15][23] =
          {{1.125,1.072,1.051,1.047,1.047,1.050,1.052,1.054,
            1.054,1.057,1.062,1.069,1.075,1.090,1.105,1.111,
            1.112,1.108,1.100,1.093,1.089,1.087,0.7235     },
           {1.408,1.246,1.143,1.096,1.077,1.059,1.053,1.051,
            1.052,1.053,1.058,1.065,1.072,1.087,1.101,1.108,
            1.109,1.105,1.097,1.090,1.086,1.082,0.7925     },
           {2.833,2.268,1.861,1.612,1.486,1.309,1.204,1.156,
            1.136,1.114,1.106,1.106,1.109,1.119,1.129,1.132,
            1.131,1.124,1.113,1.104,1.099,1.098,0.9147     },
           {3.879,3.016,2.380,2.007,1.818,1.535,1.340,1.236,
            1.190,1.133,1.107,1.099,1.098,1.103,1.110,1.113,
            1.112,1.105,1.096,1.089,1.085,1.098,0.9700     },
           {6.937,4.330,2.886,2.256,1.987,1.628,1.395,1.265,
            1.203,1.122,1.080,1.065,1.061,1.063,1.070,1.073,
            1.073,1.070,1.064,1.059,1.056,1.056,1.0022     },
           {9.616,5.708,3.424,2.551,2.204,1.762,1.485,1.330,
            1.256,1.155,1.099,1.077,1.070,1.068,1.072,1.074,
            1.074,1.070,1.063,1.059,1.056,1.052,1.0158     },
           {11.72,6.364,3.811,2.806,2.401,1.884,1.564,1.386,
            1.300,1.180,1.112,1.082,1.073,1.066,1.068,1.069,
            1.068,1.064,1.059,1.054,1.051,1.050,1.0284     },
           {18.08,8.601,4.569,3.183,2.662,2.025,1.646,1.439,
            1.339,1.195,1.108,1.068,1.053,1.040,1.039,1.039,
            1.039,1.037,1.034,1.031,1.030,1.036,1.0515     },
           {18.22,10.48,5.333,3.713,3.115,2.367,1.898,1.631,
            1.498,1.301,1.171,1.105,1.077,1.048,1.036,1.033,
            1.031,1.028,1.024,1.022,1.021,1.024,1.0834     },
           {14.14,10.65,5.710,3.929,3.266,2.453,1.951,1.669,
            1.528,1.319,1.178,1.106,1.075,1.040,1.027,1.022,
            1.020,1.017,1.015,1.013,1.013,1.020,1.0937     },
           {14.11,11.73,6.312,4.240,3.478,2.566,2.022,1.720,
            1.569,1.342,1.186,1.102,1.065,1.022,1.003,0.997,
            0.995,0.993,0.993,0.993,0.993,1.011,1.1140     },
           {22.76,20.01,8.835,5.287,4.144,2.901,2.219,1.855,
            1.677,1.410,1.224,1.121,1.073,1.014,0.986,0.976,
            0.974,0.972,0.973,0.974,0.975,0.987,1.1410     },
           {50.77,40.85,14.13,7.184,5.284,3.435,2.520,2.059,
            1.837,1.512,1.283,1.153,1.091,1.010,0.969,0.954,
            0.950,0.947,0.949,0.952,0.954,0.963,1.1750     },
           {65.87,59.06,15.87,7.570,5.567,3.650,2.682,2.182,
            1.939,1.579,1.325,1.178,1.108,1.014,0.965,0.947,
            0.941,0.938,0.940,0.944,0.946,0.954,1.1922     },
          // {45.60,47.34,15.92,7.810,5.755,3.767,2.760,2.239, // paper.....
           {55.60,47.34,15.92,7.810,5.755,3.767,2.760,2.239,
            1.985,1.609,1.343,1.188,1.113,1.013,0.960,0.939,
            0.933,0.930,0.933,0.936,0.939,0.949,1.2026     }};
           G4double cpositron[15][23] = {
           {2.589,2.044,1.658,1.446,1.347,1.217,1.144,1.110,
            1.097,1.083,1.080,1.086,1.092,1.108,1.123,1.131,
            1.131,1.126,1.117,1.108,1.103,1.100,0.7235     },
           {3.904,2.794,2.079,1.710,1.543,1.325,1.202,1.145,
            1.122,1.096,1.089,1.092,1.098,1.114,1.130,1.137,
            1.138,1.132,1.122,1.113,1.108,1.102,0.7925     },
           {7.970,6.080,4.442,3.398,2.872,2.127,1.672,1.451,
            1.357,1.246,1.194,1.179,1.178,1.188,1.201,1.205,
            1.203,1.190,1.173,1.159,1.151,1.145,0.9147     },
           {9.714,7.607,5.747,4.493,3.815,2.777,2.079,1.715,
            1.553,1.353,1.253,1.219,1.211,1.214,1.225,1.228,
            1.225,1.210,1.191,1.175,1.166,1.174,0.9700     },
           {17.97,12.95,8.628,6.065,4.849,3.222,2.275,1.820,
            1.624,1.382,1.259,1.214,1.202,1.202,1.214,1.219,
            1.217,1.203,1.184,1.169,1.160,1.151,1.0022     },
           {24.83,17.06,10.84,7.355,5.767,3.707,2.546,1.996,
            1.759,1.465,1.311,1.252,1.234,1.228,1.238,1.241,
            1.237,1.222,1.201,1.184,1.174,1.159,1.0158     },
           {23.26,17.15,11.52,8.049,6.375,4.114,2.792,2.155,
            1.880,1.535,1.353,1.281,1.258,1.247,1.254,1.256,
            1.252,1.234,1.212,1.194,1.183,1.170,1.0284     },
           {22.33,18.01,12.86,9.212,7.336,4.702,3.117,2.348,
            2.015,1.602,1.385,1.297,1.268,1.251,1.256,1.258,
            1.254,1.237,1.214,1.195,1.185,1.179,1.0515     },
           {33.91,24.13,15.71,10.80,8.507,5.467,3.692,2.808,
            2.407,1.873,1.564,1.425,1.374,1.330,1.324,1.320,
            1.312,1.288,1.258,1.235,1.221,1.205,1.0834     },
           {32.14,24.11,16.30,11.40,9.015,5.782,3.868,2.917,
            2.490,1.925,1.596,1.447,1.391,1.342,1.332,1.327,
            1.320,1.294,1.264,1.240,1.226,1.214,1.0937     },
           {29.51,24.07,17.19,12.28,9.766,6.238,4.112,3.066,
            2.602,1.995,1.641,1.477,1.414,1.356,1.342,1.336,
            1.328,1.302,1.270,1.245,1.231,1.233,1.1140     },
           {38.19,30.85,21.76,15.35,12.07,7.521,4.812,3.498,
            2.926,2.188,1.763,1.563,1.484,1.405,1.382,1.371,
            1.361,1.330,1.294,1.267,1.251,1.239,1.1410     },
           {49.71,39.80,27.96,19.63,15.36,9.407,5.863,4.155,
            3.417,2.478,1.944,1.692,1.589,1.480,1.441,1.423,
            1.409,1.372,1.330,1.298,1.280,1.258,1.1750     },
           {59.25,45.08,30.36,20.83,16.15,9.834,6.166,4.407,
            3.641,2.648,2.064,1.779,1.661,1.531,1.482,1.459,
            1.442,1.400,1.354,1.319,1.299,1.272,1.1922     },
           {56.38,44.29,30.50,21.18,16.51,10.11,6.354,4.542,
            3.752,2.724,2.116,1.817,1.692,1.554,1.499,1.474,
            1.456,1.412,1.364,1.328,1.307,1.282,1.2026     }};

   G4double sigma;

   G4double Z23 = 2.*log(AtomicNumber)/3.; Z23 = exp(Z23);

   G4double ParticleMass = particle->GetPDGMass();

  // correction if particle .ne. e-/e+
  // compute equivalent kinetic energy
  // lambda depends on p*beta ....

   if((particle->GetParticleName() != "e-") &&
      (particle->GetParticleName() != "e+") )
   {
     G4double TAU = KineticEnergy/ParticleMass ;
     G4double c = ParticleMass*TAU*(TAU+2.)/(electron_mass_c2*(TAU+1.)) ;
     G4double w = c-2. ;
     G4double tau = 0.5*(w+sqrt(w*w+4.*c)) ;
     KineticEnergy = electron_mass_c2*tau ;
   }

   G4double Charge = particle->GetPDGCharge();
   G4double ChargeSquare = Charge*Charge/(eplus*eplus);

   G4double TotalEnergy = KineticEnergy + electron_mass_c2 ;
   G4double beta2 = KineticEnergy*(TotalEnergy+electron_mass_c2)
                                 /(TotalEnergy*TotalEnergy);
   G4double bg2   = KineticEnergy*(TotalEnergy+electron_mass_c2)
                                 /(electron_mass_c2*electron_mass_c2);

   G4double eps = epsfactor*bg2/Z23;

   if     (eps<epsmin)  sigma = 2.*eps*eps;
   else if(eps<epsmax)  sigma = log(1.+2.*eps)-2.*eps/(1.+2.*eps);
   else                 sigma = log(2.*eps)-1.+1./eps;

   sigma *= ChargeSquare*AtomicNumber*AtomicNumber/(beta2*bg2);

  // nuclear size effect correction for high energy
  // ( a simple approximation at present)
  G4double corrnuclsize,a,w1,w2,w;

  G4double x0 = 1. - NuclCorrPar*ParticleMass/(KineticEnergy*
               exp(log(AtomicWeight/(g/mole))/3.));
  if ( x0 < -1. || KineticEnergy  <= 10.*MeV)
      {
        x0 = -1.;
	corrnuclsize = 1.;
      }
  else
      {
        a = 1.+1./eps;
        if (eps > epsmax) w1=log(2.*eps)+1./eps-3./(8.*eps*eps);
        else              w1=log((a+1.)/(a-1.))-2./(a+1.);
        w = 1./((1.-x0)*eps);
        if (w < epsmin)   w2=-log(w)-1.+2.*w-1.5*w*w;
        else              w2 = log((a-x0)/(a-1.))-(1.-x0)/(a-x0);
        corrnuclsize = w1/w2;
        corrnuclsize = exp(-FactPar*ParticleMass/KineticEnergy)*
                      (corrnuclsize-1.)+1.;
      }

  // interpolate in AtomicNumber and beta2
  // get bin number in Z
  G4int iZ = 14;
  while ((iZ>=0)&&(Zdat[iZ]>=AtomicNumber)) iZ -= 1;
  if (iZ==14)                               iZ = 13;
  if (iZ==-1)                               iZ = 0 ;

  G4double Z1 = Zdat[iZ];
  G4double Z2 = Zdat[iZ+1];
  G4double ratZ = (AtomicNumber-Z1)/(Z2-Z1);

  // get bin number in T (beta2)
  G4int iT = 22;
  while ((iT>=0)&&(Tdat[iT]>=KineticEnergy)) iT -= 1;
  if(iT==22)                                 iT = 21;
  if(iT==-1)                                 iT = 0 ;

  //  calculate betasquare values
  G4double T = Tdat[iT],   E = T + electron_mass_c2;
  G4double b2small = T*(E+electron_mass_c2)/(E*E);
  T = Tdat[iT+1]; E = T + electron_mass_c2;
  G4double b2big = T*(E+electron_mass_c2)/(E*E);
  G4double ratb2 = (beta2-b2small)/(b2big-b2small);

  G4double c1,c2,cc1,cc2,corr;
  if (Charge < 0.)
    {
       c1 = celectron[iZ][iT];
       c2 = celectron[iZ+1][iT];
       cc1 = c1+ratZ*(c2-c1);

       c1 = celectron[iZ][iT+1];
       c2 = celectron[iZ+1][iT+1];
       cc2 = c1+ratZ*(c2-c1);

       corr = cc1+ratb2*(cc2-cc1);
       sigma /= corr;
    }

  if (Charge > 0.)
    {
       c1 = cpositron[iZ][iT];
       c2 = cpositron[iZ+1][iT];
       cc1 = c1+ratZ*(c2-c1);

       c1 = cpositron[iZ][iT+1];
       c2 = cpositron[iZ+1][iT+1];
       cc2 = c1+ratZ*(c2-c1);

       corr = cc1+ratb2*(cc2-cc1);
       sigma /= corr;
    }

  //  nucl. size correction for particles other than e+/e- only at present !!!!
  if((particle->GetParticleName() != "e-") &&
     (particle->GetParticleName() != "e+")   )
     sigma /= corrnuclsize;

  return sigma;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LewisModel::GeomPathLength(
                          G4PhysicsTable* theLambdaTable,
                    const G4MaterialCutsCouple* couple,
		    const G4ParticleDefinition* particle,
		          G4double& T0,
			  G4double lambda,
			  G4double range,
			  G4double truePathLength)
{
  //  do the true -> geom transformation
  const G4double ztmin = 1./3., ztmax = 0.98 ;

  lambda0 = lambda;
  lambda1 = -1.;
  lambdam = -1.;
  alam    = range;
  zm      = 1.;
  cthm    = 1.;
  tPathLength = truePathLength;
  G4double tau   = tPathLength/lambda0 ;

  if (tau <= tausmall) return tPathLength;

  G4double zmean = tPathLength;

  if (tPathLength < range*dtrl) {
    zmean = lambda0*(1.-exp(-tau));

  } else {
    G4LossTableManager* theManager = G4LossTableManager::Instance();
    G4double T1 = theManager->GetEnergy(particle,range-tPathLength,couple);

    if (theLambdaTable) {
      G4bool b;
      lambda1 = ((*theLambdaTable)[couple->GetIndex()])->GetValue(T1,b);
    } else {
      lambda1 = CrossSection(couple->GetMaterial(),particle,T1,0.0,1.0);
    }
    if (T0 > particle->GetPDGMass()) alam = lambda0*tPathLength/(lambda0-lambda1) ;
    G4double blam = 1.+alam/lambda0 ;
    if (tPathLength < 2.*dtrl*range) {
      zmean = alam*(1.-exp(blam*log(1.-tPathLength/alam)))/blam ;

    } else {
      G4double w = 1.-0.5*tPathLength/alam ;
      lambdam = lambda0*w ;
      G4double clam = 1.+alam/lambdam ;
      cthm = exp(alam*log(w)/lambda0) ;
      zm = alam*(1.-exp(blam*log(w)))/blam ;
      zmean = zm + alam*(1.-exp(clam*log(w)))*cthm/clam ;
    }
  }
    //  sample z
  G4double zPathLength = zmean ;
  G4double zt = zmean/tPathLength ;
  if (tPathLength >= stepmin && samplez && zt > ztmin && zt < ztmax) {
    G4double cz = 0.5*(3.*zt-1.)/(1.-zt) ;
    G4double cz1 = 1.+cz ;
    G4double u0 = cz/cz1 ;
    G4double u,grej ;
    G4double grej0 = exp(cz*log(u0))*(1.-u0) ;
    do {
        u = exp(log(G4UniformRand())/cz1) ;
        grej  = exp(cz*log(u))*(1.-u) ;
        if (grej > grej0)
            G4cout << "G4LewisModel: Warning! majorant "
                   << grej0 << " < " << grej << G4endl;
      } while (grej < grej0*G4UniformRand()) ;
   zPathLength = tPathLength*u ;
  }

  return zPathLength ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LewisModel::TrueStepLength(G4double geomStepLength)
{
  G4double trueLength = geomStepLength;
  if (geomStepLength > lambda0*tausmall) {
    G4double blam = 1.+alam/lambda0;
//    G4cout << "alam= " << alam << " blam= " << blam << " lambda1= " << lambda1 << G4endl;
    if (lambda1 < 0.) {
      trueLength = -lambda0*log(1.-geomStepLength/lambda0) ;

    } else if (lambdam < 0.) {
      if (blam*geomStepLength < alam)
        trueLength = alam*(1.-exp(log(1.-blam*geomStepLength/alam)/blam));
      else
        trueLength = tPathLength;

    } else {
      if (geomStepLength <= zm) {

         // ??? here is a problems - redefine variable
        lambdam = -1.;

	if(blam*geomStepLength < alam )
          trueLength = alam*(1.-exp(log(1.-blam*geomStepLength/alam)/blam));
        else
          trueLength = 0.5*tPathLength;

      } else {
        G4double clam = 1.+alam/lambdam;
//    G4cout << "clam= " << clam << " zm= " << zm << " cthm= " << cthm << G4endl;
        if(clam*(geomStepLength-zm)/(alam*cthm) < 1.)
          trueLength = 0.5*tPathLength + alam*(1.-
                exp(log(1.-clam*(geomStepLength-zm)/(alam*cthm)))/clam) ;
        else
          trueLength = tPathLength;
      }
    }
//    G4cout << "tLenth= " << trueLength << " tpl= " << tPathLength << G4endl;
    if(trueLength > tPathLength) trueLength = tPathLength;
    if(trueLength < geomStepLength) trueLength = geomStepLength;
  }

  return trueLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LewisModel::SampleCosineTheta(G4double trueStepLength)
{
  G4double cth = 1.;
  currentTau = trueStepLength/lambda0;
  if(trueStepLength < stepmin)
    cth = exp(-currentTau) ;
  else
  {       
//  G4cout << "tau= " << currentTau << " lambda1= " << lambda1 << " lambdam= " << lambdam << G4endl;
    if (currentTau > taubig) cth = -1.+2.*G4UniformRand();
    else if (currentTau >= tausmall)
    {
      if(lambda1 > 0.)
      {
        if (lambdam < 0.)
          currentTau = -alam*log(1.-trueStepLength/alam)/lambda0 ;
        else
          currentTau = -log(cthm)-alam*
                       log(1.-(trueStepLength-0.5*tPathLength)/alam)/lambdam ;
      }
      if(currentTau > taubig) cth = -1.+2.*G4UniformRand();
      else
      {
        const G4double amax = 25. ;
        const G4double tau0 = 0.02  ;

        G4double a;

        G4double w = log(currentTau/tau0) ;
        if (currentTau < tau0) a = (alfa1-alfa2*w)/currentTau ;
        else            a = (alfa1+alfa3*w)/currentTau ;

        G4double x0 = 1.-xsi/a ;
        if(x0 < 0.) x0 = 0. ;

        // from continuity of the 1st derivatives
        G4double c = a*(b-x0) ;
        if(a*currentTau < c0) c = c0*(b-x0)/currentTau ;

        if(c == 1.) c=1.000001 ;
        if(c == 2.) c=2.000001 ;
        if(c == 3.) c=3.000001 ;

        G4double ea = 0.0;
        if (a*(1.-x0) < amax) ea = exp(-a*(1.-x0)) ;

        G4double eaa = 1.-ea ;
        G4double xmean1 = 1.-1./a+(1.-x0)*ea/eaa ;

        G4double b1 = b+1. ;
        G4double bx = b-x0 ;
        G4double eb1= exp((c-1.)*log(b1)) ;
        G4double ebx= exp((c-1.)*log(bx)) ;

        G4double xmean2  = (x0*eb1+ebx+(eb1*bx-b1*ebx)/(2.-c))/(eb1-ebx) ;
        G4double xmeanth = exp(-currentTau) ;

        G4double cnorm1 = a/eaa ;
        G4double cnorm2 = (c-1.)*eb1*ebx/(eb1-ebx) ;
        G4double f1x0 = cnorm1*exp(-a*(1.-x0)) ;
        G4double f2x0 = cnorm2/exp(c*log(b-x0)) ;

        // from continuity at x=x0
        G4double prob = f2x0/(f1x0+f2x0) ;
        // from xmean = xmeanth
        G4double qprob = (f1x0+f2x0)*xmeanth/(f2x0*xmean1+f1x0*xmean2) ;

        // protection against qprob > 1
        if(qprob > 1.)
        {
          qprob = 1. ;
          prob = (xmeanth-xmean2)/(xmean1-xmean2) ;
        }
/*
        G4cout << "tau= " << currentTau << " prob= " << prob << " qprob= " << qprob << G4endl;
        G4cout << "ea= " << ea << " eaa= " << eaa << " a= " << a
               << " b= " << b << " b1= " << b1 << " bx= " << bx
	       << " ebx= " << ebx << " eb1= " << eb1 << " c= " << c
	       << G4endl;
*/
        // sampling of costheta
        if (G4UniformRand() < qprob)
        {
          if (G4UniformRand() < prob)
             cth = 1.+log(ea+G4UniformRand()*eaa)/a ;
          else
             cth = b-b1*bx/exp(log(ebx-G4UniformRand()*(ebx-eb1))/(c-1.)) ;
        }
        else
        {
          cth = -1.+2.*G4UniformRand();
        }
      }
    }
  }
//   G4cout << "cth= " << cth << G4endl;
  return cth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LewisModel::SampleDisplacement()
{
  const G4double kappa = 2.5;
  const G4double kappapl1 = kappa+1.;
  const G4double kappami1 = kappa-1.;
  G4double rmean = 0.0;
  if (currentTau >= tausmall) {
    if (currentTau < taulim) {
      rmean = kappa*currentTau*currentTau*currentTau*(1.-kappapl1*currentTau*0.25)/6. ;

    } else {
      G4double etau = 0.0;
      if (currentTau<taubig) etau = exp(-currentTau);
      rmean = -kappa*currentTau;
      rmean = -exp(rmean)/(kappa*kappami1);
      rmean += currentTau-kappapl1/kappa+kappa*etau/kappami1;
   // G4cout << "tau= " << currentTau << " lambda0= " << lambda0 << " etau= " << etau << " kappa= " << kappa << G4endl;
    }
    if (rmean>0.) rmean = 2.*lambda0*sqrt(rmean/3.0);
    else          rmean = 0.;
    //     G4cout << "rmean= " << rmean << G4endl;
 }
  return rmean;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

