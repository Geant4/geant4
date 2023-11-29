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
// File name:     G4IonFluctuation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 28-12-02 add method Dispersion (V.Ivanchenko)
// 07-02-03 change signature (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 23-05-03 Add control on parthalogical cases (V.Ivanchenko)
// 16-10-03 Changed interface to Initialisation (V.Ivanchenko)
// 27-09-07 Use FermiEnergy from material, add cut dependence (V.Ivanchenko)
// 01-02-08 Add protection for small energies and optimise the code (V.Ivanchenko)
// 01-06-08 Added initialisation of effective charge prestep (V.Ivanchenko)
//
// Class Description:
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4Pow.hh"
#include "G4Log.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4IonFluctuations::G4IonFluctuations(const G4String& nam)
  : G4VEmFluctuationModel(nam),
    particleMass(CLHEP::proton_mass_c2),
    parameter(10.0*CLHEP::MeV/CLHEP::proton_mass_c2),
    theBohrBeta2(50.0*keV/CLHEP::proton_mass_c2),
    minLoss(0.001*CLHEP::eV)
{
  uniFluct = new G4UniversalFluctuation();
  g4calc = G4Pow::GetInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4IonFluctuations::~G4IonFluctuations() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonFluctuations::InitialiseMe(const G4ParticleDefinition* part)
{
  particle = part;
  particleMass = part->GetPDGMass();
  charge = part->GetPDGCharge()/eplus;
  effChargeSquare = chargeSquare = charge*charge;
  uniFluct->InitialiseMe(part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4IonFluctuations::SampleFluctuations(const G4MaterialCutsCouple* couple,
                                      const G4DynamicParticle* dp,
                                      const G4double tcut,
                                      const G4double tmax,
                                      const G4double length,
                                      const G4double meanLoss)
{
  //  G4cout << "### meanLoss= " << meanLoss << G4endl;
  if(meanLoss <= minLoss) return meanLoss;

  //G4cout << "G4IonFluctuations::SampleFluctuations E(MeV)= " 
  //     << dp->GetKineticEnergy()
  //         << "  Elim(MeV)= " << parameter*charge*particleMass << G4endl; 

  // Vavilov fluctuations above energy threshold
  if(dp->GetKineticEnergy() > parameter*charge*particleMass) {
    return uniFluct->SampleFluctuations(couple,dp,tcut,tmax,length,meanLoss);
  }

  const G4Material* material = couple->GetMaterial();
  G4double siga = Dispersion(material,dp,tcut,tmax,length);
  G4double loss = meanLoss;
  
  //G4cout << "### siga= " << sqrt(siga) << "  navr= " << navr << G4endl;

  // Gaussian fluctuation

  // Increase fluctuations for big fractional energy loss
  //G4cout << "siga= " << siga << G4endl;
  if ( meanLoss > minFraction*kineticEnergy ) {
    G4double gam = (kineticEnergy - meanLoss)/particleMass + 1.0;
    G4double b2  = 1.0 - 1.0/(gam*gam);
    if(b2 < xmin*beta2) b2 = xmin*beta2;
    G4double x   = b2/beta2;
    G4double x3  = 1.0/(x*x*x);
    siga *= 0.25*(1.0 + x)*(x3 + (1.0/b2 - 0.5)/(1.0/beta2 - 0.5) );
  }
  siga = std::sqrt(siga);
  G4double sn = meanLoss/siga;
  G4double twomeanLoss = meanLoss + meanLoss;
  //  G4cout << "siga= " << siga << "  sn= " << sn << G4endl;
  
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  // thick target case  
  if (sn >= 2.0) {

    do {
      loss = G4RandGauss::shoot(rndmEngine,meanLoss,siga);
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while (0.0 > loss || twomeanLoss < loss);

    // Gamma distribution
  } else if(sn > 0.1) {

    G4double neff = sn*sn;
    loss = meanLoss*G4RandGamma::shoot(rndmEngine,neff,1.0)/neff;

    // uniform distribution for very small steps
  } else {
    loss = twomeanLoss*rndmEngine->flat();
  }

  //G4cout << "meanLoss= " << meanLoss << " loss= " << loss << G4endl;
  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonFluctuations::Dispersion(const G4Material* material,
                                       const G4DynamicParticle* dp,
                                       const G4double tcut,
                                       const G4double tmax,
                                       const G4double length)
{
  if(dp->GetDefinition() != particle) { InitialiseMe(dp->GetDefinition()); }

  const G4double beta = dp->GetBeta();
  kineticEnergy = dp->GetKineticEnergy();
  beta2 = beta*beta;

  G4double siga = (tmax/beta2 - 0.5*tcut)*CLHEP::twopi_mc2_rcl2*length*
    material->GetElectronDensity()*effChargeSquare;

  // Low velocity - additional ion charge fluctuations according to
  // Q.Yang et al., NIM B61(1991)149-155.
  //G4cout << "sigE= " << sqrt(siga) << " charge= " << charge <<G4endl;

  G4double Z = material->GetIonisation()->GetZeffective();
  G4double fac = Factor(material, Z);

  // heavy ion correction
  //  G4double f1 = 1.065e-4*chargeSquare;
  //  if(beta2 > theBohrBeta2)  f1/= beta2;
  //  else                      f1/= theBohrBeta2;
  //  if(f1 > 2.5) f1 = 2.5;
  //  fac *= (1.0 + f1);

  // taking into account the cutg
  G4double fac_cut = 1.0 + (fac - 1.0)*2.0*CLHEP::electron_mass_c2*beta2
    /(tmax*(1.0 - beta2));
  if(fac_cut > 0.01 && fac > 0.01) {
    siga *= fac_cut;
  }
  /*
  G4cout << "siga(keV)= " << sqrt(siga)/keV << " fac= " << fac 
         << "  f1= " << fac_cut << G4endl;
  */
  return siga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonFluctuations::Factor(const G4Material* material, G4double Z)
{
  // The aproximation of energy loss fluctuations
  // Q.Yang et al., NIM B61(1991)149-155.

  // Reduced energy in MeV/AMU
  G4double energy = kineticEnergy*CLHEP::amu_c2/(particleMass*CLHEP::MeV);

  // simple approximation for higher beta2
  G4double s1 = RelativisticFactor(material, Z);

  // tabulation for lower beta2
  if( beta2 < 3.0*theBohrBeta2*Z ) {

  static const G4double a[96][4] = {
 {-0.3291, -0.8312,  0.2460, -1.0220},
 {-0.5615, -0.5898,  0.5205, -0.7258},
 {-0.5280, -0.4981,  0.5519, -0.5865},
 {-0.5125, -0.4625,  0.5660, -0.5190},
 {-0.5127, -0.8595,  0.5626, -0.8721},
 {-0.5174, -1.1930,  0.5565, -1.1980},
 {-0.5179, -1.1850,  0.5560, -1.2070},
 {-0.5209, -0.9355,  0.5590, -1.0250},
 {-0.5255, -0.7766,  0.5720, -0.9412},

 {-0.5776, -0.6665,  0.6598, -0.8484},
 {-0.6013, -0.6045,  0.7321, -0.7671},
 {-0.5781, -0.5518,  0.7605, -0.6919},
 {-0.5587, -0.4981,  0.7835, -0.6195},
 {-0.5466, -0.4656,  0.7978, -0.5771},
 {-0.5406, -0.4690,  0.8031, -0.5718},
 {-0.5391, -0.5061,  0.8024, -0.5974},
 {-0.5380, -0.6483,  0.7962, -0.6970},
 {-0.5355, -0.7722,  0.7962, -0.7839},
 {-0.5329, -0.7720,  0.7988, -0.7846},

 {-0.5335, -0.7671,  0.7984, -0.7933},
 {-0.5324, -0.7612,  0.7998, -0.8031},
 {-0.5305, -0.7300,  0.8031, -0.7990},
 {-0.5307, -0.7178,  0.8049, -0.8216},
 {-0.5248, -0.6621,  0.8165, -0.7919},
 {-0.5180, -0.6502,  0.8266, -0.7986},
 {-0.5084, -0.6408,  0.8396, -0.8048},
 {-0.4967, -0.6331,  0.8549, -0.8093},
 {-0.4861, -0.6508,  0.8712, -0.8432},
 {-0.4700, -0.6186,  0.8961, -0.8132},

 {-0.4545, -0.5720,  0.9227, -0.7710},
 {-0.4404, -0.5226,  0.9481, -0.7254},
 {-0.4288, -0.4778,  0.9701, -0.6850},
 {-0.4199, -0.4425,  0.9874, -0.6539},
 {-0.4131, -0.4188,  0.9998, -0.6332},
 {-0.4089, -0.4057,  1.0070, -0.6218},
 {-0.4039, -0.3913,  1.0150, -0.6107},
 {-0.3987, -0.3698,  1.0240, -0.5938},
 {-0.3977, -0.3608,  1.0260, -0.5852},
 {-0.3972, -0.3600,  1.0260, -0.5842},

 {-0.3985, -0.3803,  1.0200, -0.6013},
 {-0.3985, -0.3979,  1.0150, -0.6168},
 {-0.3968, -0.3990,  1.0160, -0.6195},
 {-0.3971, -0.4432,  1.0050, -0.6591},
 {-0.3944, -0.4665,  1.0010, -0.6825},
 {-0.3924, -0.5109,  0.9921, -0.7235},
 {-0.3882, -0.5158,  0.9947, -0.7343},
 {-0.3838, -0.5125,  0.9999, -0.7370},
 {-0.3786, -0.4976,  1.0090, -0.7310},
 {-0.3741, -0.4738,  1.0200, -0.7155},

 {-0.3969, -0.4496,  1.0320, -0.6982},
 {-0.3663, -0.4297,  1.0430, -0.6828},
 {-0.3630, -0.4120,  1.0530, -0.6689},
 {-0.3597, -0.3964,  1.0620, -0.6564},
 {-0.3555, -0.3809,  1.0720, -0.6454},
 {-0.3525, -0.3607,  1.0820, -0.6289},
 {-0.3505, -0.3465,  1.0900, -0.6171},
 {-0.3397, -0.3570,  1.1020, -0.6384},
 {-0.3314, -0.3552,  1.1130, -0.6441},
 {-0.3235, -0.3531,  1.1230, -0.6498},

 {-0.3150, -0.3483,  1.1360, -0.6539},
 {-0.3060, -0.3441,  1.1490, -0.6593},
 {-0.2968, -0.3396,  1.1630, -0.6649},
 {-0.2935, -0.3225,  1.1760, -0.6527},
 {-0.2797, -0.3262,  1.1940, -0.6722},
 {-0.2704, -0.3202,  1.2100, -0.6770},
 {-0.2815, -0.3227,  1.2480, -0.6775},
 {-0.2880, -0.3245,  1.2810, -0.6801},
 {-0.3034, -0.3263,  1.3270, -0.6778},
 {-0.2936, -0.3215,  1.3430, -0.6835},

 {-0.3282, -0.3200,  1.3980, -0.6650},
 {-0.3260, -0.3070,  1.4090, -0.6552},
 {-0.3511, -0.3074,  1.4470, -0.6442},
 {-0.3501, -0.3064,  1.4500, -0.6442},
 {-0.3490, -0.3027,  1.4550, -0.6418},
 {-0.3487, -0.3048,  1.4570, -0.6447},
 {-0.3478, -0.3074,  1.4600, -0.6483},
 {-0.3501, -0.3283,  1.4540, -0.6669},
 {-0.3494, -0.3373,  1.4550, -0.6765},
 {-0.3485, -0.3373,  1.4570, -0.6774},

 {-0.3462, -0.3300,  1.4630, -0.6728},
 {-0.3462, -0.3225,  1.4690, -0.6662},
 {-0.3453, -0.3094,  1.4790, -0.6553},
 {-0.3844, -0.3134,  1.5240, -0.6412},
 {-0.3848, -0.3018,  1.5310, -0.6303},
 {-0.3862, -0.2955,  1.5360, -0.6237},
 {-0.4262, -0.2991,  1.5860, -0.6115},
 {-0.4278, -0.2910,  1.5900, -0.6029},
 {-0.4303, -0.2817,  1.5940, -0.5927},
 {-0.4315, -0.2719,  1.6010, -0.5829},

 {-0.4359, -0.2914,  1.6050, -0.6010},
 {-0.4365, -0.2982,  1.6080, -0.6080},
 {-0.4253, -0.3037,  1.6120, -0.6150},
 {-0.4335, -0.3245,  1.6160, -0.6377},
 {-0.4307, -0.3292,  1.6210, -0.6447},
 {-0.4284, -0.3204,  1.6290, -0.6380},
 {-0.4227, -0.3217,  1.6360, -0.6438}
    } ;

    G4int iz = G4lrint(Z) - 2;
    if( 0 > iz )      { iz = 0; }
    else if(95 < iz ) { iz = 95; }

    const G4double ss = 1.0 + a[iz][0]*g4calc->powA(energy,a[iz][1])+
      + a[iz][2]*g4calc->powA(energy,a[iz][3]);
  
    // protection for the validity range for low beta
    static const G4double slim = 0.001;
    if(ss < slim) { s1 = 1.0/slim; }
    // for high value of beta
    else if(s1*ss < 1.0) { s1 = 1.0/ss; }
  }
  G4int i = 0 ;
  G4double factor = 1.0 ;

  // The index of set of parameters i = 0 for protons(hadrons) in gases
  //                                    1 for protons(hadrons) in solids
  //                                    2 for ions in atomic gases
  //                                    3 for ions in molecular gases
  //                                    4 for ions in solids
  static const G4double b[5][4] = {
  {0.1014,  0.3700,  0.9642,  3.987},
  {0.1955,  0.6941,  2.522,   1.040},
  {0.05058, 0.08975, 0.1419, 10.80},
  {0.05009, 0.08660, 0.2751,  3.787},
  {0.01273, 0.03458, 0.3951,  3.812}
  } ;

  // protons (hadrons)
  if(1.5 > charge) {
    if( kStateGas != material->GetState() ) { i = 1; }

  // ions
  } else {

    factor = charge * g4calc->A13(charge/Z);

    if( kStateGas == material->GetState() ) {
      energy /= (charge * std::sqrt(charge)) ;

      if(1 == (material->GetNumberOfElements())) {
        i = 2 ;
      } else {
        i = 3 ;
      }

    } else {
      energy /= (charge * std::sqrt(charge*Z)) ;
      i = 4 ;
    }
  }

  G4double x = b[i][2];
  G4double y = energy * b[i][3];
  if(y <= 0.2) x *= (y*(1.0 - 0.5*y));
  else         x *= (1.0 - g4calc->expA(-y));

  y = energy - b[i][1];

  const G4double s2 = factor * x * b[i][0] / (y*y + x*x);
  /*  
  G4cout << "s1= " << s1 << " s2= " << s2 << " q^2= " << effChargeSquare 
         << " e= " << energy << G4endl;
  */
  return s1*effChargeSquare/chargeSquare + s2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonFluctuations::RelativisticFactor(const G4Material* mat, 
                                               G4double Z)
{
  G4double eF = mat->GetIonisation()->GetFermiEnergy();
  G4double I  = mat->GetIonisation()->GetMeanExcitationEnergy();

  // H.Geissel et al. NIM B, 195 (2002) 3.
  G4double bF2= 2.0*eF/CLHEP::electron_mass_c2;
  G4double f  = 0.4*(1.0 - beta2)/((1.0 - 0.5*beta2)*Z);
  if(beta2 > bF2) f *= G4Log(2.0*CLHEP::electron_mass_c2*beta2/I)*bF2/beta2;
  else            f *= G4Log(4.0*eF/I);

  //  G4cout << "f= " << f << " beta2= " << beta2 
  //         << " bf2= " << bF2 << " q^2= " << chargeSquare << G4endl;

  return 1.0 + f;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonFluctuations::SetParticleAndCharge(const G4ParticleDefinition* part,
                                             G4double q2)
{
  if(part != particle) {
    particle       = part;
    particleMass   = part->GetPDGMass();
    charge         = part->GetPDGCharge()/eplus;
    chargeSquare   = charge*charge;
  }
  effChargeSquare = q2;
  uniFluct->SetParticleAndCharge(part, q2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
