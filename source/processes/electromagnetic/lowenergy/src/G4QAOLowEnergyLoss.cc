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
// -------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: New Implementation
//
//      ---------- G4QAOLowEnergyLoss physics process -------
//                  by Stephane Chauvie, 5 May 2000
// Modified:
//
// 24/05/2000 MGP  Modified to remove compilation warnings on Linux and DEC
//                 Introduced sizes of L0, L1, L2 arrays
// 23/05/2000 MGP  Made compliant to design
// 02/08/2000 V.Ivanchenko Clean up according new design
// 16/09/2000 S. Chauvie  Oscillator for all materials
// 03/10/2000 V.Ivanchenko CodeWizard clean up
// 05/11/2000 V.Ivanchenko "Aluminum" - correct name, end of cycle
//            over shells, and two bugs from previous edition
// 10/05/2001 V.Ivanchenko Clean up againist Linux compilation with -Wall
// 13/05/2001 S. Chauvie corrected bugs
// 01/06/2001 V.Ivanchenko replace names by Z, change the validity range
//                         from 50 keV to 5 KeV and change sign of the
//                         Barkas term
// 4/06/2001 S. Chauvie  Corrected small bugs
//
// ************************************************************
// It is the Quantal Harmonic Oscillator Model for energy loss
// of slow antiproton
// ************************************************************
// --------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4QAOLowEnergyLoss.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4AntiProton.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4QAOLowEnergyLoss::G4QAOLowEnergyLoss(const G4String& name)
  : G4VLowEnergyModel(name)
{
  numberOfMaterials = 6;
  sizeL0 = 67;
  sizeL1 = 22;
  sizeL2 = 14;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4QAOLowEnergyLoss::~G4QAOLowEnergyLoss()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::HighEnergyLimit(const G4ParticleDefinition* ,
					     const G4Material* ) const
{
  return 2.0*MeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::LowEnergyLimit(const G4ParticleDefinition* ,
					    const G4Material* ) const
{
  return 5.0*keV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::HighEnergyLimit(const G4ParticleDefinition* ) const
{
  return 2.0*MeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::LowEnergyLimit(const G4ParticleDefinition* ) const
{
  return 5.0*keV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4QAOLowEnergyLoss::IsInCharge(const G4DynamicParticle* particle,
				      const G4Material* material) const
{
  G4bool isInCharge = false;
  G4bool hasMaterial = false;

  if (material->GetNumberOfElements() == 1) hasMaterial = true;

  if ((particle->GetDefinition()) == (G4AntiProton::AntiProtonDefinition())
               && hasMaterial) isInCharge = true;

  return isInCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4QAOLowEnergyLoss::IsInCharge(const G4ParticleDefinition* aParticle,
				      const G4Material* material) const
{

  G4bool isInCharge = false;
  G4bool hasMaterial = false;

  if (material->GetNumberOfElements() == 1) hasMaterial = true;

  if (aParticle == (G4AntiProton::AntiProtonDefinition())
                && hasMaterial) isInCharge = true;
  return isInCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::TheValue(const G4DynamicParticle* particle,
	       	                      const G4Material* material)
{
  G4double zParticle = (G4int)(particle->GetCharge())/eplus;

  G4double energy = particle->GetKineticEnergy() ;
  G4double eloss  = EnergyLoss(material,energy,zParticle) ;

  return eloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::TheValue(const G4ParticleDefinition* aParticle,
       		                      const G4Material* material,
				      G4double kineticEnergy)
{
  G4double zParticle = (aParticle->GetPDGCharge())/eplus;

  G4double eloss  = EnergyLoss(material,kineticEnergy,zParticle) ;

  return eloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::EnergyLoss(const G4Material* material,
					G4double kineticEnergy,
					G4double zParticle) const
{
  G4int nbOfShell = GetNumberOfShell(material);
  if(nbOfShell < 1) nbOfShell = 1; 
  G4double dedx=0;

  G4double v = c_light * std::sqrt( 2.0 * kineticEnergy / proton_mass_c2 );
  G4double coeff = twopi * proton_mass_c2 *
                  (material-> GetTotNbOfElectPerVolume()) /
                   electron_mass_c2 ;
  G4double fBetheVelocity = fine_structure_const * c_light / v;
  coeff *= fine_structure_const * fine_structure_const * hbarc_squared /
           kineticEnergy ;

  G4double l0Term = 0, l1Term = 0, l2Term = 0;

  for (G4int nos = 0 ; nos < nbOfShell ; nos++){

    G4double l0 = 0, l1 = 0, l2 = 0;
    G4double NormalizedEnergy = ( 2.0 * electron_mass_c2 * v * v  ) /
                          ( c_squared * GetShellEnergy(material,nos) );

    G4double shStrength = GetShellStrength(material,nos);

    l0 = GetL0(NormalizedEnergy);
    l0Term += shStrength  * l0;

    l1 = GetL1(NormalizedEnergy);
    l1Term += shStrength * l1;

    l2 = GetL2(NormalizedEnergy);
    l2Term += shStrength * l2;
  }

  dedx = coeff * zParticle * zParticle * (l0Term
       + zParticle * fBetheVelocity * l1Term
       + zParticle * zParticle * fBetheVelocity * fBetheVelocity * l2Term);

  return dedx ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4QAOLowEnergyLoss::GetNumberOfShell(const G4Material* material) const
{
  // Set return value from table
  G4int Z = (G4int)(material->GetZ());
  G4int nShell = 0;

  // Set return value if in material available from Aahrus
  for(G4int i=0; i<numberOfMaterials; i++) {

    if(materialAvailable[i] == Z){
    	nShell = nbofShellForMaterial[i];
	break;
    }
    else 
      nShell = fNumberOfShells[Z];
  }
   return nShell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetShellEnergy(const G4Material* material,
                                            G4int nbOfTheShell) const
{
    //
    G4double shellEnergy = alShellEnergy[0];

    if(material->GetZ() == 13) shellEnergy =  alShellEnergy[nbOfTheShell];
    else if(material->GetZ() == 14)shellEnergy = siShellEnergy[nbOfTheShell];
    else if(material->GetZ() == 29)shellEnergy = cuShellEnergy[nbOfTheShell];
    else if(material->GetZ() == 73)shellEnergy =  taShellEnergy[nbOfTheShell];
    else if(material->GetZ() == 79)shellEnergy =  auShellEnergy[nbOfTheShell];
    else if(material->GetZ() == 78)shellEnergy =  ptShellEnergy[nbOfTheShell];
    else if  (material->GetNumberOfElements() == 1)
      shellEnergy = GetOscillatorEnergy(material, nbOfTheShell);
    else 
      {	
	  G4ExceptionDescription ed;	  
	  ed << "The model is not available for "
	     << material->GetName()
	     << G4endl;
	  G4Exception("G4QAOLowEnergyLoss::GetShellEnergy()",
		      "em2638",JustWarning,ed);
      }

  return  shellEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetOscillatorEnergy(const G4Material* material,
                                                 G4int nbOfTheShell) const
{

  const G4Element* element = material->GetElement(0);

  G4int Z = (G4int)(element->GetZ());

  G4double squaredPlasmonEnergy = 28.816 * 28.816  * 1e-6
				* material->GetDensity()/g/cm3
				* (Z/element->GetN()) ;

  G4double plasmonTerm = 0.66667 * GetOccupationNumber(Z,nbOfTheShell)
                       * squaredPlasmonEnergy / (Z*Z) ;

  G4double ionTerm = G4Exp(0.5) * (element->GetAtomicShell(nbOfTheShell)) ;
  ionTerm = ionTerm*ionTerm ;

  G4double oscShellEnergy = std::sqrt( ionTerm + plasmonTerm );
  return  oscShellEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetShellStrength(const G4Material* material,
                                              G4int nbOfTheShell) const
{
  G4double shellStrength = alShellStrength[0];

  if(material->GetZ() == 13) shellStrength = alShellStrength[nbOfTheShell];
  else if(material->GetZ() == 14)shellStrength =siShellStrength[nbOfTheShell];
  else if(material->GetZ() == 29)shellStrength =cuShellStrength[nbOfTheShell];
  else if(material->GetZ() == 73)shellStrength =taShellStrength[nbOfTheShell];
  else if(material->GetZ() == 79)shellStrength =auShellStrength[nbOfTheShell];
  else if(material->GetZ() == 78)shellStrength =ptShellStrength[nbOfTheShell];
  else if  (material->GetNumberOfElements() == 1){
    G4int Z = (G4int)(material->GetZ());
    shellStrength = GetOccupationNumber(Z,nbOfTheShell) / Z ;}
  else
    {
       G4ExceptionDescription ed;       
       ed << "The model is not available for "
	  << material->GetName()
	  << G4endl;
       G4Exception("G4QAOLowEnergyLoss::GetShellStrength()",
		   "em2639",JustWarning,ed);
    }
  return shellStrength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetOccupationNumber(G4int Z, G4int ShellNb) const
{

  G4int indice = ShellNb ;
  for (G4int z = 1 ; z < Z ; z++) {indice += fNumberOfShells[z];}

  return nbOfElectronPerSubShell[indice+1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetL0(G4double normEnergy) const
{
  G4int n;

  for(n = 0; n < sizeL0; n++) {
    if( normEnergy < L0[n][0] ) break;
  }
  if(0 == n) n = 1 ;
  if(n >= sizeL0) n = sizeL0 - 1 ;

  G4double l0    = L0[n][1];
  G4double l0p   = L0[n-1][1];
  G4double bethe = l0p + (l0 - l0p) * ( normEnergy - L0[n-1][0]) /
                  (L0[n][0] - L0[n-1][0]);
  return bethe ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetL1(G4double normEnergy) const
{
  G4int n;

  for(n = 0; n < sizeL1; n++) {
    if( normEnergy < L1[n][0] ) break;
  }
  if(0 == n) n = 1 ;
  if(n >= sizeL1) n = sizeL1 - 1 ;

  G4double l1    = L1[n][1];
  G4double l1p   = L1[n-1][1];
  G4double barkas= l1p + (l1 - l1p) * ( normEnergy - L1[n-1][0]) /
                  (L1[n][0] - L1[n-1][0]);
  return barkas;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetL2(G4double normEnergy) const
{
  G4int n;
  for(n = 0; n < sizeL2; n++) {
    if( normEnergy < L2[n][0] ) break;
  }
  if(0 == n) n = 1 ;
  if(n >= sizeL2) n = sizeL2 - 1 ;

  G4double l2    = L2[n][1];
  G4double l2p   = L2[n-1][1];
  G4double bloch = l2p + (l2 - l2p) * ( normEnergy - L2[n-1][0]) /
                  (L2[n][0] - L2[n-1][0]);
  return bloch;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4int G4QAOLowEnergyLoss::materialAvailable[6] = {13,14,29,73,79,78};
const G4int G4QAOLowEnergyLoss::nbofShellForMaterial[6] = {3,3,4,6,6,6 };

const G4double G4QAOLowEnergyLoss::alShellEnergy[3]  ={ 2795e-6, 202e-6,  16.9e-6};
const G4double G4QAOLowEnergyLoss::alShellStrength[3]={ 0.1349, 0.6387, 0.2264};
const G4double G4QAOLowEnergyLoss::siShellEnergy[3]  ={ 3179e-6, 249e-6, 20.3e-6 };
const G4double G4QAOLowEnergyLoss::siShellStrength[3]={ 0.1222, 0.5972, 0.2806};
const G4double G4QAOLowEnergyLoss::cuShellEnergy[4]  ={ 16931e-6, 1930e-6, 199e-6, 39.6e-6};
const G4double G4QAOLowEnergyLoss::cuShellStrength[4]={ 0.0505, 0.2561, 0.4913, 0.2021};
const G4double G4QAOLowEnergyLoss::taShellEnergy[6]  ={ 88926e-6, 18012e-6, 3210e-6, 575e-6, 108.7e-6, 30.8e-6};
const G4double G4QAOLowEnergyLoss::taShellStrength[6]={ 0.0126, 0.0896, 0.2599, 0.3413, 0.2057, 0.0908};
const G4double G4QAOLowEnergyLoss::auShellEnergy[6]={ 96235e-6, 25918e-6, 4116e-6, 599e-6, 87.3e-6, 36.9e-6};
const G4double G4QAOLowEnergyLoss::auShellStrength[6]={ 0.0139, 0.0803, 0.2473, 0.423, 0.1124, 0.1231};
const G4double G4QAOLowEnergyLoss::ptShellEnergy[6]={ 95017e-6, 25590e-6, 4063e-6, 576e-6, 81.9e-6, 31.4e-6};
const G4double G4QAOLowEnergyLoss::ptShellStrength[6]={ 0.0129, 0.0745, 0.2295, 0.4627, 0.1324, 0.0879};


const G4double G4QAOLowEnergyLoss::L0[67][2] =
{
  {0.00,        0.000001},
  {0.10,	0.000001},
  {0.12,	0.00001},
  {0.14,	0.00005},
  {0.16,	0.00014},
  {0.18,	0.00030},
  {0.20,	0.00057},
  {0.25,	0.00189},
  {0.30,	0.00429},
  {0.35,	0.00784},
  {0.40,	0.01248},
  {0.45,	0.01811},
  {0.50,	0.02462},
  {0.60,	0.03980},
  {0.70,	0.05731},
  {0.80,	0.07662},
  {0.90,	0.09733},
  {1.00,	0.11916},
  {1.20,	0.16532},
  {1.40,	0.21376},
  {1.60,	0.26362},
  {1.80,	0.31428},
  {2.00,	0.36532},
  {2.50,	0.49272},
  {3.00,	0.61765},
  {3.50,	0.73863},
  {4.00,	0.85496},
  {4.50,	0.96634},
  {5.00,	1.07272},
  {6.00,	1.27086},
  {7.00,	1.45075},
  {8.00,	1.61412},
  {9.00,	1.76277},
  {10.00,       1.89836},
  {12.00,       2.13625},
  {14.00,       2.33787},
  {16.00,       2.51093},
  {18.00,       2.66134},
  {20.00,       2.79358},
  {25.00, 3.06539},
  {30.00, 3.27902},
  {35.00, 3.45430},
  {40.00, 3.60281},
  {45.00, 3.73167},
  {50.00, 3.84555},
  {60.00, 4.04011},
  {70.00, 4.20264},
  {80.00, 4.34229},
  {90.00,  4.46474},
  {100.00, 4.57378},
  {120.00, 4.76155},
  {140.00, 4.91953},
  {160.00, 5.05590},
  {180.00, 5.17588},
  {200.00, 5.28299},
  {250.00, 5.50925},
  {300.00, 5.69364},
  {350.00, 5.84926},
  {400.00, 5.98388},
  {450.00, 6.10252},
  {500.00, 6.20856},
  {600.00, 6.39189},
  {700.00, 6.54677},
  {800.00, 6.68084},
  {900.00, 6.79905},
  {1000.00, 6.90474}
};

const G4double G4QAOLowEnergyLoss::L1[22][2] =
{
  {0.00,       -0.000001},
  {0.10,       -0.00001},
  {0.20,       -0.00049},
  {0.30,       -0.00084},
  {0.40,	0.00085},
  {0.50,	0.00519},
  {0.60,	0.01198},
  {0.70,	0.02074},
  {0.80,	0.03133},
  {0.90,	0.04369},
  {1.00,	0.06035},
  {2.00,	0.24023},
  {3.00,	0.44284},
  {4.00,	0.62012},
  {5.00,	0.77031},
  {6.00,	0.90390},
  {7.00,	1.02705},
  {8.00,	1.10867},
  {9.00,	1.17546},
  {10.00,       1.21599},
  {15.00,	1.24349},
  {20.00,	1.16752}
};

const G4double G4QAOLowEnergyLoss::L2[14][2] =
{
  {0.00,	0.000001},
  {0.10,	0.00001},
  {0.20,	0.00000},
  {0.40,       -0.00120},
  {0.60,       -0.00036},
  {0.80,	0.00372},
  {1.00,	0.01298},
  {2.00,	0.08296},
  {4.00,	0.21953},
  {6.00,	0.23903},
  {8.00,	0.20893},
  {10.00,	0.10879},
  {20.00,      -0.88409},
  {40.00,      -1.13902}
};

const G4int G4QAOLowEnergyLoss::nbOfElectronPerSubShell[1540] =
{
  0, // consistency with G4AtomicShells
  1,//------ H
  2,//------ He
  2,  1,//------ Li
  2,  2,//------ Be
  2,  2,  1,//------ B
  2,  2,  2,//------ C
  2,  2,  2,  1,//------ N
  2,  2,  2,  2,//------ O
  2,  2,  5,//------ F
  2,  2,  2,  4,//------ Ne
  2,  2,  2,  4,  1,//------ Na
  2,  2,  2,  4,  2,//------ Mg
  2,  2,  2,  4,  2,  1,//------ Al
  2,  2,  2,  4,  2,  2,//------ Si
  2,  2,  2,  4,  2,  3,//------ P
  2,  2,  2,  4,  2,  4,//------
  2,  2,  2,  4,  2,  5,//------
  2,  2,  2,  4,  2,  2,  4,//------
  2,  2,  2,  4,  2,  2,  4,  1,//------
  2,  2,  2,  4,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  2,  2,//------
  2,  2,  2,  4,  2,  2,  4,  3,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  5,  2,//------
  2,  2,  2,  4,  2,  2,  4,  6,  2,//------
  2,  2,  2,  4,  2,  2,  4,  7,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  5,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  1,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  3,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  4,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  5,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  1,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  2,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  3,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  5,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  6,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  7,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  5,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  1,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  3,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  4,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  5,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  2,  4,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  2,  4,  1,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  2,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  3,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  4,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  5,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  7,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  7,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  9,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6, 10,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6, 11,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6, 12,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6, 13,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  2,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  3,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  5,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  6,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  7,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  9,  1,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  1,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  1,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  3,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  4,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  2,  3,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  2,  4,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  2,  4,  1,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  2,  4,  2,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  2,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  3,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  4,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  6,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  7,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  7,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6,  8,  2,  2,  4,  1,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6, 10,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6, 11,  2,  2,  4,  2,//------
  2,  2,  2,  4,  2,  2,  4,  4,  6,  2,  2,  4,  4,  6,  6,  8,  2,  2,  4,  4,  6, 12,  2,  2,  4,  2 //-----
};

const G4int G4QAOLowEnergyLoss::fNumberOfShells[101] =
{
 0 ,  // nonexisting zero element

 1 ,  1 ,  2 ,  2 ,  3 ,    3 ,  4 ,  4 ,  3 ,  4 ,  //  1 - 10

 5 ,  5 ,  6 ,  6 ,  6 ,    6 ,  6 ,  7 ,  8 ,  8 ,  // 11 - 20

 9 ,  9 ,  9 ,  9 ,  9 ,    9 ,  9 , 10 , 10 , 10 ,  // 21 - 30

11 , 11 , 11 , 11 , 11 ,   12 , 13 , 13 , 14 , 14 ,  // 31 - 40

14 , 14 , 14 , 14 , 14 ,   15 , 15 , 15 , 16 , 16 ,  // 41 - 50

// ----------------------------------------------------------

16 , 16 , 16 , 17 , 18 ,   18 , 19 , 19 , 19 , 19 ,  // 51 - 60

19 , 19 , 19 , 20 , 19 ,   19 , 19 , 19 , 19 , 20 ,  // 61 - 70

21 , 21 , 21 , 21 , 21 ,   21 , 21 , 21 , 22 , 22 ,  // 71 - 80

23 , 23 , 23 , 23 , 24 ,   24 , 25 , 25 , 26 , 26 ,  // 81 - 90

27 , 27 , 27 , 26 , 26 ,   27 , 27 , 26 , 26 , 26    // 91 - 100

};
