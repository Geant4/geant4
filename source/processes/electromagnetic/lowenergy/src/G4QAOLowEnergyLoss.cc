// This code implementation is the intellectual property of
// the GEANT4 collaboration
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: New Implementation
//    
//      ---------- G4QAOLowEnergyLoss physics process -------
//                  by Stephane Chauvie, 5 May 2000 
// Modified:
// 24/05/2000 MGP  Modified to remove compilation warnings on Linux and DEC
//                 Introduced sizes of L0, L1, L2 arrays
// 23/05/2000 MGP  Made compliant to design
//  
// 02/08/2000 V.Ivanchenko Clean up according new design
//
// ************************************************************
// It is the Quantal Harmonic Oscillator Model for energy loss
// of slow antiproton 
// ************************************************************
// --------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4QAOLowEnergyLoss.hh"
#include "PhysicalConstants.h"
#include "SystemOfUnits.h"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4AntiProton.hh"

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

G4double G4QAOLowEnergyLoss::HighEnergyLimit(
                             const G4ParticleDefinition* aParticle,
                             const G4Material* material) const
{
  return 2.0*MeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::LowEnergyLimit(
                             const G4ParticleDefinition* aParticle,
                             const G4Material* material) const
{
  return 50.0*keV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::HighEnergyLimit(
                             const G4ParticleDefinition* aParticle) const
{
  return 2.0*MeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::LowEnergyLimit(
                             const G4ParticleDefinition* aParticle) const
{
  return 50.0*keV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4QAOLowEnergyLoss::IsInCharge(
			    const G4DynamicParticle* particle,
			    const G4Material* material) const
{
  G4bool isInCharge = false;

  G4bool hasMaterial = false;

  for (G4int m = 0; m < numberOfMaterials; m++)
    {
      G4String matName = material->GetName();
      if (matName == materialAvailable[m]){ 
	hasMaterial = true;
	break;}
    }
  
  if ((particle->GetDefinition()) == (G4AntiProton::AntiProtonDefinition())
               && hasMaterial) isInCharge = true;
  
  return isInCharge;
      
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4QAOLowEnergyLoss::IsInCharge(
			    const G4ParticleDefinition* aParticle,
			    const G4Material* material) const
{
  G4bool isInCharge = false;

  G4bool hasMaterial = false;

  for (G4int m = 0; m < numberOfMaterials; m++)
    {
      G4String matName = material->GetName();
      if (matName == materialAvailable[m]){ 
	hasMaterial = true;
	break;}
    }
  
  if (aParticle == (G4AntiProton::AntiProtonDefinition())
                && hasMaterial) isInCharge = true;
  
  return isInCharge;
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::TheValue(const G4DynamicParticle* particle,
	       	                      const G4Material* material) 
{
  G4int zParticle = G4int(particle->GetCharge());

  G4double energy = particle->GetKineticEnergy() ;
  G4double eloss  = EnergyLoss(material,energy,zParticle) ;

  return eloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::TheValue(const G4ParticleDefinition* aParticle,
       		                      const G4Material* material,
                                            G4double kineticEnergy) 
{
  G4int zParticle = G4int (aParticle->GetPDGCharge());

  G4double eloss  = EnergyLoss(material,kineticEnergy,zParticle) ;

  return eloss ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::EnergyLoss(const G4Material* material,
                                              G4double kineticEnergy,
                                              G4int zParticle) const 
{
  G4int nbOfShell = GetNumberOfShell(material);
  G4double ionisationEnergy = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double dedx=0;
  G4double v=0;
  v= c_light * sqrt( 2 * kineticEnergy / proton_mass_c2 );
  G4double coeff=0;
  coeff= (twopi * proton_mass_c2 * material-> GetTotNbOfElectPerVolume()) / ( electron_mass_c2);
  coeff*= ( fine_structure_const * fine_structure_const * hbarc_squared ) / ( kineticEnergy );
  G4double fractionOfBetheVelocity = 0;
  fractionOfBetheVelocity = ( fine_structure_const * c_light) / v;
  
  G4double stoppingNumber = 0, l0Term = 0, l1Term = 0, l2Term = 0;
  
  for (G4int nos = 0 ; nos < nbOfShell ; nos++){
    
    G4double l0 = 0, l1 = 0, l2 = 0;
    G4double NormalizedEnergy = 0;
    NormalizedEnergy = ( 2 * electron_mass_c2 * v * v  ) / ( c_squared * GetShellEnergy(material,nos) );
    l0 = GetL0(NormalizedEnergy);
    l0Term += GetShellStrength(material,nos)  * l0; 
    
    l1 = GetL1(NormalizedEnergy);
    l1Term += GetShellStrength(material,nos) * l1; 
    
    l2 = GetL2(NormalizedEnergy);
    l2Term += GetShellStrength(material,nos) * l2; 
    
  }

       
  stoppingNumber = zParticle * zParticle * ( l0Term + zParticle * fractionOfBetheVelocity * l1Term + zParticle * zParticle * fractionOfBetheVelocity * fractionOfBetheVelocity * l2Term);
  dedx = ( coeff * stoppingNumber);
              
  return dedx ; 
                            
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4QAOLowEnergyLoss::GetNumberOfShell(const G4Material* material) const
{
  // Set default return value
  G4int nShell = nbofShellForMaterial[0];

  if(material->GetName() == "Aluminium")  nShell =  nbofShellForMaterial[0];
  else if  (material->GetName() == "Silicon"  ) nShell = nbofShellForMaterial[1] ;
  else if  (material->GetName()== "Copper") nShell = nbofShellForMaterial[2];  
  else if  (material->GetName() == "Tantalum") nShell = nbofShellForMaterial[3];
  else if  (material->GetName() == "Gold" )  nShell = nbofShellForMaterial[4];  
  else if  (material->GetName() == "Platinum") nShell = nbofShellForMaterial[5];
  else G4cout << "WARNING - G4QAOLowEnergyLoss::GetNumberOfShell - "
	      << "The model is not available for "
	      << material->GetName() 
	      << G4endl;
  
  return nShell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetShellEnergy(const G4Material* material,G4int nbOfTheShell) const
{
  // 
  G4double shellEnergy = alShellEnergy[0];

  if(material->GetName() == "Aluminium") shellEnergy =  alShellEnergy[nbOfTheShell];
  else if  (material->GetName() == "Silicon"  ) shellEnergy =  siShellEnergy[nbOfTheShell];
  else if  (material->GetName() == "Copper") shellEnergy =  cuShellEnergy[nbOfTheShell];  
  else if  (material->GetName() == "Tantalum") shellEnergy =  taShellEnergy[nbOfTheShell];
  else if  (material->GetName() == "Gold" )  shellEnergy =  auShellEnergy[nbOfTheShell];   
  else if  (material->GetName() == "Platinum") shellEnergy =  ptShellEnergy[nbOfTheShell];
  else G4cout << "WARNING - G4QAOLowEnergyLoss::GetShellEnergy - "
	      << "The model is not available for "
	      << material->GetName() 
	      << G4endl;

  return  shellEnergy;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetShellStrength(const G4Material* material,G4int nbOfTheShell) const
{
  G4double shellStrength = alShellStrength[0];
  
  if(material->GetName() == "Aluminium") shellStrength = alShellStrength[nbOfTheShell];
  else if  (material->GetName() == "Silicon"  ) shellStrength = siShellStrength[nbOfTheShell];
  else if  (material->GetName() == "Copper") shellStrength = cuShellStrength[nbOfTheShell];  
  else if  (material->GetName() == "Tantalum") shellStrength = taShellStrength[nbOfTheShell];
  else if  (material->GetName() == "Gold" )  shellStrength = auShellStrength[nbOfTheShell];   
  else if  (material->GetName() == "Platinum") shellStrength = ptShellStrength[nbOfTheShell];
  else G4cout << "WARNING - G4QAOLowEnergyLoss::GetShellEnergy - "
	      << "The model is not available for "
	      << material->GetName() 
	      << G4endl;

  return shellStrength;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetL0(G4double normEnergy) const 
{
  G4double l0 = 0, l0p = 0;
  G4double bethe = 0;
  G4int n = 0;
  do{
    n++;
    if ( n >= sizeL0 ) break;
    l0 = L0[n][1];
    l0p = L0[n-1][1];
    bethe = (l0 - l0p) * ( normEnergy - L0[n-1][0]) / (L0[n][0] - L0[n-1][0]);
    bethe+= l0p;
  } while( normEnergy >= L0[n][0] );
  
  return bethe ;
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetL1(G4double normEnergy) const
{
  G4double l1 = 0, l1p = 0 ;
  G4double barkas = 0;
  G4int n = 0;
  do{
    n++;
    if ( n >= sizeL1 ) break;
    l1 = L1[n][1];
    l1p = L1[n-1][1];
    barkas = (l1 - l1p) * ( normEnergy - L1[n-1][0]) / (L1[n][0] - L1[n-1][0]);
    barkas+= l1p;
  } while( normEnergy >= L1[n][0]);
  
  return barkas;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetL2(G4double normEnergy) const
{
  G4double l2 = 0, l2p = 0;
  G4double bloch = 0;
  G4int n = 0;
  do{
    n++;
    if ( n >= sizeL2 ) break;
    l2 = L2[n][1];
    l2p = L2[n-1][1];
    bloch = (l2 - l2p) * ( normEnergy - L2[n-1][0]) / (L2[n][0] - L2[n-1][0]);
    bloch+= l2p;
  } while( normEnergy >= L2[n][0] );
  
  return bloch;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4String G4QAOLowEnergyLoss::materialAvailable[6] = {
  "Aluminium",
  "Silicon",
  "Copper",
  "Tantalum",
  "Gold",
  "Platinum"};

const G4int G4QAOLowEnergyLoss::nbofShellForMaterial[6] = {3,3,4,6,6,6 };
              
G4double G4QAOLowEnergyLoss::alShellEnergy[3]  ={  2795e-6,  202e-6,  16.9e-6};
G4double G4QAOLowEnergyLoss::alShellStrength[3]={ 0.1349, 0.6387, 0.2264};
G4double G4QAOLowEnergyLoss::siShellEnergy[3]  ={  3179e-6, 249e-6, 20.3e-6 };
G4double G4QAOLowEnergyLoss::siShellStrength[3]={ 0.1222, 0.5972, 0.2806};
G4double G4QAOLowEnergyLoss::cuShellEnergy[4]  ={ 16931e-6, 1930e-6, 199e-6, 39.6e-6};
G4double G4QAOLowEnergyLoss::cuShellStrength[4]={ 0.0505, 0.2561, 0.4913, 0.2021};
G4double G4QAOLowEnergyLoss::taShellEnergy[6]  ={ 88926e-6, 18012e-6, 3210e-6, 575e-6, 108.7e-6, 30.8e-6};
G4double G4QAOLowEnergyLoss::taShellStrength[6]={ 0.0126, 0.0896, 0.2599, 0.3413, 0.2057, 0.0908};
G4double G4QAOLowEnergyLoss::auShellEnergy[6]={ 96235e-6, 25918e-6, 4116e-6, 599e-6, 87.3e-6, 36.9e-6};
G4double G4QAOLowEnergyLoss::auShellStrength[6]={ 0.0139, 0.0803, 0.2473, 0.423, 0.1124, 0.1231};
G4double G4QAOLowEnergyLoss::ptShellEnergy[6]={ 95017e-6, 25590e-6, 4063e-6, 576e-6, 81.9e-6, 31.4e-6};
G4double G4QAOLowEnergyLoss::ptShellStrength[6]={ 0.0129, 0.0745, 0.2295, 0.4627, 0.1324, 0.0879};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4double G4QAOLowEnergyLoss::L0[67][2] =
{
  0.00,  0.000001,
  0.10,	0.000001,
  0.12,	0.00001,
  0.14,	0.00005,
  0.16,	0.00014,
  0.18,	0.00030,
  0.20,	0.00057,
  0.25,	0.00189,
  0.30,	0.00429,
  0.35,	0.00784,
  0.40,	0.01248,
  0.45,	0.01811,
  0.50,	0.02462,
  0.60,	0.03980,
  0.70,	0.05731,
  0.80,	0.07662,
  0.90,	0.09733,
  1.00,	0.11916,
  1.20,	0.16532,
  1.40,	0.21376,
  1.60,	0.26362,
  1.80,	0.31428,
  2.00,	0.36532,
  2.50,	0.49272,
  3.00,	0.61765,
  3.50,	0.73863,
  4.00,	0.85496,
  4.50,	0.96634,
  5.00,	1.07272,
  6.00,	1.27086,
  7.00,	1.45075,
  8.00,	1.61412,
  9.00,	1.76277,
  10.00,  1.89836,
  12.00,  2.13625,
  14.00,  2.33787,
  16.00,  2.51093,
  18.00,  2.66134,
  20.00,  2.79358,
  25.00,  3.06539,
  30.00,  3.27902,
  35.00,  3.45430,
  40.00,  3.60281,
  45.00,  3.73167,
  50.00,  3.84555,
  60.00,  4.04011,
  70.00,  4.20264,
  80.00,  4.34229,
  90.00,   4.46474,
  100.00,  4.57378,
  120.00,  4.76155,
  140.00,  4.91953,
  160.00,  5.05590,
  180.00,  5.17588,
  200.00,  5.28299,
  250.00,  5.50925,
  300.00,  5.69364,
  350.00,  5.84926,
  400.00,  5.98388,
  450.00,  6.10252,
  500.00,  6.20856,
  600.00,  6.39189,
  700.00,  6.54677,
  800.00,  6.68084,
  900.00,  6.79905,
  1000.00,  6.90474
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4double G4QAOLowEnergyLoss::L1[22][2] =
{
  0.00,  -0.000001,
  0.10,  -0.00001,
  0.20,  -0.00049,
  0.30,  -0.00084,
  0.40,	0.00085,
  0.50,	0.00519,
  0.60,	0.01198,
  0.70,	0.02074,
  0.80,	0.03133,
  0.90,	0.04369,
  1.00,	0.06035,
  2.00,	0.24023,
  3.00,	0.44284,
  4.00,	0.62012,
  5.00,	0.77031,
  6.00,	0.90390,
  7.00,	1.02705,
  8.00,	1.10867,
  9.00,	1.17546,
  10.00,   1.21599,
  15.00,	1.24349,
  20.00,	1.16752
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4double G4QAOLowEnergyLoss::L2[14][2] =
{
  0.00,	0.000001,
  0.10,	0.00001,
  0.20,	0.00000,
  0.40,  -0.00120,
  0.60,  -0.00036,
  0.80,	0.00372,
  1.00,	0.01298,
  2.00,	0.08296,
  4.00,	0.21953,
  6.00,	0.23903,
  8.00,	0.20893,
  10.00,	0.10879,
  20.00,  -0.88409,	
  40.00,  -1.13902
};

