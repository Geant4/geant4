// This code implementation is the intellectual property of
// the GEANT4 collaboration.
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
// 23/05/2000 MGP  Made compliant to design
//  
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

G4QAOLowEnergyLoss::G4QAOLowEnergyLoss()
{
  numberOfMaterials = 6;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4QAOLowEnergyLoss::~G4QAOLowEnergyLoss()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::EnergyLoss(const G4DynamicParticle* particle,
					const G4Material* material) const
{
  G4double kinen = particle->GetKineticEnergy();
  G4double zpart = particle->GetCharge();


  SetShellMaterial(material->GetName());
  G4double pot_ionizz = 0;
  pot_ionizz = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double dedx = 0;
  G4double v = 0;
  v =  c_light * sqrt( 2 * kinen / proton_mass_c2 );
  G4double coeff = 0;
  coeff =  (twopi * proton_mass_c2 * material-> GetTotNbOfElectPerVolume()) / ( electron_mass_c2);
  coeff* =  ( fine_structure_const * fine_structure_const * hbarc_squared ) / ( kinen );
  G4double enorm = 0;  //dueemmeviquadrosuaccatagliatoomega
  enorm = ( 2 * electron_mass_c2 * v * v  ) / ( c_squared * pot_ionizz );
  G4double v0suv = 0;
  v0suv = ( fine_structure_const * c_light) / v;
  
  G4double l = 0, l_0 = 0, l_1 = 0, l_2 = 0;
  
  for (G4int nos = 0 ; nos < nbofshell ; nos++){
    
    G4double l_0_ = 0, l_1_ = 0, l_2_ = 0;
    
    enorm = ( 2 * electron_mass_c2 * v * v  ) / ( c_squared * *(shellenergy+nos) );
    
    l_0_ = GetL0(enorm);
    l_0 + =  ( *(shellstrength+nos) ) * l_0_; 
    
    l_1_ = GetL1(enorm);
    l_1 + =  ( *(shellstrength+nos) ) * l_1_; 
    
    l_2_ = GetL2(enorm);
    l_2 + =  ( *(shellstrength+nos) ) * l_2_; 
    
  }
  
  
  l = zpart * zpart * ( l_0 + zpart * v0suv * l_1 + zpart * zpart * v0suv * v0suv * l_2);
  dedx = ( coeff * l);
  
  return dedx ; 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4QAOLowEnergyLoss::SetShellMaterial(G4String materialsymbol)
{
  if       (materialsymbol  ==  "Aluminium"){
    nbofshell = 3;
    shellenergy = al_en;
    shellstrength = al_st;         }
  else if  (materialsymbol  ==  "Silicon"  ){
    nbofshell=3;
    shellenergy = si_en;
    shellstrength = si_st;         }
  else if  (materialsymbol == "Copper")   {
    nbofshell=4;
    shellenergy = cu_en;
    shellstrength = cu_st;         }
  else if  (materialsymbol == "Tantalum"){
    nbofshell=6;
    shellenergy = ta_en;
    shellstrength = ta_st;         }
  else if  (materialsymbol == "Gold" )     {
    nbofshell=6;
    shellenergy = au_en;
    shellstrength = au_st;          }
  else if  (materialsymbol == "Platinum"){
    nbofshell=6;
    shellenergy = pt_en;
    shellstrength = pt_st;          }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetL0(G4double _enorm0) const
{
  G4double _l0 = 0, _l0p = 0;
  G4double bethe = 0;
  G4int nl0 = 0;
  do{
    nl0++;
    if ( nl0 >= 67 ) break;
    _l0 = L0[nl0][1];
    _l0p = L0[nl0-1][1];
    bethe = (_l0 - _l0p) * ( _enorm0 - L0[nl0-1][0]) / (L0[nl0][0] - L0[nl0-1][0]);
    bethe+= _l0p;
  } while( _enorm0 >= L0[nl0][0] );
  
  return bethe ;
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetL1(G4double _enorm1) const
{
  G4double _l1 = 0, _l1p = 0 ;
  G4double barkas = 0;
  G4int nl1 = 0;
  do{
    nl1++;
    if ( nl1 >= 22 ) break;
    _l1 = L1[nl1][1];
    _l1p = L1[nl1-1][1];
    barkas = (_l1 - _l1p) * ( _enorm1 - L1[nl1-1][0]) / (L1[nl1][0] - L1[nl1-1][0]);
    barkas+= _l1p;
  } while( _enorm1 >= L1[nl1][0]);
  
  return barkas;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::GetL2(G4double _enorm2) const
{
  G4double _l2 = 0, _l2p = 0;
  G4double bloch = 0;
  G4int nl2 = 0;
  do{
    nl2++;
    if ( nl2 >= 14 ) break;
    _l2 = L2[nl2][1];
    _l2p = L2[nl2-1][1];
    bloch = (_l2 - _l2p) * ( _enorm2 - L2[nl2-1][0]) / (L2[nl2][0] - L2[nl2-1][0]);
    bloch+= _l2p;
  } while( _enorm2 >= L2[nl2][0] );
  
  return bloch;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::L0[67][2] =
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
  10.00,	1.89836,
  12.00,	2.13625,
  14.00,	2.33787,
  16.00,	2.51093,
  18.00,	2.66134,
  20.00,	2.79358,
  25.00,	3.06539,
  30.00,	3.27902,
  35.00,	3.45430,
  40.00,	3.60281,
  45.00,	3.73167,
  50.00,	3.84555,
  60.00,	4.04011,
  70.00,	4.20264,
  80.00,	4.34229,
  90.00,	4.46474,
  100.00,	4.57378,
  120.00,	4.76155,
  140.00,	4.91953,
  160.00,	5.05590,
  180.00,	5.17588,
  200.00,	5.28299,
  250.00,	5.50925,
  300.00,	5.69364,
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

G4double G4QAOLowEnergyLoss::L1[22][2] =
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

G4double G4QAOLowEnergyLoss::L2[14][2] =
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//G4bool G4QAOLowEnergyLoss::IsMaterial(G4String matname){
//  for(int m=0;m<6;m++){
//    if( matname==materialavailable[m]) return true;
//  }
//  return false;
//};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String G4QAOLowEnergyLoss::materialavailable[6] = {"Aluminium",
						     "Silicon",
						     "Copper",
						     "Tantalum",
						     "Gold",
						     "Platinum"};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4QAOLowEnergyLoss::al_en[3]={ 2795e-6, 202e-6, 16.9e-6};
G4double G4QAOLowEnergyLoss::al_st[3]={ 0.1349, 0.6387, 0.2264};
G4double G4QAOLowEnergyLoss::si_en[3]={ 3179e-6,249e-6,20.3e-6 };
G4double G4QAOLowEnergyLoss::si_st[3]={ 0.1222, 0.5972,	0.2806};
G4double G4QAOLowEnergyLoss::cu_en[4]={ 16931e-6, 1930e-6, 199e-6,39.6e-6};
G4double G4QAOLowEnergyLoss::cu_st[4]={ 0.0505, 0.2561, 0.4913, 0.2021};
G4double G4QAOLowEnergyLoss::ta_en[6]={ 88926e-6, 18012e-6, 3210e-6, 575e-6, 108.7e-6, 30.8e-6};
G4double G4QAOLowEnergyLoss::ta_st[6]={ 0.0126, 0.0896, 0.2599, 0.3413, 0.2057, 0.0908};
G4double G4QAOLowEnergyLoss::au_en[6]={ 96235e-6, 25918e-6, 4116e-6, 599e-6, 87.3e-6, 36.9e-6};
G4double G4QAOLowEnergyLoss::au_st[6]={ 0.0139, 0.0803, 0.2473, 0.423, 0.1124, 0.1231};
G4double G4QAOLowEnergyLoss::pt_en[6]={ 95017e-6, 25590e-6, 4063e-6, 576e-6, 81.9e-6, 31.4e-6};
G4double G4QAOLowEnergyLoss::pt_st[6]={ 0.0129, 0.0745, 0.2295, 0.4627, 0.1324, 0.0879};


G4bool G4QAOLowEnergyLoss::IsInCharge(G4double energy, 
			    const G4ParticleDefinition* particleDefinition,
			    const G4Material* material) const
{
  G4bool isInCharge = false;

  G4bool hasMaterial = false;

  for (int m = 0; m < numberOfMaterials; m++)
    {
      G4String matName = material->GetName();
      if (matName == materialavailable[m]) hasMaterial = true;
      break;
    }
  
  if (particleDefinition == G4AntiProton::AntiProtonDefinition()
      &&
      hasMaterial
      && energy >= LowEnergyLimit() && energy <= HighEnergyLimit() )
    isInCharge = true;
  
  return isInCharge;
      
}


G4double G4QAOLowEnergyLoss::HighEnergyLimit() const
{
  G4double eMax = 2. * MeV;
  return eMax;
}
