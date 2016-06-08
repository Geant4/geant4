// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// For information related to this code contact:
// Geant4 Collaboration
//
// File name:     G4IonChuFluctuationModel
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 18 August 2000
//
// Modifications: 
// 18/08/2000  V.Ivanchenko First implementation
// 04/09/2000  V.Ivanchenko Rename fluctuations            
// 03/10/2000  V.Ivanchenko CodeWizard clean up
//
// -------------------------------------------------------------------
//
// Class Description: 
//
// The aproximation of additional ion energy loss fluctuations 
// W.K.Chu, In: Ion Beam Handbook for Material Analysis.
// eds. J.W. Mayer and E. Rimini (Academic Press, New York, 1977).
// Q.Yang et al., NIM B61(1991)149-155.
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4IonChuFluctuationModel.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4IonChuFluctuationModel::G4IonChuFluctuationModel(const G4String& name)
  : G4VLowEnergyModel(name) 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4IonChuFluctuationModel::~G4IonChuFluctuationModel() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonChuFluctuationModel::TheValue(const G4DynamicParticle* particle,
                	                    const G4Material* material) 
{
  G4double energy = particle->GetKineticEnergy() ;
  G4double particleMass = particle->GetMass() ;

  G4double q = ChuFluctuationModel(material,energy,particleMass) ;

  return q ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonChuFluctuationModel::TheValue(
                                   const G4ParticleDefinition* aParticle,
       		                   const G4Material* material,
                                         G4double kineticEnergy) 
{
  G4double particleMass = aParticle->GetPDGMass() ;

  G4double q = ChuFluctuationModel(material,kineticEnergy,particleMass);

  return q ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonChuFluctuationModel::HighEnergyLimit(
                             const G4ParticleDefinition* aParticle,
                             const G4Material* material) const
{
  return 1.0*TeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonChuFluctuationModel::LowEnergyLimit(
                             const G4ParticleDefinition* aParticle,
                             const G4Material* material) const
{
  return 0.0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonChuFluctuationModel::HighEnergyLimit(
                             const G4ParticleDefinition* aParticle) const
{
  return 1.0*TeV ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonChuFluctuationModel::LowEnergyLimit(
                             const G4ParticleDefinition* aParticle) const
{
  return 0.0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4bool G4IonChuFluctuationModel::IsInCharge(
                           const G4DynamicParticle* particle,
    	                   const G4Material* material) const
{
  return true ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4bool G4IonChuFluctuationModel::IsInCharge(
                           const G4ParticleDefinition* aParticle,
      	                   const G4Material* material) const
{
  return true ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonChuFluctuationModel::ChuFluctuationModel(
                                   const G4Material* material,
                                         G4double kineticEnergy,
                                         G4double particleMass) const
{
  // The aproximation of energy loss fluctuations 
  // Q.Yang et al., NIM B61(1991)149-155.

  // Reduced energy in MeV/AMU
  G4double energy = kineticEnergy * amu_c2/(particleMass*MeV) ;

  G4double zeff = (material->GetElectronDensity())/
                  (material->GetTotNbOfAtomsPerVolume()) ;

  static G4double a[96][4] = {
 -0.3291, -0.8312,  0.2460, -1.0220,
 -0.5615, -0.5898,  0.5205, -0.7258,
 -0.5280, -0.4981,  0.5519, -0.5865,
 -0.5125, -0.4625,  0.5660, -0.5190,
 -0.5127, -0.8595,  0.5626, -0.8721,
 -0.5174, -1.1930,  0.5565, -1.1980,
 -0.5179, -1.1850,  0.5560, -1.2070,
 -0.5209, -0.9355,  0.5590, -1.0250,
 -0.5255, -0.7766,  0.5720, -0.9412,

 -0.5776, -0.6665,  0.6598, -0.8484,
 -0.6013, -0.6045,  0.7321, -0.7671,
 -0.5781, -0.5518,  0.7605, -0.6919,
 -0.5587, -0.4981,  0.7835, -0.6195,
 -0.5466, -0.4656,  0.7978, -0.5771,
 -0.5406, -0.4690,  0.8031, -0.5718,
 -0.5391, -0.5061,  0.8024, -0.5974,
 -0.5380, -0.6483,  0.7962, -0.6970,
 -0.5355, -0.7722,  0.7962, -0.7839,
 -0.5329, -0.7720,  0.7988, -0.7846,

 -0.5335, -0.7671,  0.7984, -0.7933,
 -0.5324, -0.7612,  0.7998, -0.8031,
 -0.5305, -0.7300,  0.8031, -0.7990,
 -0.5307, -0.7178,  0.8049, -0.8216,
 -0.5248, -0.6621,  0.8165, -0.7919,
 -0.5180, -0.6502,  0.8266, -0.7986,
 -0.5084, -0.6408,  0.8396, -0.8048,
 -0.4967, -0.6331,  0.8549, -0.8093,
 -0.4861, -0.6508,  0.8712, -0.8432,
 -0.4700, -0.6186,  0.8961, -0.8132,

 -0.4545, -0.5720,  0.9227, -0.7710,
 -0.4404, -0.5226,  0.9481, -0.7254,
 -0.4288, -0.4778,  0.9701, -0.6850,
 -0.4199, -0.4425,  0.9874, -0.6539,
 -0.4131, -0.4188,  0.9998, -0.6332,
 -0.4089, -0.4057,  1.0070, -0.6218,
 -0.4039, -0.3913,  1.0150, -0.6107,
 -0.3987, -0.3698,  1.0240, -0.5938,
 -0.3977, -0.3608,  1.0260, -0.5852,
 -0.3972, -0.3600,  1.0260, -0.5842,

 -0.3985, -0.3803,  1.0200, -0.6013,
 -0.3985, -0.3979,  1.0150, -0.6168,
 -0.3968, -0.3990,  1.0160, -0.6195,
 -0.3971, -0.4432,  1.0050, -0.6591,
 -0.3944, -0.4665,  1.0010, -0.6825,
 -0.3924, -0.5109,  0.9921, -0.7235,
 -0.3882, -0.5158,  0.9947, -0.7343,
 -0.3838, -0.5125,  0.9999, -0.7370,
 -0.3786, -0.4976,  1.0090, -0.7310,
 -0.3741, -0.4738,  1.0200, -0.7155,

 -0.3969, -0.4496,  1.0320, -0.6982,
 -0.3663, -0.4297,  1.0430, -0.6828,
 -0.3630, -0.4120,  1.0530, -0.6689,
 -0.3597, -0.3964,  1.0620, -0.6564,
 -0.3555, -0.3809,  1.0720, -0.6454,
 -0.3525, -0.3607,  1.0820, -0.6289,
 -0.3505, -0.3465,  1.0900, -0.6171,
 -0.3397, -0.3570,  1.1020, -0.6384,
 -0.3314, -0.3552,  1.1130, -0.6441,
 -0.3235, -0.3531,  1.1230, -0.6498,

 -0.3150, -0.3483,  1.1360, -0.6539,
 -0.3060, -0.3441,  1.1490, -0.6593,
 -0.2968, -0.3396,  1.1630, -0.6649,
 -0.2935, -0.3225,  1.1760, -0.6527,
 -0.2797, -0.3262,  1.1940, -0.6722,
 -0.2704, -0.3202,  1.2100, -0.6770,
 -0.2815, -0.3227,  1.2480, -0.6775,
 -0.2880, -0.3245,  1.2810, -0.6801,
 -0.3034, -0.3263,  1.3270, -0.6778,
 -0.2936, -0.3215,  1.3430, -0.6835,

 -0.3282, -0.3200,  1.3980, -0.6650,
 -0.3260, -0.3070,  1.4090, -0.6552,
 -0.3511, -0.3074,  1.4470, -0.6442,
 -0.3501, -0.3064,  1.4500, -0.6442,
 -0.3490, -0.3027,  1.4550, -0.6418,
 -0.3487, -0.3048,  1.4570, -0.6447,
 -0.3478, -0.3074,  1.4600, -0.6483,
 -0.3501, -0.3283,  1.4540, -0.6669,
 -0.3494, -0.3373,  1.4550, -0.6765,
 -0.3485, -0.3373,  1.4570, -0.6774,

 -0.3462, -0.3300,  1.4630, -0.6728,
 -0.3462, -0.3225,  1.4690, -0.6662,
 -0.3453, -0.3094,  1.4790, -0.6553,
 -0.3844, -0.3134,  1.5240, -0.6412,
 -0.3848, -0.3018,  1.5310, -0.6303,
 -0.3862, -0.2955,  1.5360, -0.6237,
 -0.4262, -0.2991,  1.5860, -0.6115,
 -0.4278, -0.2910,  1.5900, -0.6029,
 -0.4303, -0.2817,  1.5940, -0.5927,
 -0.4315, -0.2719,  1.6010, -0.5829,

 -0.4359, -0.2914,  1.6050, -0.6010,
 -0.4365, -0.2982,  1.6080, -0.6080,
 -0.4253, -0.3037,  1.6120, -0.6150,
 -0.4335, -0.3245,  1.6160, -0.6377,
 -0.4307, -0.3292,  1.6210, -0.6447,
 -0.4284, -0.3204,  1.6290, -0.6380,
 -0.4227, -0.3217,  1.6360, -0.6438,
  } ;

  G4int iz = (G4int)zeff - 2 ;
  if( 0 > iz ) iz = 0 ;
  if(95 < iz ) iz = 95 ;

  G4double q = 1.0 / (1.0 + a[iz][0]*pow(energy,a[iz][1])+
                          + a[iz][2]*pow(energy,a[iz][3])) ; 

  return q ;    
}


