//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestSensorHit.hh"
#include "G4VDataSetAlgorithm.hh"
#include "FluoTestDataSet.hh"
#include "G4LogLogInterpolation.hh"


G4Allocator<FluoTestSensorHit> FluoTestSensorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSensorHit::FluoTestSensorHit()
{
   EdepTot = 0.;  
   EdepDetect = 0.;
   Efficiency = 0.9;
   F = 0.15;
   epsilon = 2.96 * eV;
   deltaE = 220 * eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSensorHit::FluoTestSensorHit(const FluoTestSensorHit& right)
{
   EdepTot = right.EdepTot ; 
  EdepDetect = right.EdepDetect;
  Efficiency = right.Efficiency;
   F = right.F;
  epsilon = right.epsilon;
  deltaE = right.deltaE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestSensorHit::~FluoTestSensorHit()
{;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const FluoTestSensorHit& FluoTestSensorHit::operator=(const FluoTestSensorHit& right)
{
  EdepTot = right.EdepTot ; 
  EdepDetect = right.EdepDetect;
  Efficiency = right.Efficiency;
 F = right.F;
  epsilon = right.epsilon;
  deltaE = right.deltaE;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int FluoTestSensorHit::operator==(const FluoTestSensorHit& right) const
{
  return 0;
}


void FluoTestSensorHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestSensorHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double FluoTestSensorHit::RandomCut()
{ 

  G4String fileName = "efficienza";
  G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();
  FluoTestDataSet* dataSet = new FluoTestDataSet(1,fileName,interpolation,keV,1);

  G4double id = 0;
 
  Efficiency = dataSet->FindValue(EdepTot,id);
  
  G4double  Random = G4UniformRand(); 
  
    if ( Random<Efficiency )
      {
	G4double sigma = sqrt(F*epsilon*EdepTot+pow(deltaE/2355,2));

	RandEngine FluoTestEngine;

	EdepDetect = RandGaussQ::shoot(&FluoTestEngine, EdepTot, sigma );

  }
    else {EdepDetect = 0.;}
    return   EdepDetect;
    
};












