// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IonisParamMat.cc,v 1.4 2001-01-11 10:37:48 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// 18-07-98, bug corrected in ComputeDensityEffect() for gas
// 09-07-98, data moved from G4Material, M.Maire

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#include "G4IonisParamMat.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamMat::G4IonisParamMat(G4Material* material)
:fMaterial(material)
{
  ComputeMeanParameters();
  ComputeDensityEffect();
  ComputeFluctModel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4IonisParamMat::ComputeMeanParameters()
{
  // compute mean excitation energy and shell correction vector

  fTaul = (*(fMaterial->GetElementVector()))[0]->GetIonisation()->GetTaul();
  fLogMeanExcEnergy = 0.;

  for (G4int i=0; i < fMaterial->GetNumberOfElements(); i++)
    fLogMeanExcEnergy += (fMaterial->GetVecNbOfAtomsPerVolume())[i]
                   *((*(fMaterial->GetElementVector()))[i]->GetZ())
                   *log((*(fMaterial->GetElementVector()))[i]->GetIonisation()
                                                             ->GetMeanExcitationEnergy());

  fLogMeanExcEnergy /= fMaterial->GetTotNbOfElectPerVolume();
  fMeanExcitationEnergy = exp(fLogMeanExcEnergy);
  fShellCorrectionVector = new G4double[3];

  for (G4int j=0; j<=2; j++)
  {
    fShellCorrectionVector[j] = 0.;

    for (G4int k=0; k<fMaterial->GetNumberOfElements(); k++)
      fShellCorrectionVector[j] += (fMaterial->GetVecNbOfAtomsPerVolume())[k] 
              *((*(fMaterial->GetElementVector()))[k]->GetIonisation()
                                                     ->GetShellCorrectionVector()[j]);
     
    fShellCorrectionVector[j] /= fMaterial->GetTotNbOfElectPerVolume();
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
                    
void G4IonisParamMat::ComputeDensityEffect()
{
  // Compute parameters for the density effect correction in DE/Dx formula.
  // The parametrization is from R.M. Sternheimer, Phys. Rev.B,3:3681 (1971)

  const G4double Cd2 = 4*pi*hbarc_squared*classic_electr_radius;
  const G4double twoln10 = 2.*log(10.);

  G4int icase;
  
  fCdensity = 1. + log(fMeanExcitationEnergy*fMeanExcitationEnergy
              /(Cd2*fMaterial->GetTotNbOfElectPerVolume()));

  //
  // condensed materials
  //
  G4State State = fMaterial->GetState();
  
  if ((State == kStateSolid)||(State == kStateLiquid)) {

      const G4double E100eV  = 100.*eV; 
      const G4double ClimiS[] = {3.681 , 5.215 };
      const G4double X0valS[] = {1.0   , 1.5   };
      const G4double X1valS[] = {2.0   , 3.0   };
                                
      if(fMeanExcitationEnergy < E100eV) icase = 0 ;
         else                             icase = 1 ;

      if(fCdensity < ClimiS[icase]) fX0density = 0.2;
         else                       fX0density = 0.326*fCdensity-X0valS[icase];

      fX1density = X1valS[icase] ; fMdensity = 3.0;
      
      //special: Hydrogen
      if ((fMaterial->GetNumberOfElements()==1)&&(fMaterial->GetZ()==1.)) {
         fX0density = 0.425; fX1density = 2.0; fMdensity = 5.949;
      }
  }

  //
  // gases
  //
  if (State == kStateGas) { 

      const G4double ClimiG[] = { 10. , 10.5 , 11. , 11.5 , 12.25 , 13.804};
      const G4double X0valG[] = { 1.6 , 1.7 ,  1.8 ,  1.9 , 2.0   ,  2.0 };
      const G4double X1valG[] = { 4.0 , 4.0 ,  4.0 ,  4.0 , 4.0   ,  5.0 };

      icase = 5;
      fX0density = 0.326*fCdensity-2.5 ; fX1density = 5.0 ; fMdensity = 3. ; 
      while((icase > 0)&&(fCdensity < ClimiG[icase])) icase-- ;
      fX0density = X0valG[icase]  ; fX1density = X1valG[icase] ;
      
      //special: Hydrogen
      if ((fMaterial->GetNumberOfElements()==1)&&(fMaterial->GetZ()==1.)) {
         fX0density = 1.837; fX1density = 3.0; fMdensity = 4.754;
      }
      
      //special: Helium
      if ((fMaterial->GetNumberOfElements()==1)&&(fMaterial->GetZ()==2.)) {
         fX0density = 2.191; fX1density = 3.0; fMdensity = 3.297;
      }

      // change parameters if the gas is not in STP.
      // For the correction the density(STP) is needed. Density(STP) is calculated here : 
      
      G4double Density  = fMaterial->GetDensity();
      G4double Pressure = fMaterial->GetPressure();
      G4double Temp     = fMaterial->GetTemperature();
      
      G4double DensitySTP = Density*STP_Pressure*Temp/(Pressure*STP_Temperature);

      G4double ParCorr = log(Density/DensitySTP) ;
  
      fCdensity  -= ParCorr;
      fX0density -= ParCorr/twoln10 ;
      fX1density -= ParCorr/twoln10 ;
  }

  G4double Xa = fCdensity/twoln10 ;
  fAdensity = twoln10*(Xa-fX0density)
              /pow((fX1density-fX0density),fMdensity);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4IonisParamMat::ComputeFluctModel()
{
  // compute parameters for the energy loss fluctuation model

  // need an 'effective Z' ?????
  G4double Zeff = 0.;
  for (G4int i=0;i<fMaterial->GetNumberOfElements();i++) 
     Zeff += (fMaterial->GetFractionVector())[i]
             *((*(fMaterial->GetElementVector()))[i]->GetZ());

  if (Zeff > 2.) fF2fluct = 2./Zeff ;
    else         fF2fluct = 0.;

  fF1fluct         = 1. - fF2fluct;
  fEnergy2fluct    = 10.*Zeff*Zeff*eV;
  fLogEnergy2fluct = log(fEnergy2fluct);
  fLogEnergy1fluct = (fLogMeanExcEnergy - fF2fluct*fLogEnergy2fluct)
                     /fF1fluct;
  fEnergy1fluct    = exp(fLogEnergy1fluct);
  fEnergy0fluct    = 10.*eV;
  fRateionexcfluct = 0.4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamMat::~G4IonisParamMat()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IonisParamMat::G4IonisParamMat(const G4IonisParamMat &right)
{
    *this = right;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

const G4IonisParamMat& G4IonisParamMat::operator=(const G4IonisParamMat& right)
{
  if (this != &right)
    {
      fMaterial              = right.fMaterial;
      fMeanExcitationEnergy  = right.fMeanExcitationEnergy;
      fLogMeanExcEnergy      = right.fLogMeanExcEnergy;
      fShellCorrectionVector = right.fShellCorrectionVector;
      fTaul                  = right.fTaul;
      fCdensity              = right.fCdensity;
      fMdensity              = right.fMdensity;
      fAdensity              = right.fAdensity;
      fX0density             = right.fX0density;
      fX1density             = right.fX1density;
      fF1fluct               = right.fF1fluct;
      fF2fluct               = right.fF2fluct;
      fEnergy1fluct          = right.fEnergy1fluct;
      fLogEnergy1fluct       = right.fLogEnergy1fluct;      
      fEnergy2fluct          = right.fEnergy2fluct;
      fLogEnergy2fluct       = right.fLogEnergy2fluct;      
      fEnergy0fluct          = right.fEnergy0fluct;
      fRateionexcfluct       = right.fRateionexcfluct;
     } 
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4int G4IonisParamMat::operator==(const G4IonisParamMat &right) const
{
  return (this == (G4IonisParamMat *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4int G4IonisParamMat::operator!=(const G4IonisParamMat &right) const
{
  return (this != (G4IonisParamMat *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

