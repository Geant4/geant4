#ifndef G4EMModelBuilder_h
#define G4EMModelBuilder_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MultipleScatteringSTD.hh"
#include "G4eIonisationSTD.hh"
#include "G4eBremsstrahlungSTD.hh"

class G4EMModelBuilder 
{
  public: 
    G4EMModelBuilder();
    virtual ~G4EMModelBuilder();

  public: 
    void Build();

  protected:
    G4PhotoElectricEffect thePhotoEffect;
    G4ComptonScattering theComptonEffect;
    G4GammaConversion thePairProduction;
  
    G4MultipleScatteringSTD theElectronMultipleScattering;
    G4eIonisationSTD theElectronIonisation;
    G4eBremsstrahlungSTD theElectronBremsStrahlung;
  
    G4MultipleScatteringSTD thePositronMultipleScattering;
    G4eIonisationSTD thePositronIonisation; 
    G4eBremsstrahlungSTD thePositronBremsStrahlung;  
    G4eplusAnnihilation theAnnihilation;
};
// 2002 by J.P. Wellisch
#endif





