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
#ifndef G4EMModelBuilder_h
#define G4EMModelBuilder_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MultipleScattering52.hh"
#include "G4eIonisation52.hh"
#include "G4eBremsstrahlung52.hh"

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
  
    G4MultipleScattering52 theElectronMultipleScattering;
    G4eIonisation52 theElectronIonisation;
    G4eBremsstrahlung52 theElectronBremsStrahlung;
  
    G4MultipleScattering52 thePositronMultipleScattering;
    G4eIonisation52 thePositronIonisation; 
    G4eBremsstrahlung52 thePositronBremsStrahlung;  
    G4eplusAnnihilation theAnnihilation;
};
// 2002 by J.P. Wellisch
#endif





