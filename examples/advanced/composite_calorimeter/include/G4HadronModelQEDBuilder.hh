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
#ifndef G4HadronModelQEDBuilder_h
#define G4HadronModelQEDBuilder_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4MultipleScattering52.hh"
#include "G4hIonisation52.hh"
#include "G4ProcessManager.hh"

class G4HadronModelQEDBuilder 
{
  public: 
    G4HadronModelQEDBuilder();
    virtual ~G4HadronModelQEDBuilder();

  public: 
    virtual void Build();

  private:  
    void RegisterOne(G4ProcessManager* aP, G4MultipleScattering52 * aM, G4hIonisation52* aI);

  private:
  
   // Pi + 
   G4MultipleScattering52 thePionPlusMult;
   G4hIonisation52 thePionPlusIonisation;

   // Pi -
   G4MultipleScattering52 thePionMinusMult;
   G4hIonisation52 thePionMinusIonisation;

   // K + 
   G4MultipleScattering52 theKaonPlusMult;
   G4hIonisation52 theKaonPlusIonisation;

   // K -
   G4MultipleScattering52 theKaonMinusMult;
   G4hIonisation52 theKaonMinusIonisation;

   // Proton
   G4MultipleScattering52 theProtonMult;
   G4hIonisation52 theProtonIonisation;
 
   // anti-proton
   G4MultipleScattering52 theAntiProtonMult;
   G4hIonisation52 theAntiProtonIonisation; 
      
   // SigmaMinus
   G4MultipleScattering52 theSigmaMinusMult;
   G4hIonisation52 theSigmaMinusIonisation;
  
   // AntiSigmaMinus
   G4MultipleScattering52 theAntiSigmaMinusMult;
   G4hIonisation52 theAntiSigmaMinusIonisation;
   
   // SigmaPlus
   G4MultipleScattering52 theSigmaPlusMult;
   G4hIonisation52 theSigmaPlusIonisation;
  
   // AntiSigmaPlus
   G4MultipleScattering52 theAntiSigmaPlusMult;
   G4hIonisation52 theAntiSigmaPlusIonisation;
    
   // XiMinus
   G4MultipleScattering52 theXiMinusMult;
   G4hIonisation52 theXiMinusIonisation;

   // AntiXiMinus
   G4MultipleScattering52 theAntiXiMinusMult;
   G4hIonisation52 theAntiXiMinusIonisation;
  
   // OmegaMinus
   G4MultipleScattering52 theOmegaMinusMult;
   G4hIonisation52 theOmegaMinusIonisation;
   
   // AntiOmegaMinus
   G4MultipleScattering52 theAntiOmegaMinusMult;
   G4hIonisation52  theAntiOmegaMinusIonisation;

   
};

// 2002 by J.P. Wellisch

#endif

