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

#include "G4MultipleScatteringSTD.hh"
#include "G4hIonisationSTD.hh"
#include "G4ProcessManager.hh"

class G4HadronModelQEDBuilder 
{
  public: 
    G4HadronModelQEDBuilder();
    virtual ~G4HadronModelQEDBuilder();

  public: 
    virtual void Build();

  private:  
    void RegisterOne(G4ProcessManager* aP, G4MultipleScatteringSTD * aM, G4hIonisationSTD* aI);

  private:
  
   // Pi + 
   G4MultipleScatteringSTD thePionPlusMult;
   G4hIonisationSTD thePionPlusIonisation;

   // Pi -
   G4MultipleScatteringSTD thePionMinusMult;
   G4hIonisationSTD thePionMinusIonisation;

   // K + 
   G4MultipleScatteringSTD theKaonPlusMult;
   G4hIonisationSTD theKaonPlusIonisation;

   // K -
   G4MultipleScatteringSTD theKaonMinusMult;
   G4hIonisationSTD theKaonMinusIonisation;

   // Proton
   G4MultipleScatteringSTD theProtonMult;
   G4hIonisationSTD theProtonIonisation;
 
   // anti-proton
   G4MultipleScatteringSTD theAntiProtonMult;
   G4hIonisationSTD theAntiProtonIonisation; 
      
   // SigmaMinus
   G4MultipleScatteringSTD theSigmaMinusMult;
   G4hIonisationSTD theSigmaMinusIonisation;
  
   // AntiSigmaMinus
   G4MultipleScatteringSTD theAntiSigmaMinusMult;
   G4hIonisationSTD theAntiSigmaMinusIonisation;
   
   // SigmaPlus
   G4MultipleScatteringSTD theSigmaPlusMult;
   G4hIonisationSTD theSigmaPlusIonisation;
  
   // AntiSigmaPlus
   G4MultipleScatteringSTD theAntiSigmaPlusMult;
   G4hIonisationSTD theAntiSigmaPlusIonisation;
    
   // XiMinus
   G4MultipleScatteringSTD theXiMinusMult;
   G4hIonisationSTD theXiMinusIonisation;

   // AntiXiMinus
   G4MultipleScatteringSTD theAntiXiMinusMult;
   G4hIonisationSTD theAntiXiMinusIonisation;
  
   // OmegaMinus
   G4MultipleScatteringSTD theOmegaMinusMult;
   G4hIonisationSTD theOmegaMinusIonisation;
   
   // AntiOmegaMinus
   G4MultipleScatteringSTD theAntiOmegaMinusMult;
   G4hIonisationSTD theAntiOmegaMinusIonisation;

   
};

// 2002 by J.P. Wellisch

#endif

