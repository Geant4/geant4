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
#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"

#include "G4MultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ProcessManager.hh"

class G4HadronQEDBuilder 
{
  public: 
    G4HadronQEDBuilder();
    virtual ~G4HadronQEDBuilder();

  public: 
    virtual void Build();

  private:  
    void RegisterOne(G4ProcessManager* aP, G4MultipleScattering * aM, G4hIonisation* aI);

  private:
  
   // Pi + 
   G4MultipleScattering thePionPlusMult;
   G4hIonisation thePionPlusIonisation;

   // Pi -
   G4MultipleScattering thePionMinusMult;
   G4hIonisation thePionMinusIonisation;

   // K + 
   G4MultipleScattering theKaonPlusMult;
   G4hIonisation theKaonPlusIonisation;

   // K -
   G4MultipleScattering theKaonMinusMult;
   G4hIonisation theKaonMinusIonisation;

   // Proton
   G4MultipleScattering theProtonMult;
   G4hIonisation theProtonIonisation;
 
   // anti-proton
   G4MultipleScattering theAntiProtonMult;
   G4hIonisation theAntiProtonIonisation; 
      
   // SigmaMinus
   G4MultipleScattering theSigmaMinusMult;
   G4hIonisation theSigmaMinusIonisation;
  
   // AntiSigmaMinus
   G4MultipleScattering theAntiSigmaMinusMult;
   G4hIonisation theAntiSigmaMinusIonisation;
   
   // SigmaPlus
   G4MultipleScattering theSigmaPlusMult;
   G4hIonisation theSigmaPlusIonisation;
  
   // AntiSigmaPlus
   G4MultipleScattering theAntiSigmaPlusMult;
   G4hIonisation theAntiSigmaPlusIonisation;
    
   // XiMinus
   G4MultipleScattering theXiMinusMult;
   G4hIonisation theXiMinusIonisation;

   // AntiXiMinus
   G4MultipleScattering theAntiXiMinusMult;
   G4hIonisation theAntiXiMinusIonisation;
  
   // OmegaMinus
   G4MultipleScattering theOmegaMinusMult;
   G4hIonisation theOmegaMinusIonisation;
   
   // AntiOmegaMinus
   G4MultipleScattering theAntiOmegaMinusMult;
   G4hIonisation theAntiOmegaMinusIonisation;

   
};

// 2002 by J.P. Wellisch

#endif

