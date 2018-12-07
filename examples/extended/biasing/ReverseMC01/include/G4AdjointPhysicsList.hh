//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file biasing/ReverseMC01/include/G4AdjointPhysicsList.hh
/// \brief Definition of the G4AdjointPhysicsList class
//
//
//////////////////////////////////////////////////////////////
//  Class Name:             G4AdjointPhysicsList
//        Author:               L. Desorgher
//        Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//        Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////
// CHANGE HISTORY
//--------------
//      ChangeHistory:
//                 17-11-2009 creation by L. Desorgher
//
//-------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4AdjointPhysicsList_h
#define G4AdjointPhysicsList_h 1
#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "G4eIonisation.hh"
#include "G4hIonisation.hh"
class G4AdjointPhysicsMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4AdjointPhysicsList: public G4VUserPhysicsList
{
  public:
   G4AdjointPhysicsList();
   virtual ~G4AdjointPhysicsList();
   void SetLossFluctuationFlag(bool aBool);
   inline void SetUseIonisation(bool aBool){fUse_eionisation = aBool;}
   inline void SetUseProtonIonisation(bool aBool){fUse_pionisation = aBool;}
   inline void SetUseBrem(bool aBool){fUse_brem = aBool;}
   inline void SetUseCompton(bool aBool){fUse_compton = aBool;}
   inline void SetUseMS(bool aBool){fUse_ms = aBool;}
   inline void SetUsePEEffect(bool aBool){fUse_peeffect = aBool;}
   inline void SetUseGammaConversion(bool aBool){ fUse_gamma_conversion
                       = aBool;}
   inline void SetUseEgainFluctuation(bool aBool){ fUse_egain_fluctuation
                        = aBool;}
   inline void SetEminAdjModels(G4double aVal){ fEmin_adj_models = aVal;}
   inline void SetEmaxAdjModels(G4double aVal){ fEmax_adj_models = aVal;}
  
  protected:
    // Construct particle and physics
   virtual void ConstructParticle();
   virtual void ConstructProcess();
   virtual void SetCuts();

    // these methods Construct particles 
   void ConstructBosons();
   void ConstructLeptons();
   void ConstructMesons();
   void ConstructBaryons();
   void ConstructAdjointParticles();

    // these methods Construct physics processes and register them
   void ConstructGeneral();
   void ConstructEM();
   G4eIonisation* fEminusIonisation;
   G4hIonisation* fPIonisation;

  private:
   G4AdjointPhysicsMessenger* fPhysicsMessenger;
   G4bool fUse_forced_interaction;
   G4bool fUse_eionisation;
   G4bool fUse_pionisation;
   G4bool fUse_brem;
   G4bool fUse_compton;
   G4bool fUse_ms;
   G4bool fUse_egain_fluctuation;
   G4bool fUse_peeffect;
   G4bool fUse_gamma_conversion;
   G4double fEmin_adj_models;
   G4double fEmax_adj_models;
   G4double fCS_biasing_factor_compton;
   G4double fCS_biasing_factor_brem;
   G4double fCS_biasing_factor_ionisation;
   G4double fCS_biasing_factor_PEeffect;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

