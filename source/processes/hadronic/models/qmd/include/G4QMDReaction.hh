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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name: G4QMDReaction.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 02 April 2007
// -----------------------------------------------------------------------------
//
// 081107 Add UnUseGEM (then use the default channel of G4Evaporation)
//            UseFrag (chage criterion of a inelastic reaction)
// 090331 Change member shenXS and genspaXS object to pointer 
//

#ifndef G4QMDReaction_hh
#define G4QMDReaction_hh

#include "G4QMDSystem.hh"
#include "G4QMDCollision.hh"
#include "G4QMDMeanField.hh"
#include "G4QMDParticipant.hh"

#include "G4IonsShenCrossSection.hh"
//#include "G4GeneralSpaceNNCrossSection.hh"
#include "G4PiNuclearCrossSection.hh"

#include "G4HadronicInteraction.hh"

#include "G4Evaporation.hh"
#include "G4ExcitationHandler.hh"
//#include "G4PreCompoundModel.hh"

class G4QMDReaction : public G4HadronicInteraction
{
   public:
      G4QMDReaction();
      ~G4QMDReaction();

      std::vector< G4QMDSystem* > GetFinalStates(); 

      G4HadFinalState *ApplyYourself( const G4HadProjectile &aTrack, G4Nucleus & targetNucleus );

      G4ExcitationHandler* GetExcitationHandler(){ return excitationHandler; };

      void UnUseGEM(){ gem = false; setEvaporationCh(); };
      void UseFRAG(){ frag = true; };

      void SetTMAX( G4int i ){ maxTime = i; };
      void SetDT( G4double t ){ deltaT = t; };
      void SetEF( G4double x ){ envelopF = x; };

      virtual void ModelDescription(std::ostream& outFile) const;

   private:
      //copy is unexpeced
      void operator =( const G4QMDReaction& ){};
      G4QMDReaction( const G4QMDReaction& ):G4HadronicInteraction("QMDModel"){};

      void setEvaporationCh();
      void setHighEnergyModel();

      G4QMDMeanField* meanField;

      G4QMDCollision* collision;

      void doCollision();
      std::vector< G4QMDSystem* > doClusterJudgment();
      
      G4QMDSystem* system;
      G4double deltaT;
      G4int maxTime;
      G4double envelopF;

      G4Evaporation* evaporation;
      G4ExcitationHandler* excitationHandler;

//      G4VPreCompoundModel* preco; 


      //                            b        pd_proj                pd_targ                  z_p     a_p     z_t     a_t      plab       elab
//      G4double offSetOfCollision( G4double , G4ParticleDefinition* , G4ParticleDefinition* , G4int , G4int , G4int , G4int , G4double , G4double  );  
      //                           b          pd_proj                 pd_targ                 plab       elab      bmax boostToCM
      void calcOffSetOfCollision( G4double , const G4ParticleDefinition* , const G4ParticleDefinition* , G4double , G4double , G4double , G4ThreeVector );
      G4double coulomb_collision_gamma_proj;
      G4double coulomb_collision_rx_proj;
      G4double coulomb_collision_rz_proj;
      G4double coulomb_collision_px_proj;
      G4double coulomb_collision_pz_proj;

      G4double coulomb_collision_gamma_targ;
      G4double coulomb_collision_rx_targ;
      G4double coulomb_collision_rz_targ;
      G4double coulomb_collision_px_targ;
      G4double coulomb_collision_pz_targ;

//090331
      G4IonsShenCrossSection* shenXS;
      //G4GeneralSpaceNNCrossSection* genspaXS;

      G4PiNuclearCrossSection* piNucXS;

      G4bool gem;
      G4bool frag;

};

#endif
