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
// $Id: G4AdjointPhysicsList.hh,v 1.1 2009-11-19 22:41:18 ldesorgh Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointPhysicsList
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
//////////////////////////////////////////////////////////////
// CHANGE HISTORY
//--------------
//      ChangeHistory:
//	 	17-11-2009 creation by L. Desorgher
//
//-------------------------------------------------------------
#ifndef G4AdjointPhysicsList_h
#define G4AdjointPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "G4eIonisation.hh"
#include "G4hIonisation.hh"
class G4AdjointPhysicsMessenger;

class G4AdjointPhysicsList: public G4VUserPhysicsList
{
  public:
    G4AdjointPhysicsList();
   ~G4AdjointPhysicsList();
   
   
    void SetLossFluctuationFlag(bool aBool);
    inline void SetUseIonisation(bool aBool){use_eionisation = aBool;}
    inline void SetUseProtonIonisation(bool aBool){use_pionisation = aBool;}
    inline void SetUseBrem(bool aBool){use_brem = aBool;}
    inline void SetUseCompton(bool aBool){use_compton = aBool;}
    inline void SetUseMS(bool aBool){use_ms = aBool;}
    inline void SetUsePEEffect(bool aBool){use_peeffect = aBool;}
    inline void SetUseGammaConversion(bool aBool){use_gamma_conversion = aBool;}
    inline void SetUseEgainFluctuation(bool aBool){use_egain_fluctuation = aBool;}
    inline void SetEminAdjModels(G4double aVal){emin_adj_models = aVal;}
    inline void SetEmaxAdjModels(G4double aVal){emax_adj_models = aVal;}
    
  
  protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
 
    void SetCuts();

   
  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();
    void ConstructAdjointParticles();

  protected:
    // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    
    
    
//    G4eIonisationWithAdjoint* theeminusIonisation; 
   G4eIonisation* theeminusIonisation;
   G4hIonisation* thepIonisation;



  private:
  
     G4AdjointPhysicsMessenger* thePhysicsMessenger;
  
  
     G4bool use_eionisation;
     G4bool use_pionisation;
     G4bool use_brem;
     G4bool use_compton;
     G4bool use_ms; 
     G4bool use_egain_fluctuation;
     G4bool use_peeffect;
     G4bool use_gamma_conversion;
     
     G4double emin_adj_models;
     G4double emax_adj_models;
     
     G4double CS_biasing_factor_compton;
     G4double CS_biasing_factor_brem;
     G4double CS_biasing_factor_ionisation;
     G4double CS_biasing_factor_PEeffect;
};


#endif



