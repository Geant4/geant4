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
// *                                                                  *
// * Parts of this code which have been  developed by Abdel-Waged     *
// * et al under contract (31-465) to the King Abdul-Aziz City for    *
// * Science and Technology (KACST), the National Centre of           *
// * Mathematics and Physics (NCMP), Saudi Arabia.                    *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/include/G4UrQMD1_3Model.hh
/// \brief Definition of the G4UrQMD1_3Model class
//
#ifndef G4UrQMD1_3Model_hh
#define G4UrQMD1_3Model_hh
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:           G4UrQMD1_3Model.hh
//
// Version:          0.B
// Date:             20/10/12
// Author:           Kh. Abdel-Waged and Nuha Felemban
// Revised by:       V.V. Uzhinskii        
//                   SPONSERED BY
// Customer:         KAUST/NCMP
// Contract:         31-465
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
//
//
// Class Description - End
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4Track.hh"
#include <fstream>                
#include <string>
//--------------------------------
//AND->
//This should be included in the .cc to 
//avoid multiple definition issues
//#include "G4UrQMD1_3Interface.hh"
//AND<-


////////////////////////////////////////////////////////////////////////////////

class G4UrQMD1_3Model : public G4VIntraNuclearTransportModel {



public:

    G4UrQMD1_3Model(const G4String& name = "UrQMD1_3");

    
    virtual ~G4UrQMD1_3Model ();

//    G4double GetMinEnergy( const G4Material *aMaterial,
//                                  const G4Element *anElement ) const;
//    G4double GetMaxEnergy( const G4Material *aMaterial,
//                                  const G4Element *anElement ) const;

    
   G4ReactionProductVector* Propagate(G4KineticTrackVector* 
   theSecondaries, G4V3DNucleus* theTarget);



    virtual G4HadFinalState *ApplyYourself
       (const G4HadProjectile &, G4Nucleus &);



    
   private: 

   G4int operator==(G4UrQMD1_3Model& right) {
    return (this == &right);
   }

   G4int operator!=(G4UrQMD1_3Model& right) {
    return (this != &right);
   }

  
   G4int verbose;
                                  
   void InitialiseDataTables();
   
  G4int CurrentEvent;
    
  private:


    void WelcomeMessage () const;                         

    G4HadFinalState theResult; 
   

};

// inline G4double G4UrQMD1_3Model::GetMinEnergy( const G4Material *,
//  const G4Element * ) const
//  {return theMinEnergy;}
////////////////////////////////////////////////////////////////////////////////
//
// inline G4double G4UrQMD1_3Model::GetMaxEnergy( const G4Material *,
//  const G4Element * ) const
//  {return theMaxEnergy;}



////////////////////////////////////////////////////////////////////////////////
//
#endif
