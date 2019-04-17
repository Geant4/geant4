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
#ifndef G4HIJING_Model_hh
#define G4HIJING_Model_hh
//
// MODULE:           G4HIJING_Model.hh
//
// Version:          1.B
// Date:      10/90/2013
// Author:    Khaled Abdel-Waged 
// Institute: Umm Al-Qura University
//Country:   SAUDI ARABIA
//
// Class Description
//
//
// Class Description - End
//
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
//---------------

class G4HadFinalState;   //khaled new
class HistoManager;
//-------------------------------------------------------------

class G4HIJING_Model : public G4VIntraNuclearTransportModel {



public:

    G4HIJING_Model(const G4String& name = "HIJING");

    
    virtual ~G4HIJING_Model ();

//    G4double GetMinEnergy( const G4Material *aMaterial,
//                             const G4Element *anElement ) const;

//    G4double GetMaxEnergy( const G4Material *aMaterial,
//                           const G4Element *anElement ) const;

    
   G4ReactionProductVector* Propagate(G4KineticTrackVector* 
   theSecondaries, G4V3DNucleus* theTarget);



    virtual G4HadFinalState *ApplyYourself
       (const G4HadProjectile &, G4Nucleus &);



    
   private: 

   G4bool operator==(G4HIJING_Model& right) {
    return (this == &right);
   }

   G4bool operator!=(G4HIJING_Model& right) {
    return (this != &right);
   }
   

// -------------initialize HIJING ------------------------
float efrm;        //CM energy per nucleon (GeV/A)
float bmin;        //impact parameter
float bmax;
G4int Nproduce;    //Number of produced particles
//-------------------------------------------------------
  
   G4int verbose;   //print options

// 
                                  
  void InitialiseDataTables();
  G4double Eplab(G4double, G4double);
   
  G4int CurrentEvent;

//  void ConvertToFortran(char*, const char*); 

  private:


    void WelcomeMessage () const;                         
// -------------------------------------------
    G4HadFinalState theResult;   //khaled new
    HistoManager *fHistoManager;
//-------------------------------------------   

};

#endif
