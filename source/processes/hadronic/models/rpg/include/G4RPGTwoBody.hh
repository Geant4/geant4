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
//
// Author: D.H. Wright
// Date:   30 May 2007
//
 
#ifndef G4RPGTwoBody_h
#define G4RPGTwoBody_h 1

// Class Description:
//
// Calculation of the hadron-nucleus interaction using simple two body
// scattering.  This method is used if the fragmentation or two cluster
// methods fail or are inappropriate for the given energy. 


#include "G4RPGReaction.hh"
#include "G4DynamicParticle.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4FastVector.hh"
#include "G4HadProjectile.hh"


 class G4RPGTwoBody : public G4RPGReaction 
 {
 public:  // with description

   G4RPGTwoBody();
    
   G4bool ReactionStage(const G4HadProjectile* /*originalIncident*/,
                        G4ReactionProduct& /*modifiedOriginal*/,
                        G4bool& /*incidentHasChanged*/,
                        const G4DynamicParticle* /*originalTarget*/,
                        G4ReactionProduct& /*targetParticle*/,
                        G4bool& /*targetHasChanged*/,
                        const G4Nucleus& /*targetNucleus*/,
                        G4ReactionProduct& /*currentParticle*/,
                        G4FastVector<G4ReactionProduct,256>& /*vec*/,
                        G4int& /*vecLen*/,
                        G4bool /*leadFlag*/,
                        G4ReactionProduct& /*leadingStrangeParticle*/);


 private:
    
 };
 
#endif
 
