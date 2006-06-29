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
#ifndef G4ParametrizedHadronicVertex_h
#define G4ParametrizedHadronicVertex_h 1

#include "globals.hh"
#include "G4ParticleChange.hh"
#include "G4Nucleus.hh"
#include "G4ReactionProduct.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4Track.hh"

class G4ParametrizedHadronicVertex
{
  public:
   G4ParametrizedHadronicVertex()
   {
     theLowEPionPlus = new G4LEPionPlusInelastic;
     theLowEPionMinus = new G4LEPionMinusInelastic;
     theHighEPionPlus = new G4HEPionPlusInelastic;
     theHighEPionMinus = new G4HEPionMinusInelastic;
   }
   G4VParticleChange * ApplyYourself(G4Nucleus & theTarget, 
                                     const G4Track &thePhoton);

  private:
   G4LEPionPlusInelastic  *theLowEPionPlus;
   G4LEPionMinusInelastic *theLowEPionMinus;
   G4HEPionPlusInelastic  *theHighEPionPlus;
   G4HEPionMinusInelastic *theHighEPionMinus;
  
};
#endif
