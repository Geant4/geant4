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
