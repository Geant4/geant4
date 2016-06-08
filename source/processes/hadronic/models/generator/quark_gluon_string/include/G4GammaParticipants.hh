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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4GammaParticipants_h
#define G4GammaParticipants_h 1

// An implementation of a parton string model, re-using
// the most parts of quark gluon string model. The calculation
// of collision probabilities was implemented as a transparent 
// nucleus approximation. This implies that exactely ont nucleon
// in the nucleus participates in the reaction.

#include "G4QGSParticipants.hh" 

class G4GammaParticipants : public G4QGSParticipants
{
  public:
    virtual ~G4GammaParticipants(){}

  private:
    virtual G4VSplitableHadron* SelectInteractions(const G4ReactionProduct  &thePrimary);
          
};

#endif


