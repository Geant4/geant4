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
//
// $Id: G4VDataSetAlgorithm.hh,v 1.2 2001-10-08 07:45:36 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Base class for the generation of final state in electromagnetic processes
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4VFINALSTATE_HH
#define G4VFINALSTATE_HH 1

#include "globals.hh"

class G4VFinalState {
 
public:

  G4VFinalState() { }

  virtual ~G4VFinalState() { }
 
  virtual G4std::vector<G4DynamicParticle*>* Generate(const G4Track& aTrack,
						      const G4Step& aStep) const = 0;

private:
  
  // Hide copy constructor and assignment operator
  G4VFinalState(const G4VFinalState&);
  G4VFinalState & operator=(const G4VFinalState &right);

};
 
#endif
