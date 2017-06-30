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
// $Id: G4ImportanceAlgorithm.hh 102994 2017-03-07 16:31:28Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4ImportanceAlgorithm
//
// Class description:
//
// An implementation of a G4VImportanceAlgorithm (see description in
// G4VImportanceAlgorithm.hh).

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ImportanceAlgorithm_hh
#define G4ImportanceAlgorithm_hh G4ImportanceAlgorithm_hh

#include "G4VImportanceAlgorithm.hh"

class G4ImportanceAlgorithm : public G4VImportanceAlgorithm
{

public:  // with description

  G4ImportanceAlgorithm();
    // simple construction

  virtual ~G4ImportanceAlgorithm();
    // repeate warning if triggered

  virtual G4Nsplit_Weight Calculate(G4double ipre, 
				    G4double ipost, 
				    G4double init_w) const;
    // calculate the number of tracks and their weight according 
    // to the pre and post importance value and the weight of
    // the mother track.  
private:

  void Error(const G4String &m) const;
  void Warning(const G4String &m) const;

private:

  mutable G4bool fWorned;

#ifdef G4MULTITHREADED
  static G4Mutex ImportanceMutex;
#endif

};

#endif
