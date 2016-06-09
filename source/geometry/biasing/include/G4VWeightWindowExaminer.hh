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
// $Id: G4VWeightWindowExaminer.hh,v 1.3 2006/06/29 18:16:54 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4VWeightWindowExaminer
//
// Class description:
//
// Interface class to a weight window examiner. The examiner
// has to return the number of tracks and their weight according to
// the initial track weight and the energy of the track.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VWeightWindowExaminer_hh
#define G4VWeightWindowExaminer_hh G4VWeightWindowExaminer_hh

#include "globals.hh"
#include "G4Nsplit_Weight.hh"

class G4VWeightWindowExaminer
{

public:  // with description
  
  // constructor and destructor
  G4VWeightWindowExaminer();
  virtual ~G4VWeightWindowExaminer();
  
  

  virtual G4Nsplit_Weight Examine(G4double w, G4double energy) const = 0; 
    // evaluate the number of tracks and their weight according
    // to the initial track weight and energy

};

#endif
