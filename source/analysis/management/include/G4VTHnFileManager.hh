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

// Templated base class for histogram/profile file output.
//
// Author: Ivana Hrivnacova, 15/09/2020  (ivana@ipno.in2p3.fr)

#ifndef G4VTHnFileManager_h
#define G4VTHnFileManager_h 1

#include "G4BaseFileManager.hh"
#include "globals.hh"

class G4AnalysisManagerState;

template <typename HT>
class G4VTHnFileManager
{
  public:
    G4VTHnFileManager() {}
    virtual ~G4VTHnFileManager() {}

    // Methods for writing objects
    // Write a single object to an extra file (the file is closed after write)
    virtual G4bool WriteExtra(HT* ht, const G4String& htName, const G4String& fileName) = 0;
    // Write to an open file  (handled with OpenFile()/CloseFile methods)
    virtual G4bool Write(HT* ht, const G4String& htName, G4String& fileName) = 0;
};

#endif
