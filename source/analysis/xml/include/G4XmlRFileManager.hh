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

// The manager for Xml file input operations.

// Author: Ivana Hrivnacova, 10/09/2014  (ivana@ipno.in2p3.fr)

#ifndef G4XmlRFileManager_h
#define G4XmlRFileManager_h 1

#include "G4BaseFileManager.hh"
#include "globals.hh"

#include <tools/raxml>

#include <fstream>

class G4AnalysisManagerState;

class G4XmlRFileManager : public G4BaseFileManager
{
  public:
    explicit G4XmlRFileManager(const G4AnalysisManagerState& state);
    ~G4XmlRFileManager();

    // Methods to manipulate input files
    virtual G4bool OpenRFile(const G4String& fileName);

    // Specific methods for files per objects
    tools::raxml* GetRFile(const G4String& fileName) const;
    
   private:
    // data members
    tools::xml::default_factory*  fReadFactory;
    std::map<G4String, tools::raxml*> fRFiles;
};

#endif
