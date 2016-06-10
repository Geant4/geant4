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
// $Id: G4XmlFileManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// The manager for Xml file output operations.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4XmlFileManager_h
#define G4XmlFileManager_h 1

#include "G4VFileManager.hh"
#include "G4TNtupleDescription.hh"
#include "globals.hh"

#include "tools/waxml/ntuple"

#include <fstream>
#include <memory>

class G4AnalysisManagerState;

class G4XmlFileManager : public G4VFileManager
{
  public:
    explicit G4XmlFileManager(const G4AnalysisManagerState& state);
    ~G4XmlFileManager();

    // Type aliases
    using NtupleType = tools::waxml::ntuple;
    using NtupleDescriptionType = G4TNtupleDescription<NtupleType>;

    // Methods to manipulate output files
    virtual G4bool OpenFile(const G4String& fileName) final;
    virtual G4bool WriteFile() final;
    virtual G4bool CloseFile() final; 

    // Specific methods for files per objects
    G4bool CreateHnFile();
    G4bool CloseHnFile(); 
    G4bool CreateNtupleFile(NtupleDescriptionType* ntupleDescription);
    G4bool CloseNtupleFile(NtupleDescriptionType* ntupleDescription); 
    
    // Get methods
    std::shared_ptr<std::ofstream> GetHnFile() const;
    
   private:
    // data members
    std::shared_ptr<std::ofstream> fHnFile;
};

// inline functions

inline std::shared_ptr<std::ofstream> G4XmlFileManager::GetHnFile() const
{ return fHnFile; }

#endif
