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
// $Id$
//
/// \file hbook/include/ExG4HbookFileManager.hh
/// \brief Definition of the ExG4HbookFileManager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#ifndef ExG4HbookFileManager_h
#define ExG4HbookFileManager_h 1

#include "G4VFileManager.hh"
#include "globals.hh"

#include <tools/hbook/wfile>

/// Manager class for HBook file operations.
///
/// The class implements the G4VFileManager manager for HBook.
/// It is provided separately from geant4/source/analysis in order
/// to avoid a need of linking Geant4 kernel libraries with cerblib.

class ExG4HbookFileManager : public G4VFileManager
{
  public:
    ExG4HbookFileManager(const G4AnalysisManagerState& state);
    virtual ~ExG4HbookFileManager();

    // Methods to manipulate files
    using G4VFileManager::OpenFile;
    virtual G4bool OpenFile(const G4String& fileName);
    virtual G4bool WriteFile();
    virtual G4bool CloseFile(); 
    
    // Specific methods for files per objects
    void CreateNtupleDirectory();

    // Get method
    G4bool IsFile() const;
   
  private:
    static const G4String fgkDefaultNtupleDirectoryName;

    // data members
    tools::hbook::wfile*  fFile;
};

// inline functions

inline G4bool ExG4HbookFileManager::IsFile() const
{ return fFile; }

#endif 

#endif
