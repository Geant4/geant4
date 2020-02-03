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

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef g4hdf5_defs_h
#define g4hdf5_defs_h

#include "tools/hdf5/ntuple"
#include "G4Hdf5AnalysisManager.hh"
#include "G4Hdf5AnalysisReader.hh"
#include "g4hntools_defs.hh"

namespace G4Hdf5
{
  // Ntuple types
  using G4Ntuple = tools::hdf5::ntuple; 
  using G4NtupleIterator = std::vector<tools::hdf5::ntuple*>::iterator;
  using G4NtupleConstIterator = std::vector<tools::hdf5::ntuple*>::const_iterator;

  // RNtuple types
  using G4RNtuple = tools::hdf5::ntuple; 

  // Managers
  using G4AnalysisManager = G4Hdf5AnalysisManager; 
  using G4AnalysisReader = G4Hdf5AnalysisReader; 
}  

#endif

  
