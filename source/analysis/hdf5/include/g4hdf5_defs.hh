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

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef g4hdf5_defs_h
#define g4hdf5_defs_h

#include "tools/histo/h1d"
#include "tools/histo/h2d"
#include "tools/histo/h3d"
#include "tools/histo/p1d"
#include "tools/histo/p2d"
#include "tools/hdf5/ntuple"
#include "G4Hdf5AnalysisManager.hh"
#include "G4Hdf5AnalysisReader.hh"

namespace G4Hdf5 {

  // H1 types
  using G4AnaH1 = tools::histo::h1d; // keep for backward compatibility
  using G4H1 = tools::histo::h1d;
  using G4H1Iterator  = std::vector<tools::histo::h1d*>::iterator;
  using G4H1ConstIterator = std::vector<tools::histo::h1d*>::const_iterator;

  // H2 types
  using G4AnaH2 = tools::histo::h2d; // keep for backward compatibility    
  using G4H2 = tools::histo::h2d;
  using G4H2Iterator = std::vector<tools::histo::h2d*>::iterator;
  using G4H2ConstIterator = std::vector<tools::histo::h2d*>::const_iterator;

  // H3 types
  using G4H3 = tools::histo::h3d;    
  using G4H3Iterator = std::vector<tools::histo::h3d*>::iterator;
  using G4H3ConstIterator = std::vector<tools::histo::h3d*>::const_iterator;

  // P1 types
  using G4P1 = tools::histo::p1d;
  using G4P1Iterator = std::vector<tools::histo::p1d*>::iterator;
  using G4P1ConstIterator = std::vector<tools::histo::p1d*>::const_iterator;

  // P2 types
  using G4P2 = tools::histo::p2d;    
  using G4P2Iterator = std::vector<tools::histo::p2d*>::iterator;
  using G4P2ConstIterator = std::vector<tools::histo::p2d*>::const_iterator;

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

  
