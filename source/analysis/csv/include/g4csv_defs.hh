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
// $Id: g4csv_defs.hh 84303 2014-10-13 13:14:26Z gcosmo $

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifndef g4csv_defs_h
#define g4csv_defs_h

#include "tools/histo/h1d"
#include "tools/histo/h2d"
#include "tools/histo/h3d"
#include "tools/histo/p1d"
#include "tools/histo/p2d"
#include "tools/wcsv_ntuple"
#include "G4CsvAnalysisManager.hh"
#include "G4CsvAnalysisReader.hh"

namespace G4Csv {

  // H1 types
  typedef tools::histo::h1d  G4AnaH1; // keep for backward compatibility
  typedef tools::histo::h1d  G4H1;
  typedef std::vector<tools::histo::h1d*>::iterator  G4H1Iterator;
  typedef std::vector<tools::histo::h1d*>::const_iterator  G4H1ConstIterator;

  // H2 types
  typedef tools::histo::h2d  G4AnaH2; // keep for backward compatibility    
  typedef tools::histo::h2d  G4H2;
  typedef std::vector<tools::histo::h2d*>::iterator  G4H2Iterator;
  typedef std::vector<tools::histo::h2d*>::const_iterator  G4H2ConstIterator;

  // H3 types
  typedef tools::histo::h3d  G4H3;    
  typedef std::vector<tools::histo::h3d*>::iterator  G4H3Iterator;
  typedef std::vector<tools::histo::h3d*>::const_iterator  G4H3ConstIterator;

  // P1 types
  typedef tools::histo::p1d  G4P1;
  typedef std::vector<tools::histo::p1d*>::iterator  G4P1Iterator;
  typedef std::vector<tools::histo::p1d*>::const_iterator  G4P1ConstIterator;

  // P2 types
  typedef tools::histo::p2d  G4P2;    
  typedef std::vector<tools::histo::p2d*>::iterator  G4P2Iterator;
  typedef std::vector<tools::histo::p2d*>::const_iterator  G4P2ConstIterator;

  // Ntuple types
  typedef tools::wcsv::ntuple  G4Ntuple; 
  typedef std::vector<tools::wcsv::ntuple*>::iterator  G4NtupleIterator;
  typedef std::vector<tools::wcsv::ntuple*>::const_iterator  G4NtupleConstIterator;

  // Managers
  typedef G4CsvAnalysisManager G4AnalysisManager; 
  typedef G4CsvAnalysisReader  G4AnalysisReader; 
}  

#endif

  
