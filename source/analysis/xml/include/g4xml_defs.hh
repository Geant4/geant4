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
// $Id: g4xml_defs.hh 92688 2015-09-14 07:01:13Z gcosmo $

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifndef g4xml_defs_h
#define g4xml_defs_h

#include "tools/histo/h1d"
#include "tools/histo/h2d"
#include "tools/histo/h3d"
#include "tools/histo/p1d"
#include "tools/histo/p2d"
#include "tools/waxml/ntuple"
#include "G4XmlAnalysisManager.hh"
#include "G4XmlAnalysisReader.hh"

namespace G4Xml {

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
  using G4Ntuple = tools::waxml::ntuple; 
  using G4NtupleIterator = std::vector<tools::waxml::ntuple*>::iterator;
  using G4NtupleConstIterator = std::vector<tools::waxml::ntuple*>::const_iterator;

  // RNtuple types
  using G4RNtuple = tools::aida::ntuple; 

  // Managers
  using G4AnalysisManager = G4XmlAnalysisManager; 
  using G4AnalysisReader = G4XmlAnalysisReader; 
}  

#endif

  
