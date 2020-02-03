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

// Author: Ivana Hrivnacova, 10/09/2019  (ivana@ipno.in2p3.fr)

#ifndef g4hntools_defs_h
#define g4hntools_defs_h

#include "tools/histo/h1d"
#include "tools/histo/h2d"
#include "tools/histo/h3d"
#include "tools/histo/p1d"
#include "tools/histo/p2d"

// H1 types
using G4H1 = tools::histo::h1d;
using G4H1Iterator = std::vector<tools::histo::h1d*>::iterator;
using G4H1ConstIterator = std::vector<tools::histo::h1d*>::const_iterator;

// H2 types
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

#endif

  
