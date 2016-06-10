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
// $Id: G4HadNucl3BodyMomDst.hh 67796 2013-03-08 06:18:39Z mkelsey $
// Author:  Michael Kelsey (SLAC)
// Date:    7 March 2013
//
// Description: class containing parametrized momentum distributions
//              in the CM for hadron/nucleon 3-body final states
//

#ifndef G4HadNucl3BodyMomDst_h
#define G4HadNucl3BodyMomDst_h 1

#include "G4InuclParamMomDst.hh"

class G4HadNucl3BodyMomDst : public G4InuclParamMomDst {
public:
  G4HadNucl3BodyMomDst(G4int verbose = 0);
  virtual ~G4HadNucl3BodyMomDst() {;}
};

#endif
