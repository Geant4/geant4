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
#include "G4VisTrajContext.hh"

// Default configuration
G4VisTrajContext::G4VisTrajContext(const G4String& name)
  :fName(name)
  ,fLineColour(G4Colour::Grey())
  ,fLineVisible(true)
  ,fDrawLine(true)
  ,fDrawAuxPts(false)
  ,fAuxPtsType(G4Polymarker::squares)
  ,fAuxPtsSize(2)
  ,fAuxPtsSizeType(G4VMarker::screen)
  ,fAuxPtsFillStyle(G4VMarker::filled)
  ,fAuxPtsColour(G4Colour::Magenta())
  ,fAuxPtsVisible(true)
  ,fDrawStepPts(false)
  ,fStepPtsType(G4Polymarker::circles)
  ,fStepPtsSize(2)
  ,fStepPtsSizeType(G4VMarker::screen)
  ,fStepPtsFillStyle(G4VMarker::filled)
  ,fStepPtsColour(G4Colour::Yellow())
  ,fStepPtsVisible(true)
  ,fTimeSliceInterval(0.)
{}

// Destructor
G4VisTrajContext::~G4VisTrajContext() {}
