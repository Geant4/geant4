//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
  ,fAuxPtsFillStyle(G4VMarker::filled)
  ,fAuxPtsColour(G4Colour::Magenta())
  ,fAuxPtsVisible(true)
  ,fDrawStepPts(false)
  ,fStepPtsType(G4Polymarker::circles)
  ,fStepPtsSize(2)
  ,fStepPtsFillStyle(G4VMarker::filled)
  ,fStepPtsColour(G4Colour::Yellow())
  ,fStepPtsVisible(true)
{}

// Destructor
G4VisTrajContext::~G4VisTrajContext() {}
