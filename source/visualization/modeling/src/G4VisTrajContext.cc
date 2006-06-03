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
