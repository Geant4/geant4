// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4JasSystem.cc,v 1.8 2000/11/16 13:44:43 barrand Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
// Guy Barrand 14 September 2000

#ifdef G4ANALYSIS_BUILD_JAS

#include "JasHistogramFactory.h"

#include "G4ios.hh"

#include "G4JasSystem.hh"

G4JasSystem::G4JasSystem (
 const G4String& aName
)
:fName(aName)
,fHistogramFactory(0)
{
}
G4JasSystem::~G4JasSystem (
)
{
  delete fHistogramFactory;
}
const G4String& G4JasSystem::GetName() const {
  return fName;
}
IHistogramFactory* G4JasSystem::GetHistogramFactory() {
  if(!fHistogramFactory) {
    G4cout << "Activate jas analysis system" << G4std::endl;
    fHistogramFactory = new JasHistogramFactory();
  }
  return fHistogramFactory;
}
void G4JasSystem::Store(IHistogram*,const G4String&) {
}
void G4JasSystem::Plot(IHistogram*) {
}

#endif
