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
//
// $Id: G4JasSystem.cc,v 1.8.4.1 2001/06/28 19:07:47 gunter Exp $
// GEANT4 tag $Name:  $
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
