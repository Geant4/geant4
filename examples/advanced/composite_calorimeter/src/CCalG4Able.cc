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
///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Able.cc
// Description: CCalG4Able is the base class of a Geant4 geometry factory
///////////////////////////////////////////////////////////////////////////////
//Comment/Uncomment next line to unset/set debug information printing
//#define debug
//#define ddebug

#ifdef ddebug
  #include "G4Timer.hh"
#endif
#ifdef debug
  #include "CCalutils.hh"
#endif

#include "CCalGeometryConfiguration.hh"
#include "CCalG4Able.hh"
#include "CCalSensitiveConfiguration.hh"

#include "G4Color.hh"
#include "G4VisAttributes.hh"


CCalG4Able::CCalG4Able(G4String name):
  detPhysicalVolume(0), g4ableName(name), sensitivity(false),
  visProperties(CCalSensitiveConfiguration::getInstance()->getFileName(name)+".vis") {
  //Initialize g4VisAtt pointers
  for (G4int i=0; i<CCalVisualisable::TotalVisTypes; ++i) {
    g4VisAtt[i]=0;
  }
  sensitivity = 
    CCalSensitiveConfiguration::getInstance()->getSensitiveFlag(name);
}

CCalG4Able::~CCalG4Able() {
  if (detPhysicalVolume) 
    delete[] detPhysicalVolume;
}

G4VPhysicalVolume* CCalG4Able::PhysicalVolume(G4VPhysicalVolume* pv) {
  //If detPhysicalVolume is not (nil) the volume has already been built
  //so return it. In other case, construct it and its daughters, then
  //check for sensitivity and build it if set.
#ifdef ddebug
  G4Timer timer;
  timer.Start();
#endif
  if (CCalGeometryConfiguration::getInstance()->getConstructFlag(G4Name())!=0){
    if (!detPhysicalVolume) {
      detPhysicalVolume = constructIn(pv);
      for (unsigned int i = 0; i < theG4DetectorsInside.size(); i++) {
        theG4DetectorsInside[i]->PhysicalVolume(detPhysicalVolume);
      }
      if (sensitivity) {
#ifdef debug
        G4cout << "==> Making " << detPhysicalVolume->GetName() << " sensitive..." 
               << G4endl;
#endif
        constructSensitive();
      } //if sensitivity
    } //if sensitive
  } //if construct
  else {
    G4cout << "NOTE: You decided to skip the construction of " 
           << G4Name() << G4endl;
  }
#ifdef ddebug
  timer.Stop();
  G4cout << tab << "CCalG4Able::PhysicalVolume(...) --> time spent: " 
         << timer << G4endl;
#endif
  return detPhysicalVolume;
}

void CCalG4Able::AddCCalG4Able(CCalG4Able* det) {
  theG4DetectorsInside.push_back(det);
}

void CCalG4Able::setVisType(CCalVisualisable::visType vt, G4LogicalVolume* log) {
  if (!g4VisAtt[vt]) {
#ifdef debug
    G4cout << "CCalG4Able::setVisType: Constructing G4VisAttributes for " 
           << log->GetName() << " as " << vt << G4endl;
#endif
    G4Color col(visProperties.colorRed(vt),
                visProperties.colorGreen(vt),
                visProperties.colorBlue(vt));
    G4bool wf      = visProperties.isWireFrame(vt);
    G4bool visible = visProperties.isVisible(vt);
    
#ifdef debug
    G4cout << "Color: " 
           << visProperties.colorRed(vt)   << ", " 
           << visProperties.colorGreen(vt) << ", "
           << visProperties.colorBlue(vt)  << tab
           << "Wireframe: " << wf << tab
           << "Visible: " << visible << G4endl;
#endif
    g4VisAtt[vt] = new G4VisAttributes(col);
    g4VisAtt[vt]->SetForceWireframe(wf);
    g4VisAtt[vt]->SetVisibility(visible);
  }
  log->SetVisAttributes(g4VisAtt[vt]);
}



G4bool CCalG4Able::operator==(const CCalG4Able& right) const {
  return detPhysicalVolume==right.detPhysicalVolume;
}



//========================================================================
//Protected and private methods.

//========================================================================
//Global operators
std::ostream& operator<<(std::ostream& os, const CCalG4Able& det) {
  if (det.detPhysicalVolume)
    os << "Physical volume already constructed." << G4endl;
  else
    os << "Physical volume still not constructed." << G4endl;

  if (det.isSensitive())
    os << "and it is Sensitive" << G4endl;
  else
    os << "and it is not Sensitive" << G4endl;
  
  return os;
}
