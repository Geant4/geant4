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
//
// $Id: G4VSensitiveDetector.cc 94771 2015-12-09 09:44:05Z gcosmo $
//
// G4VSensitiveDetector
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"

G4VSensitiveDetector::G4VSensitiveDetector(G4String name)
:verboseLevel(0),active(true),
 ROgeometry(nullptr),filter(nullptr)
{
  size_t sLast = name.last('/');
  if(sLast==std::string::npos)
  { // detector name only
    SensitiveDetectorName = name;
    thePathName = "/";
  }
  else
  { // name conatin the directory path
    SensitiveDetectorName = name;
    SensitiveDetectorName.remove(0,sLast+1);
    thePathName = name;
    thePathName.remove(sLast+1,name.length()-sLast);
    if(thePathName(0)!='/') thePathName.prepend("/");
  }
  fullPathName = thePathName + SensitiveDetectorName;
}

G4VSensitiveDetector::G4VSensitiveDetector(const G4VSensitiveDetector &right)
{
  SensitiveDetectorName = right.SensitiveDetectorName;
  thePathName = right.thePathName;
  fullPathName = right.fullPathName;
  verboseLevel = right.verboseLevel;
  active = right.active;
  ROgeometry = right.ROgeometry;
  filter = right.filter;
}

G4VSensitiveDetector::~G4VSensitiveDetector()
{
}

G4VSensitiveDetector* G4VSensitiveDetector::Clone() const
{
    G4ExceptionDescription msg;
    msg << "Derived class does not implement cloning,\n"
    << "but Clone method called.\n"
    << "Cannot continue;";
    G4Exception("G4VSensitiveDetector::Clone","Det0010",FatalException,msg);
    return nullptr;
}

G4VSensitiveDetector & G4VSensitiveDetector::operator=(const G4VSensitiveDetector &right)
{
  if (this == &right ) return *this;
  SensitiveDetectorName = right.SensitiveDetectorName;
  thePathName = right.thePathName;
  fullPathName = right.fullPathName;
  verboseLevel = right.verboseLevel;
  active = right.active;
  ROgeometry = right.ROgeometry;
  filter = right.filter;
  return *this;
}

G4int G4VSensitiveDetector::operator==(const G4VSensitiveDetector &right) const
{
   return (this==&right);
}

G4int G4VSensitiveDetector::operator!=(const G4VSensitiveDetector &right) const
{
   return (this!=&right);
}

G4int G4VSensitiveDetector::GetCollectionID(G4int i)
{
   return G4SDManager::GetSDMpointer()->GetCollectionID(SensitiveDetectorName+"/"+collectionName[i]); 
}

//----- following methoods are abstract methods to be
//----- implemented in the concrete classes

void G4VSensitiveDetector::Initialize(G4HCofThisEvent*)
{
}

void G4VSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
}

void G4VSensitiveDetector::clear()
{
}

void G4VSensitiveDetector::DrawAll()
{
}

void G4VSensitiveDetector::PrintAll()
{
}

