// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSensitiveDetector.cc,v 1.1 1999/01/07 16:06:27 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// G4VSensitiveDetector
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"

G4VSensitiveDetector::G4VSensitiveDetector(G4String name)
:verboseLevel(0),active(true),ROgeometry(NULL)
{
  size_t sLast = name.last('/');
  if(sLast==RW_NPOS)
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
}

G4VSensitiveDetector::~G4VSensitiveDetector()
{
}

const G4VSensitiveDetector & G4VSensitiveDetector::operator=(const G4VSensitiveDetector &right)
{
    return *this;
}

int G4VSensitiveDetector::operator==(const G4VSensitiveDetector &right) const
{
   return false;
}

int G4VSensitiveDetector::operator!=(const G4VSensitiveDetector &right) const
{
   return true;
}

G4int G4VSensitiveDetector::GetCollectionID(G4int i)
{
   return G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[i]); 
}

//----- following methoods are abstract methods to be
//----- implemented in the concrete classes

void G4VSensitiveDetector::Initialize(G4HCofThisEvent*HCE)
{
}

void G4VSensitiveDetector::EndOfEvent(G4HCofThisEvent*HCE)
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

