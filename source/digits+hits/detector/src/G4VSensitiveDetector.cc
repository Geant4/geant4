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
// $Id: G4VSensitiveDetector.cc,v 1.5 2001/07/13 15:00:09 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// G4VSensitiveDetector
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"

G4VSensitiveDetector::G4VSensitiveDetector(G4String name)
:verboseLevel(0),active(true),ROgeometry(0)
{
  size_t sLast = name.last('/');
  if(sLast==G4std::string::npos)
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

G4int G4VSensitiveDetector::operator==(const G4VSensitiveDetector &right) const
{
   return false;
}

G4int G4VSensitiveDetector::operator!=(const G4VSensitiveDetector &right) const
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

