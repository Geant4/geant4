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
// $Id: MyTrackerHitsCollection.hh,v 1.4 2001-07-11 10:09:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyTrackerHitsCollection_h
#define MyTrackerHitsCollection_h 1

#include "g4rw/tvordvec.h"

#include "G4VHitsCollection.hh"
#include "MyTrackerHit.hh"

class G4VSensitiveDetector;

class MyTrackerHitsCollection : public G4VHitsCollection
{
  public:
    MyTrackerHitsCollection();
    MyTrackerHitsCollection (G4String, G4String);
    ~MyTrackerHitsCollection();
    void DrawAllHits();
    void PrintAllHits();
  private:
    G4RWTValOrderedVector<MyTrackerHit> theCollection;
  public:
    inline G4RWTValOrderedVector<MyTrackerHit>& GetVector()
    { return theCollection; };
    inline void insert(MyTrackerHit* pHit)
    { theCollection.insert(*pHit); }
};

#endif

