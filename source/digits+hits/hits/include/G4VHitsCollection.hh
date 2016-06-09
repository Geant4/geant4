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
// $Id: G4VHitsCollection.hh,v 1.10 2002/10/30 16:20:06 asaim Exp $
// GEANT4 tag $Name: geant4-05-01 $
//

#ifndef G4VHitsCollection_h
#define G4VHitsCollection_h 1

class G4VHit;
#include "globals.hh"

// class description:
//
//  This is the base class of hits collection. The user is advised to
// use G4THitsCollection template class in case his/her collection is
// transient. While, in case the collection is persistent with ODBMS,
// the concrete collection class can be directly derived from this
// class.
//  Geant4 kernel will use this class methods.

class G4VHitsCollection 
{
  public:
      G4VHitsCollection();
      G4VHitsCollection(G4String detName,G4String colNam);
      virtual ~G4VHitsCollection();
      G4int operator==(const G4VHitsCollection &right) const;

      virtual void DrawAllHits();
      virtual void PrintAllHits();

  protected:

      // Collection name
      G4String collectionName;
      G4String SDname;

  public:
      inline G4String GetName()
      { return collectionName; }
      inline G4String GetSDname()
      { return SDname; }

  public:
      // GetHit and GetSize are given a default implementation here so
      // that the template G4THitsCollection can be used, but they
      // are re-implemented G4THitsCollection.
      virtual G4VHit* GetHit(size_t i) const { return 0; } 
      virtual size_t GetSize() const { return 0; };

};

#endif

