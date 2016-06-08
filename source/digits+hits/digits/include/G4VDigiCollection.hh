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
// $Id: G4VDigiCollection.hh,v 1.7 2001/07/13 15:00:12 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#ifndef G4VDigiCollection_h
#define G4VDigiCollection_h 1

#include "globals.hh"

// class description:
//
//  This is the base class of digi collection. The user is advised to
// use G4TDigiCollection template class in case his/her collection is
// transient. While, in case the collection is persistent with ODBMS,
// the concrete collection class can be directly derived from this
// class.
//  Geant4 kernel will use this class methods.

class G4VDigiCollection 
{
  public:
      G4VDigiCollection();
      G4VDigiCollection(G4String DMnam,G4String colNam);
      virtual ~G4VDigiCollection();
      G4int operator==(const G4VDigiCollection &right) const;

      virtual void DrawAllDigi();
      virtual void PrintAllDigi();

  protected:

      // Collection name
      G4String collectionName;
      G4String DMname;

  public:
      inline G4String GetName()
      { return collectionName; };
      inline G4String GetDMname()
      { return DMname; };
};

#endif

