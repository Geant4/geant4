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
// $Id: G4VDigiCollection.hh 67992 2013-03-13 10:59:57Z gcosmo $
//

#ifndef G4VDigiCollection_h
#define G4VDigiCollection_h 1

class G4VDigi;
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

  public:
      // GetDigi and GetSize are given a default implementation here so
      // that the template G4TDigiCollection can be used, but they
      // are re-implemented in G4TDigiCollection.
      virtual G4VDigi* GetDigi(size_t) const { return 0; }
      virtual size_t GetSize() const { return 0; }

};

#endif

