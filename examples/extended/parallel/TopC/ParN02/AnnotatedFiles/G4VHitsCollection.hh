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
/// \file G4VHitsCollection.hh
/// \brief Definition of the G4VHitsCollection class
//
//
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

//MSH_BEGIN
class G4VHitsCollection 
{
  public:
      G4VHitsCollection();
      G4VHitsCollection(G4String detName,G4String colNam);
      virtual ~G4VHitsCollection();
      G4bool operator==(const G4VHitsCollection &right) const;

      virtual void DrawAllHits();
      virtual void PrintAllHits();

  protected:

      // Collection name
      G4String collectionName;  /*MSH: predefined
      [elementGet: { $ELEMENT = $THIS->GetName(); }]
      [elementSet: { Shadowed_param->collectionName=$ELEMENT; }]
			       */

      G4String SDname;   /* MSH: predefined
      [elementGet: { $ELEMENT = $THIS->GetSDname(); }]
      [elementSet: { Shadowed_param->SDname=$ELEMENT; }]
			 */

      /* MSH_derivedclass: G4THitsCollection<ExN02TrackerHit>("","") */

  public:
      inline G4String GetName()
      { return collectionName; }
      inline G4String GetSDname()
      { return SDname; }

  public:
      // GetHit and GetSize are given a default implementation here so
      // that the template G4THitsCollection can be used, but they
      // are re-implemented G4THitsCollection.
      virtual G4VHit* GetHit(size_t) const { return 0; } 
      virtual size_t GetSize() const { return 0; };

};
//MSH_END
#endif

