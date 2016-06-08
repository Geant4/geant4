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
// $Id: PersEx02TrackerHit.ddl,v 1.5 2001/07/11 09:58:14 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#ifndef PersEx02TrackerHit_h
#define PersEx02TrackerHit_h 1

#include "G4PVHit.hh"
#include "G4ThreeVector.hh"
#include "G4PersistentTypes.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class PersEx02TrackerHit
 : public G4PVHit
{
  public:

      PersEx02TrackerHit();
      ~PersEx02TrackerHit();
      PersEx02TrackerHit(const PersEx02TrackerHit &right);
      const PersEx02TrackerHit& operator=(const PersEx02TrackerHit &right);
      int operator==(const PersEx02TrackerHit &right) const;

      void Draw();
      void Print();

  private:
      G4Pdouble edep;
      G4ThreeVector pos;

  public:
      inline void SetEdep(G4double de)
      { edep = de; };
      inline G4double GetEdep()
      { return edep; };
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; };
      inline G4ThreeVector GetPos()
      { return pos; };

};

#endif


