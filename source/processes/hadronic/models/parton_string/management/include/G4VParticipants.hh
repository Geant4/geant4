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
// $Id: G4VParticipants.hh 100828 2016-11-02 15:25:59Z gcosmo $
//

#ifndef G4VParticipants_h
#define G4VParticipants_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4VParticipants ----------------
//             by Gunter Folger, May 1998.
//      abstract class finding participants in a hadron Nucleus collision
//       in Parton String Models.
// ------------------------------------------------------------
// 20110805  M. Kelsey -- Move #include, Init() and SetNucleus() to .cc file

#include "globals.hh"

class G4V3DNucleus;


class G4VParticipants 
{
  public:
    G4VParticipants();
    G4VParticipants(const G4VParticipants &right);
    virtual ~G4VParticipants();

    const G4VParticipants & operator=(const G4VParticipants &right);
    int operator==(const G4VParticipants &right) const;
    int operator!=(const G4VParticipants &right) const;

    virtual void Init(G4int theZ, G4int theA);
    virtual void SetNucleus(G4V3DNucleus* aNucleus);
    virtual G4V3DNucleus* GetWoundedNucleus() const;

    virtual void InitProjectileNucleus(G4int theZ, G4int theA);
    virtual void SetProjectileNucleus(G4V3DNucleus* aNucleus);
    virtual G4V3DNucleus* GetProjectileNucleus() const;

    G4V3DNucleus* theNucleus;
    G4V3DNucleus* theProjectileNucleus;      
  private:
  
};

// Class G4VParticipants 

inline G4V3DNucleus * G4VParticipants::GetWoundedNucleus() const
{
  return theNucleus;
}

inline G4V3DNucleus * G4VParticipants::GetProjectileNucleus() const 
{
  return theProjectileNucleus;
}

#endif

