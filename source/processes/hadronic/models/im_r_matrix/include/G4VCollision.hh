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
// $Id: G4VCollision.hh,v 1.1 2003-10-07 12:37:31 hpw Exp $ //

#ifndef G4VCollision_h
#define G4VCollision_h

#include "globals.hh"
#include "G4CollisionVector.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"

class G4KineticTrack;

class G4VCollision 
{
public:

  G4VCollision();

  virtual ~G4VCollision();

  G4bool operator==(const G4VCollision &right) const;
  G4bool operator!=(const G4VCollision &right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, 
				const G4KineticTrack& trk2) const;

  virtual G4KineticTrackVector* FinalState(const G4KineticTrack& trk1, 
					   const G4KineticTrack& trk2) const = 0;

  virtual G4bool IsInCharge(const G4KineticTrack& trk1, 
			    const G4KineticTrack& trk2) const = 0;

  virtual G4String GetName() const = 0;

  virtual void Print() const;
  virtual void Print(const G4KineticTrack& trk1, 
		     const G4KineticTrack& trk2) const;
protected:

  G4int GetNumberOfPartons(G4ParticleDefinition * aP) const
  {
    G4int result = 0;
    for(G4int i=0; i<6; i++) 
    {
      result += aP->GetQuarkContent(i+1);
      result += aP->GetAntiQuarkContent(i+1);
    }
    return result;
  }
  
  virtual const G4CollisionVector* GetComponents() const { return 0;}

  virtual const G4VCrossSectionSource* GetCrossSectionSource() const = 0;

  virtual const G4VAngularDistribution* GetAngularDistribution() const = 0;

  virtual const std::vector<G4String>& GetListOfColliders(G4int whichOne) const = 0;


private:  

  G4VCollision(const G4VCollision &right);

  const G4VCollision& operator=(const G4VCollision &right);

};

#endif


















