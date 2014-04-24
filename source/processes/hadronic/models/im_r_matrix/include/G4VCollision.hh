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
// $Id: G4VCollision.hh,v 1.2 2006-06-29 20:36:04 gunter Exp $ //

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
  void establish_G4MT_TLS_G4VCollision();
  G4VCollision(void *s1, void *s2, void *s3, void *s4, void *s5, void *s6, void *s7);

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

  G4int GetNumberOfPartons(const G4ParticleDefinition * aP) const
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


















