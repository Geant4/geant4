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
// $Id: G4VKineticNucleon.hh,v 1.2 2003/10/07 12:37:00 hpw Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
#ifndef G4VKineticNucleon_h
#define G4VKineticNucleon_h 1

#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"

class G4KineticTrackVector;

class G4VKineticNucleon 
{

  public:

      G4VKineticNucleon();

      G4VKineticNucleon(const G4VKineticNucleon &right);

      virtual ~G4VKineticNucleon();

      const G4VKineticNucleon& operator=(const G4VKineticNucleon& right);

      int operator==(const G4VKineticNucleon& right) const;

      int operator!=(const G4VKineticNucleon& right) const;

      virtual G4KineticTrackVector* Decay();

      virtual const G4LorentzVector& Get4Momentum() const =0;

      virtual G4ParticleDefinition* GetDefinition()const =0;

      virtual const G4ThreeVector& GetPosition() const =0;

  protected:

  private:
};

// Class G4VKineticNucleon 



inline G4KineticTrackVector* G4VKineticNucleon::Decay()
{
  return NULL;
}

inline const G4VKineticNucleon& G4VKineticNucleon::operator=(const G4VKineticNucleon& )
{
	return *this;
}
#endif


