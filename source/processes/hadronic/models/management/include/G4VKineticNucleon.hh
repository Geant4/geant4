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
// $Id: G4VKineticNucleon.hh 80499 2014-04-24 13:54:26Z gcosmo $
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

      virtual const G4ParticleDefinition* GetDefinition()const =0;

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


