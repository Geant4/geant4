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
//      V. Uzhinsky Nov. 2012
//  Added method GetProjectileNucleus for simulation of nucleus-nucleus inter. 
//
#ifndef G4VHighEnergyGenerator_h
#define G4VHighEnergyGenerator_h 1

// Class Description
// Base class for high energy interaction models in geant4. 
// By merit of inheriting from this class a high energy 
// interaction model can be used in conjunction with
// any cascade, precompound model and evaporation phase in the
// generation of complete final states for inelastic scattering.
// Class Description - End

#include "G4HadronicInteraction.hh"
#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
#include "G4ReactionProduct.hh"
#include "G4V3DNucleus.hh"

class G4KineticTrackVector;

class G4VHighEnergyGenerator : public G4HadronicInteraction
{
  public:
      G4VHighEnergyGenerator(const G4String& modelName = "HighEnergyGenerator");
      ~G4VHighEnergyGenerator() override;

      G4VHighEnergyGenerator(const G4VHighEnergyGenerator &right) = delete;
      const G4VHighEnergyGenerator & operator=(const G4VHighEnergyGenerator &right) = delete;
      G4bool operator==(const G4VHighEnergyGenerator &right) const = delete;
      G4bool operator!=(const G4VHighEnergyGenerator &right) const = delete;

      virtual G4V3DNucleus * GetWoundedNucleus() const = 0;
      virtual G4V3DNucleus * GetProjectileNucleus() const;
      virtual G4KineticTrackVector * Scatter(const G4Nucleus &theNucleus, 
                                             const G4DynamicParticle &thePrimary) = 0;
      void ModelDescription(std::ostream&) const override;

};
#endif


