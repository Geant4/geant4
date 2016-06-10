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
// $Id: G4VHighEnergyGenerator.hh 69717 2013-05-13 09:47:57Z gcosmo $
//      V. Uzhinsky Nov. 2012
//          Added method GetProjectileNucleus for simulation of nucleus-nucleus inter. 
//
#ifndef G4VHighEnergyGenerator_h
#define G4VHighEnergyGenerator_h 1

// Class Description
// Base class for high energy interaction models in geant4. By merit of inheriting
// from this class a high energy interaction model can be used in conjunction with
// any cascade, precompound model and evaporation phase in the
// generation of complete final states for inelastic scattering.
// Class Description - End

#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
class G4KineticTrackVector;
#include "G4ReactionProduct.hh"
#include "G4V3DNucleus.hh"

class G4VHighEnergyGenerator 
{
  public:
      G4VHighEnergyGenerator(const G4String& modelName = "High Energy Generator");
      virtual ~G4VHighEnergyGenerator();

  private:
      G4VHighEnergyGenerator(const G4VHighEnergyGenerator &right);
      const G4VHighEnergyGenerator & operator=(const G4VHighEnergyGenerator &right);
      int operator==(const G4VHighEnergyGenerator &right) const;
      int operator!=(const G4VHighEnergyGenerator &right) const;

 public:
      virtual G4V3DNucleus * GetWoundedNucleus() const = 0;
      virtual G4V3DNucleus * GetProjectileNucleus() const;  // Uzhi Nov. 2012
      virtual G4KineticTrackVector * Scatter(const G4Nucleus &theNucleus, 
                                             const G4DynamicParticle &thePrimary) = 0;
      std::pair<G4double, G4double> GetEnergyMomentumCheckLevels() const;
      void SetEnergyMomentumCheckLevels(G4double relativeLevel, G4double AbsoluteLevel);
      virtual void ModelDescription(std::ostream&) const;
      virtual G4String GetModelName() const;

  private:
      std::pair<G4double, G4double> epCheckLevels;

  private:
      G4String theGeneratorModelName;

};
#endif


