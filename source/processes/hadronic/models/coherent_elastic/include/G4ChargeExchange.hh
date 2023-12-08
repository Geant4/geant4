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
// G4 Model: Charge and strangness exchange based on G4LightMedia model
//           28 May 2006 V.Ivanchenko
//
// Modified:
// 25-Jul-06 V.Ivanchenko add 19 MeV low energy, below which S-wave is sampled
//

#ifndef G4ChargeExchange_h
#define G4ChargeExchange_h 1

// Class Description
// Final state production model for hadron nuclear coherent charge exchange;
// Class Description - End

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"

class G4ParticleDefinition;
class G4ChargeExchangeXS;

class G4ChargeExchange : public G4HadronicInteraction
{
public:

  explicit G4ChargeExchange(G4ChargeExchangeXS*);
  ~G4ChargeExchange() override = default;

  G4ChargeExchange( const G4ChargeExchange &right ) = delete;
  const G4ChargeExchange & operator=( const G4ChargeExchange &right ) = delete;
  G4bool operator==( const G4ChargeExchange &right ) const = delete;
  G4bool operator!=( const G4ChargeExchange &right ) const = delete;

  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                 G4Nucleus& targetNucleus) override;

  G4double SampleT(const G4ParticleDefinition* theSec, const G4int A,
		   const G4double tmax) const;

private:

  G4ChargeExchangeXS* fXSection;

  G4int secID;  // Creator model ID for the secondaries created by this model
  G4double lowEnergyLimit; // lowest limit to avoid numerical problems
};

#endif
