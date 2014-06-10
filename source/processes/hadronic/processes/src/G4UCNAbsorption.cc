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
// $Id: G4UCNAbsorption.cc 69576 2013-05-08 13:48:13Z gcosmo $
//
///////////////////////////////////////////////////////////////////////
// UCN Absorption Class Implementation
///////////////////////////////////////////////////////////////////////
//
// File:        G4UCNAbsorption.cc
// Description: Discrete Process -- Absorption of Ultra Cold Neutrons
// Version:     1.0
// Created:     2014-05-12
// Author:      Peter Gumplinger
//              adopted from Geant4UCN by Peter Fierlinger (7.9.04) and
//              Marcin Kuzniak (21.4.06)
//              1/v energy dependent absorption cross section
//              inside materials
// Updated:
//
// mail:        gum@triumf.ca
//
///////////////////////////////////////////////////////////////////////

#include "G4UCNProcessSubType.hh"

#include "G4UCNAbsorption.hh"

//#include "G4Nucleus.hh"
//#include "G4ReactionProduct.hh"
//#include "G4NucleiPropertiesTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

// G4UCNAbsorption::operator=(const G4UCNAbsorption &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

G4UCNAbsorption::G4UCNAbsorption(const G4String& processName,G4ProcessType type)
               : G4VDiscreteProcess(processName, type)
{
  if (verboseLevel>0) G4cout << GetProcessName() << " is created " << G4endl;

  SetProcessSubType(fUCNAbsorption);
}

// G4UCNAbsorption::G4UCNAbsorption(const G4UCNAbsorpton &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

G4UCNAbsorption::~G4UCNAbsorption(){}

        ////////////
        // Methods
        ////////////

// PostStepDoIt
// -------------

G4VParticleChange*
G4UCNAbsorption::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  aParticleChange.ProposeTrackStatus(fStopAndKill);

  if ( verboseLevel > 0 ) G4cout << "UCNABSORPTION at: " 
     << aTrack.GetProperTime()/s << "s, "
     << aTrack.GetGlobalTime()/s << "s. "
     << ", after track length " << aTrack.GetTrackLength()/cm << "cm, "
           << "in volume "
           << aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()
           << G4endl;
  
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// GetMeanFreePath
// ---------------

G4double G4UCNAbsorption::GetMeanFreePath(const G4Track& aTrack,
                                          G4double ,
                                          G4ForceCondition* )
{
  G4double AttenuationLength = DBL_MAX;

  const G4Material* aMaterial = aTrack.GetMaterial();
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                     aMaterial->GetMaterialPropertiesTable();

  G4double losscs = 0.0;
  if (aMaterialPropertiesTable) {
     losscs = aMaterialPropertiesTable->GetConstProperty("ABSCS");
//     if (losscs == 0.0)
//       G4cout << "No UCN Absorption length specified" << G4endl;
  }
//  else G4cout << "No UCN Absorption length specified" << G4endl;
 
  if (losscs) {

    // Calculate a UCN absorption length for this cross section

    // *** Thermal boost ***

    // Prepare neutron

    //G4double theA = aMaterial->GetElement(0)->GetN();
    //G4double theZ = aMaterial->GetElement(0)->GetZ();

    //G4ReactionProduct 
    //  theNeutron(const_cast<G4ParticleDefinition *>(aTrack.GetDefinition()));
    //theNeutron.SetMomentum(aTrack.GetMomentum());
    //theNeutron.SetKineticEnergy(aTrack.GetKineticEnergy());
    //G4ThreeVector neuVelo = theNeutron.GetMomentum()/
    //                        aTrack.GetDefinition()->GetPDGMass());

    // Prepare properly biased thermal nucleus

    //G4double theA = aMaterial->GetElement(0)->GetN();
    //G4double theZ = aMaterial->GetElement(0)->GetZ();

    //G4double eps = 0.0001;

    //G4double eleMass =
    //             G4NucleiPropertiesTable::
    //               GetNuclearMass(static_cast<G4int>(theZ+eps),
    //                              static_cast<G4int>(theA+eps))) 
    //                              / G4Neutron::Neutron()->GetPDGMass();

    //G4Nucleus aNuc;

    //G4ReactionProduct aThermalNuc =
    //  aNuc.GetBiasedThermalNucleus(eleMass,
    //                               neuVelo,
    //                               aMaterial->GetTemperature());

    // Boost to rest system and return

    //G4ReactionProduct boosted;
    //boosted.Lorentz(theNeutron, aThermalNuc);

    //G4double vel = sqrt(2*boosted.GetKineticEnergy()/
    //                       neutron_mass_c2*c_squared);

    G4double density = aMaterial->GetTotNbOfAtomsPerVolume();

    // Calculate cross section for a constant loss 

    G4double vel = aTrack.GetVelocity();
    
    //G4cout << aTrack.GetVelocity()/meter*second << " "
    //       << vel/meter*second << "meters/second" << G4endl;
   
    // Input data is normally taken from the website:
    // http://rrdjazz.nist.gov/resources/n-lengths/list.html
    // and coresponds to 2200 m/s fast neutrons

    G4double crossect = losscs*barn*2200.*meter/second/vel;   
    
    // In principle, if one asks for the MaterialProperty incoherent cross
    // section, one could put the formula for inelastic up scattering here
    // and add the cross section to the absorption
    
    //    sigma inelastic = ... ignatovic, p. 174.
    
    // attenuation length in mm
    AttenuationLength = 1./density/crossect;

    if (verboseLevel>0) G4cout << "UCNABSORPTION with" 
        << " AttenuationLength: " << AttenuationLength/m << "m"
        << " CrossSection: " << crossect/barn << "barn" << G4endl;
  }

  return AttenuationLength;
}
