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
// $Id: G4HadronicInteraction.hh 96490 2016-04-19 06:57:04Z gcosmo $
//
// Hadronic Interaction  abstract base class
// This class is the base class for the model classes.
// It sorts out the energy-range for the models and provides
// class utilities.
// original by H.P. Wellisch
// Modified by J.L.Chuma, TRIUMF, 21-Mar-1997
// Last modified: 3-Apr-1997
// Added units to energy initialization: J.L. Chuma  04-Apr-97
// Modified by J.L.Chuma, 05-May-97  to Initialize theBlockedCounter
// Modified by J.L.Chuma, 08-Jul-97 to implement the Nucleus changes
// Adding a registry for memory management of hadronic models, HPW 22-Mar-99
// 23-Jan-2009 V.Ivanchenko move constructor and destructor to the body
//                          and reorder methods in the header 
// 29-Jun-2009 V.Ivanchenko add SampleInvariantT method
// 29-Aug-2009 V.Ivanchenko moved G4ReactionDynamics to G4InelasticInteraction,
//                          add const pointers, and added recoilEnergyThreshold
//                          member and accesors

// Class Description
// This is the base class for all hadronic interaction models in geant4.
// If you want to implement a new way of producing a final state, please,
//  inherit from here.
// Class Description - End
 
#ifndef G4HadronicInteraction_h
#define G4HadronicInteraction_h 1
 
#include "G4HadFinalState.hh"
#include "G4Material.hh"
#include "G4Nucleus.hh"
#include "G4Track.hh"
#include "G4HadProjectile.hh"

class G4HadronicInteractionRegistry;

class G4HadronicInteraction
{
public: // With description
    
  explicit G4HadronicInteraction(const G4String& modelName = "HadronicModel");
    
  virtual ~G4HadronicInteraction();
    
  virtual G4HadFinalState *ApplyYourself(const G4HadProjectile &aTrack, 
					 G4Nucleus & targetNucleus ) = 0;
  // The interface to implement for final state production code.

  virtual G4double SampleInvariantT(const G4ParticleDefinition* p, 
				    G4double plab,
				    G4int Z, G4int A);
  // The interface to implement sampling of scattering or change exchange
     
  virtual G4bool IsApplicable(const G4HadProjectile & aTrack, 
  			      G4Nucleus & targetNucleus);
 
  inline G4double GetMinEnergy() const
  { return theMinEnergy; }
    
  G4double GetMinEnergy( const G4Material *aMaterial,
			 const G4Element *anElement ) const;
   
  inline void SetMinEnergy( G4double anEnergy )
  { theMinEnergy = anEnergy; }
    
  void SetMinEnergy( G4double anEnergy, const G4Element *anElement );
    
  void SetMinEnergy( G4double anEnergy, const G4Material *aMaterial );
    
  inline G4double GetMaxEnergy() const
  { return theMaxEnergy; }
    
  G4double GetMaxEnergy( const G4Material *aMaterial,
			 const G4Element *anElement ) const;
    
  inline void SetMaxEnergy( const G4double anEnergy )
  { theMaxEnergy = anEnergy; }
    
  void SetMaxEnergy( G4double anEnergy, const G4Element *anElement );
    
  void SetMaxEnergy( G4double anEnergy, const G4Material *aMaterial );

  inline G4int GetVerboseLevel() const
  { return verboseLevel; }

  inline void SetVerboseLevel( G4int value )
  { verboseLevel = value; }

  inline const G4String& GetModelName() const
  { return theModelName; }

  void DeActivateFor(const G4Material* aMaterial);
 
  inline void ActivateFor( const G4Material *aMaterial ) 
  { 
    Block(); 
    SetMaxEnergy(GetMaxEnergy(), aMaterial);
    SetMinEnergy(GetMinEnergy(), aMaterial);
  }

  void DeActivateFor( const G4Element *anElement ); 
  inline void ActivateFor( const G4Element *anElement )
  { 
    Block(); 
    SetMaxEnergy(GetMaxEnergy(), anElement);
    SetMinEnergy(GetMinEnergy(), anElement);
  }

  G4bool IsBlocked( const G4Material *aMaterial ) const;
  G4bool IsBlocked( const G4Element *anElement) const;

  inline void SetRecoilEnergyThreshold(G4double val) 
  { recoilEnergyThreshold = val; }

  G4double GetRecoilEnergyThreshold() const 
  { return recoilEnergyThreshold;}

  virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;

  virtual std::pair<G4double, G4double> GetEnergyMomentumCheckLevels() const;

  inline void 
  SetEnergyMomentumCheckLevels(G4double relativeLevel, G4double absoluteLevel)
  { epCheckLevels.first = relativeLevel;
    epCheckLevels.second = absoluteLevel; }

  virtual void ModelDescription(std::ostream& outFile) const ; //=0;

  // Initialisation before run
  virtual void BuildPhysicsTable(const G4ParticleDefinition&);
  virtual void InitialiseModel();

protected:

  inline void SetModelName(const G4String& nam) 
  { theModelName = nam; }

  inline G4bool IsBlocked() const { return isBlocked;}
  inline void Block() { isBlocked = true; }
    
  G4HadFinalState theParticleChange;
  // the G4HadFinalState object which is modified and returned
  // by address by the ApplyYourself method,
  // (instead of aParticleChange as found in G4VProcess)
    
  G4int verboseLevel;
  // control flag for output messages
  // 0: silent
  // 1: warning messages
  // 2: more
  // (instead of verboseLevel as found in G4VProcess)

  // these two have global validity energy range    
  G4double theMinEnergy;
  G4double theMaxEnergy;

  G4bool isBlocked;

private:       

  G4HadronicInteraction(const G4HadronicInteraction &right ) = delete;
  const G4HadronicInteraction& operator=(const G4HadronicInteraction &right) = delete;
  G4bool operator==(const G4HadronicInteraction &right ) const = delete;
  G4bool operator!=(const G4HadronicInteraction &right ) const = delete;
    
  G4HadronicInteractionRegistry* registry;

  G4double recoilEnergyThreshold;

  G4String theModelName;

  std::pair<G4double, G4double> epCheckLevels;

  std::vector<std::pair<G4double, const G4Material *> > theMinEnergyList;
  std::vector<std::pair<G4double, const G4Material *> > theMaxEnergyList;
  std::vector<std::pair<G4double, const G4Element *> > theMinEnergyListElements;
  std::vector<std::pair<G4double, const G4Element *> > theMaxEnergyListElements;
  std::vector<const G4Material *> theBlockedList;
  std::vector<const G4Element *> theBlockedListElements;
};
 
#endif
