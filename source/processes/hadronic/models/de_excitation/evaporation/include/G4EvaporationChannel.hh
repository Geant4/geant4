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
// $Id: G4EvaporationChannel.hh 98739 2016-08-09 12:56:55Z gcosmo $
//
//
//J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// 17-11-2010 V.Ivanchenko in constructor replace G4VEmissionProbability by 
//            G4EvaporationProbability and do not new and delete probability
//            object at each call; use G4Pow

#ifndef G4EvaporationChannel_h
#define G4EvaporationChannel_h 1

#include "G4VEvaporationChannel.hh"
#include "G4EvaporationProbability.hh"
#include "G4VCoulombBarrier.hh"

class G4PairingCorrection;

class G4EvaporationChannel : public G4VEvaporationChannel
{
public:

  explicit G4EvaporationChannel(G4int A, G4int Z, 
                                const G4String & aName,
		                G4EvaporationProbability*,
	                        G4VCoulombBarrier*);

  virtual ~G4EvaporationChannel();

  void Initialise();

  virtual G4double GetEmissionProbability(G4Fragment* fragment); 
  
  virtual G4Fragment* EmittedFragment(G4Fragment* theNucleus);

private: 
  
  G4EvaporationChannel(const G4EvaporationChannel & right) = delete;
  const G4EvaporationChannel & operator=
  (const G4EvaporationChannel & right) = delete;
  G4bool operator==(const G4EvaporationChannel & right) const = delete;
  G4bool operator!=(const G4EvaporationChannel & right) const = delete;

private:

  // This data member define the channel. 
  // They are intializated at object creation (constructor) time.

  G4int theA;
  G4int theZ;

  G4double EvapMass;
  G4double CoulombBarrier;

  // For evaporation probability calcualation
  G4EvaporationProbability * theProbability;

  // For Coulomb Barrier calculation
  G4VCoulombBarrier * theCoulombBarrier;

  // For pairing correction calculation
  G4PairingCorrection* pairingCorrection;
   
  //---------------------------------------------------

  // These values depend on the nucleus that is being evaporated.
  // They are calculated through the Initialize method which 
  // takes as parameters 
  // the atomic number, charge and excitation energy of nucleus.

  G4int ResA;
  G4int ResZ;

  G4double Mass;
	
  // Emission Probability
  G4double EmissionProbability;

  // Kinetic Energy that can be carried by fragment
  G4double MinKinEnergy;
  G4double MaxKinEnergy;

};


#endif
