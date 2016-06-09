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
// $Id: G4VEmissionProbability.hh,v 1.7 2010-11-15 20:30:26 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// JMQ (06 September 2008) Also external choices have been added for 
// superimposed Coulomb barrier (if useSICB is set true, by default is false) 


#ifndef G4VEmissionProbability_h
#define G4VEmissionProbability_h 1


#include "globals.hh"
#include "G4Fragment.hh"
#include "G4PairingCorrection.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4Pow.hh"

class G4VEmissionProbability 
{
public:
  G4VEmissionProbability();
  virtual ~G4VEmissionProbability();

private:  
  G4VEmissionProbability(const G4VEmissionProbability &right);

  const G4VEmissionProbability & operator=(const G4VEmissionProbability &right);
  G4bool operator==(const G4VEmissionProbability &right) const;
  G4bool operator!=(const G4VEmissionProbability &right) const;
  
public:
  virtual G4double EmissionProbability(const G4Fragment & fragment, const G4double anEnergy) = 0;

  // for cross section selection
  inline void SetOPTxs(G4int opt) { OPTxs = opt; }
  // for superimposed Coulomb Barrier for inverse cross sections 	
  inline void UseSICB(G4bool use) { useSICB = use; }	

protected:
  G4int OPTxs;
  G4bool useSICB;

  G4Pow*   fG4pow;
  G4PairingCorrection* fPairCorr;
  G4EvaporationLevelDensityParameter * theEvapLDPptr;

};


#endif
