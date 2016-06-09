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
// $Id: G4AlphaEvaporationProbability.hh,v 1.3 2006/06/29 20:09:45 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999) 
//



#ifndef G4AlphaEvaporationProbability_h
#define G4AlphaEvaporationProbability_h 1


#include "G4EvaporationProbability.hh"


class G4AlphaEvaporationProbability : public G4EvaporationProbability
{
public:
  // Only available constructor is default constructor
  G4AlphaEvaporationProbability();

  ~G4AlphaEvaporationProbability() {}
private:  
  // Copy constructor
  G4AlphaEvaporationProbability(const G4AlphaEvaporationProbability &right);

  const G4AlphaEvaporationProbability & operator=(const G4AlphaEvaporationProbability &right);
  G4bool operator==(const G4AlphaEvaporationProbability &right) const;
  G4bool operator!=(const G4AlphaEvaporationProbability &right) const;
  

private:

  virtual G4double CalcAlphaParam(const G4Fragment & fragment) const 
  { return 1.0 + CCoeficient(static_cast<G4double>(fragment.GetZ()-GetZ()));}
	
  virtual G4double CalcBetaParam(const G4Fragment & ) const 
  { return 0.0; }

	
  G4double CCoeficient(const G4double ) const;

  // Excitation energy levels 
  std::vector<G4double> ExcitEnergies;
  // Spin of excitation energy levels 
  std::vector<G4int> ExcitSpins;
};


#endif
