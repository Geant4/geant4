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
// $Id: G4NeutronEvaporationProbability.hh,v 1.15 2010-11-17 11:06:03 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 17-11-2010 V.Ivanchenko integer Z and A
//

#ifndef G4NeutronEvaporationProbability_h
#define G4NeutronEvaporationProbability_h 1


#include "G4EvaporationProbability.hh"
#include "G4NeutronCoulombBarrier.hh"

class G4NeutronEvaporationProbability : public G4EvaporationProbability
{
public:
 
  G4NeutronEvaporationProbability();
		
  virtual ~G4NeutronEvaporationProbability();

private:  
  
  G4NeutronEvaporationProbability(const G4NeutronEvaporationProbability &right);

  const G4NeutronEvaporationProbability & operator=(const G4NeutronEvaporationProbability &right);
  G4bool operator==(const G4NeutronEvaporationProbability &right) const;
  G4bool operator!=(const G4NeutronEvaporationProbability &right) const;

private:

  virtual G4double CrossSection(const  G4Fragment & fragment, G4double K);

  G4double GetOpt12(G4double K);
  G4double GetOpt34(G4double K);

  virtual G4double CalcAlphaParam(const G4Fragment & fragment);
 
  virtual G4double CalcBetaParam(const G4Fragment & fragment);
 
  //data members

  G4NeutronCoulombBarrier theCoulombBarrier; 

  G4int ResidualA;
  G4int ResidualZ; 
  G4int theA;
  G4int theZ;
  G4double ResidualAthrd;
  G4int FragmentA;
  G4double FragmentAthrd;
};


#endif
