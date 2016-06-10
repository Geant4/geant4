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
// $Id: G4He3EvaporationProbability.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 17-11-2010 V.Ivanchenko integer Z and A
//
#ifndef G4He3EvaporationProbability_h
#define G4He3EvaporationProbability_h 1


#include "G4EvaporationProbability.hh"
#include "G4He3CoulombBarrier.hh"

class G4He3EvaporationProbability : public G4EvaporationProbability
{
public:
  // Only available constructor
  G4He3EvaporationProbability();

  virtual ~G4He3EvaporationProbability();

private:  
  // Copy constructor
  G4He3EvaporationProbability(const G4He3EvaporationProbability &right);

  const G4He3EvaporationProbability & operator=(const G4He3EvaporationProbability &right);
  G4bool operator==(const G4He3EvaporationProbability &right) const;
  G4bool operator!=(const G4He3EvaporationProbability &right) const;

private:

  virtual G4double CrossSection(const  G4Fragment & fragment, G4double K);

  G4double GetOpt0(G4double K);
  G4double GetOpt12(G4double K);
  G4double GetOpt34(G4double K);

  virtual G4double CalcAlphaParam(const G4Fragment & fragment) ;
 
  virtual G4double CalcBetaParam(const G4Fragment & fragment) ;
 
  G4double CCoeficient(G4int aZ) ;
  
  //data members
   
  G4He3CoulombBarrier theCoulombBarrier;

  G4int ResidualA;
  G4int ResidualZ; 
  G4int theA;
  G4int theZ;
  G4double ResidualAthrd;
  G4int FragmentA;
  G4double FragmentAthrd;
};


#endif




