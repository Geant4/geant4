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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 03.09.2008 (J.M.Quesada) for external choice of inverse cross section option
// 06.09.2008 (J.M.Quesada) external choices have been added for superimposed 
//                          Coulomb barrier (if useSICB is set true, by default is false) 
// 24.04.2010 (V.Ivanchenko) moved constructor and destructor to source; added two
//                          new virtual methods EmittedFragment(s) to allow more optimal
//                          work with G4Fragment objects
//                          

#ifndef G4VEvaporationChannel_h
#define G4VEvaporationChannel_h 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4VEvaporationChannel
{
public:

  G4VEvaporationChannel(const G4String & aName = "Anonymous");
  virtual ~G4VEvaporationChannel();

  // return emitted fragment, initial fragment is modified
  // and not deleted
  virtual G4Fragment* EmittedFragment(G4Fragment* theNucleus);

  // return vector of emitted fragments, initial fragment is modified
  // but not included in this vector
  virtual G4FragmentVector* BreakUpFragment(G4Fragment* theNucleus);

  // old method initial fragment is not modified, its copy included 
  // in the list of emitted fragments
  virtual G4FragmentVector * BreakUp(const G4Fragment & theNucleus) = 0;

  virtual G4double GetEmissionProbability(G4Fragment* theNucleus) = 0;

  inline G4String GetName() const {return Name;}
  inline void SetName(const G4String & aName) { Name = aName;}

  // for cross section selection
  inline void SetOPTxs(G4int opt) { OPTxs = opt; }
  // for superimposed Coulomb Barrier for inverse cross sections 	
  inline void UseSICB(G4bool use) { useSICB = use; }	

protected:

  G4int OPTxs;
  G4bool useSICB;

private:
  G4String Name;

  G4VEvaporationChannel(const G4VEvaporationChannel & right);
  const G4VEvaporationChannel & operator=(const G4VEvaporationChannel & right);
  G4bool operator==(const G4VEvaporationChannel & right) const;
  G4bool operator!=(const G4VEvaporationChannel & right) const;

};


#endif
