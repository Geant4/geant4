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
// $Id: G4VEvaporationChannel.hh,v 1.4 2008/09/19 13:32:54 ahoward Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// JMQ (06 September 2008) Also external choices have been added for 
// superimposed Coulomb barrier (if useSICB is set true, by default is false) 



#ifndef G4VEvaporationChannel_h
#define G4VEvaporationChannel_h 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4VEvaporationChannel
{
public:
  G4VEvaporationChannel() : Name("Anonymous") {};
  G4VEvaporationChannel(const G4String & aName) : Name(aName) {};
  G4VEvaporationChannel(const G4String * aName) : Name(*aName) {};
  virtual ~G4VEvaporationChannel() {};

private:
  G4VEvaporationChannel(const G4VEvaporationChannel & right);

  const G4VEvaporationChannel & operator=(const G4VEvaporationChannel & right);
public:
  G4bool operator==(const G4VEvaporationChannel & right) const;
  G4bool operator!=(const G4VEvaporationChannel & right) const;

public:
  virtual void Initialize(const G4Fragment & fragment) = 0;

  virtual G4FragmentVector * BreakUp(const G4Fragment & theNucleus) = 0;

  virtual G4double GetEmissionProbability(void) const = 0;


  G4String GetName() const {return Name;}
  void SetName(const G4String & aName) { Name = aName;}

  // for cross section selection
  inline void SetOPTxs(G4int opt) { OPTxs = opt; }
  // for superimposed Coulomb Barrier for inverse cross sections 	
  inline void UseSICB(G4bool use) { useSICB = use; }	
protected:
  G4int OPTxs;
  G4bool useSICB;


private:
  G4String Name;

};


#endif
