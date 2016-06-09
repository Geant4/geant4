//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VEvaporationChannel.hh,v 1.7 2002/12/12 19:17:14 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


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

private:
  G4String Name;

};


#endif
