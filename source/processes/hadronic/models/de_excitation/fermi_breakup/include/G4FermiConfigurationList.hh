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
// $Id: G4FermiConfigurationList.hh,v 1.2 2003/11/20 09:46:23 jwellisc Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4FermiConfigurationList_h
#define G4FermiConfigurationList_h 1

#include "globals.hh"
#include "G4FermiConfiguration.hh"
#include "Randomize.hh"
#include "G4FermiFragmentsPool.hh"

class G4FermiConfigurationList 
{
public:
  G4FermiConfigurationList();

  ~G4FermiConfigurationList()
  {};
  
private:
  G4FermiConfigurationList(const G4FermiConfigurationList &right);
  
  const G4FermiConfigurationList & operator=(const G4FermiConfigurationList &right);
  G4bool operator==(const G4FermiConfigurationList &right) const;
  G4bool operator!=(const G4FermiConfigurationList &right) const;
  
public:

  G4bool Initialize(const G4int A, const G4int Z, const G4double TotalEnergyRF);

  G4FermiConfiguration ChooseConfiguration(void);

  G4FermiFragmentsPool & GetFragmentsPoolInstance();

private:


  enum {MaxNumOfFragments = 16};

  std::vector<G4double> NormalizedWeights;
  
  std::vector<G4FermiConfiguration*> Configurations;

};


#endif


