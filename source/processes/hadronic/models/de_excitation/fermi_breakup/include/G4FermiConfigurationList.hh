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
// $Id: G4FermiConfigurationList.hh,v 1.5 2006-06-29 20:12:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  {
    std::for_each(Configurations.begin(),Configurations.end(),
                  DeleteConfiguration()); 
  }
  
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

  struct DeleteConfiguration
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };


};


#endif


