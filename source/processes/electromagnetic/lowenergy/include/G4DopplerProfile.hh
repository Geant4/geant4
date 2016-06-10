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
// $Id: G4DopplerProfile.hh 66241 2012-12-13 18:34:42Z gunter $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
//  6 Aug 2001   MGP        Created
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Shell data set: shell identifiers and binding energies

// -------------------------------------------------------------------

#ifndef G4DOPPLERPROFILE_HH
#define G4DOPPLERPROFILE_HH 1

#include "globals.hh"
#include <vector>
#include <map>

class G4VEMDataSet;

class G4DopplerProfile 
{ 
public:

  G4DopplerProfile(G4int minZ = 1, G4int maxZ = 100);

  ~G4DopplerProfile();
 
  size_t NumberOfProfiles(G4int Z) const;

  const G4VEMDataSet* Profiles(G4int Z) const;
  const G4VEMDataSet* Profile(G4int Z, G4int ShellIndex) const;

   void PrintData() const;

  // Random select a momentum value based on Doppler profiles
  G4double RandomSelectMomentum(G4int Z, G4int shellIndex) const;
 
private:

  // Hide copy constructor and assignment operator 
  G4DopplerProfile& operator=(const G4DopplerProfile& right);
  G4DopplerProfile(const G4DopplerProfile&);

  std::map<G4int,G4VEMDataSet*,std::less<G4int> > profileMap;
  std::vector<G4int> nShells;

  G4int zMin;
  G4int zMax; 
  
  size_t nBiggs;

  std::vector<G4double> biggsP;

  void LoadBiggsP(const G4String& fileName);
  void LoadProfile(const G4String& fileName, G4int Z);

};
 
#endif









