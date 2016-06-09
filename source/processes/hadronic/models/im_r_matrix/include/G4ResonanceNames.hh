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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4ResonanceNames
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4RESONANCENAMES_HH
#define G4RESONANCENAMES_HH

#include "globals.hh"
#include <vector>
#include <map>

class G4ParticleDefinition;

class G4ResonanceNames 
{

public:

  G4ResonanceNames();

  ~G4ResonanceNames();

  G4bool operator==(const G4ResonanceNames& right) const;
  G4bool operator!=(const G4ResonanceNames& right) const;

  const std::vector<G4String> NstarNames() const { return nameNstar; }

  const std::vector<G4String> DeltastarNames() const { return nameDeltastar; }

  const std::vector<G4String> DeltaNames() const { return nameDelta; }

  const std::vector<G4String> LambdaNames() const { return nameLambda; }

  const std::vector<G4String> SigmaNames() const { return nameSigma; }

  const std::vector<G4String> XiNames() const { return nameSigma; }

  const G4String ShortName(const G4String& name);

  G4double MinMass(const G4String& name);

protected:

  
private:  

  G4ResonanceNames(const G4ResonanceNames& right);
  G4ResonanceNames& operator=(const G4ResonanceNames& right);

  std::vector<G4String> nameNstar;
  std::vector<G4String> nameDeltastar;
  std::vector<G4String> nameDelta;
  std::vector<G4String> nameLambda;
  std::vector<G4String> nameSigma;
  std::vector<G4String> nameXi;

  // Lowest resonance among each category (N*, Delta, Lambda, Sigma, Xi)
  std::map<G4String, G4ParticleDefinition*, std::less<G4String> > lowResMap;

  // Resonance short names
  std::map<G4String, G4String, std::less<G4String> > shortMap;

};
  
#endif
