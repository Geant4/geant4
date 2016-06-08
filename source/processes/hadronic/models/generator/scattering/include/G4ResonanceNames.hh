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
#include "g4std/vector"
#include "g4std/map"

class G4ParticleDefinition;

class G4ResonanceNames 
{

public:

  G4ResonanceNames();

  ~G4ResonanceNames();

  G4bool operator==(const G4ResonanceNames& right) const;
  G4bool operator!=(const G4ResonanceNames& right) const;

  const G4std::vector<G4String> NstarNames() const { return nameNstar; }

  const G4std::vector<G4String> DeltastarNames() const { return nameDeltastar; }

  const G4std::vector<G4String> DeltaNames() const { return nameDelta; }

  const G4std::vector<G4String> LambdaNames() const { return nameLambda; }

  const G4std::vector<G4String> SigmaNames() const { return nameSigma; }

  const G4std::vector<G4String> XiNames() const { return nameSigma; }

  const G4String ShortName(const G4String& name);

  G4double MinMass(const G4String& name);

protected:

  
private:  

  G4ResonanceNames(const G4ResonanceNames& right);
  G4ResonanceNames& operator=(const G4ResonanceNames& right);

  G4std::vector<G4String> nameNstar;
  G4std::vector<G4String> nameDeltastar;
  G4std::vector<G4String> nameDelta;
  G4std::vector<G4String> nameLambda;
  G4std::vector<G4String> nameSigma;
  G4std::vector<G4String> nameXi;

  // Lowest resonance among each category (N*, Delta, Lambda, Sigma, Xi)
  G4std::map<G4String, G4ParticleDefinition*, G4std::less<G4String> > lowResMap;

  // Resonance short names
  G4std::map<G4String, G4String, G4std::less<G4String> > shortMap;

};
  
#endif
