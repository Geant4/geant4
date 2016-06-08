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

#ifndef G4PARTICLETYPECONVERTER_HH
#define G4PARTICLETYPECONVERTER_HH

#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"

#include "g4std/map"

class G4ParticleTypeConverter
{
public:

  enum GenericType { NUCLEON, 
		     N1440, N1520, N1535, N1650, N1675, N1680, N1700, N1710, N1720, N1900, N1990, N2090, N2190, N2220, N2250,
		     D1232, D1600, D1620, D1700, D1900, D1905, D1910, D1920, D1930, D1950,
		     L1405, L1520, L1600, L1670, L1690, L1800, L1810, L1820, L1830, L1890, L2100, L2110,
		     Sigma, S1385, S1660, S1670, S1750, S1775, S1915, S1940, S2030,
		     X1530, X1690, X1820, X1950, X2030,
		     GAMMA, PION, KAON, ETA, RHO, omega, Lambda, UNKNOWN };

  G4ParticleTypeConverter();

  GenericType GetGenericType(const G4ParticleDefinition* const aParticleDef);
  GenericType GetGenericType(const G4KineticTrack& aTrack);
  GenericType GetGenericType(const G4String& aParticleName);

  G4int GetUrqmdItyp(GenericType gType);
  G4int GetUrqmdItyp(const G4ParticleDefinition* aParticleDef);

  const G4ParticleDefinition* FindIso3State(const GenericType gType, const G4int isospin3);

private:

  typedef G4std::vector<G4std::pair<const G4ParticleDefinition*, GenericType> >::const_iterator MapIterator; 
  G4std::vector<G4std::pair<const G4ParticleDefinition*, GenericType> > defMap;

};


#endif




