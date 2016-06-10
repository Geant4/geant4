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

#ifndef G4XNNTotalLowE_h
#define G4XNNTotalLowE_h

#include <map>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include "G4CrossSectionVector.hh"
#include "G4LowEXsection.hh"

class G4KineticTrack;

class G4XNNTotalLowE : public G4VCrossSectionSource
{

public:

  G4XNNTotalLowE();

  virtual ~G4XNNTotalLowE();

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
  virtual const G4CrossSectionVector* GetComponents() const { return 0; }
  virtual G4bool IsValid(G4double e) const;
  
  virtual G4String Name() const;
  virtual G4double HighLimit() const { return 3.*CLHEP::GeV; }


protected:


private:  

  G4XNNTotalLowE(const G4XNNTotalLowE &right);
  const G4XNNTotalLowE& operator=(const G4XNNTotalLowE &right);
  
  static const G4double ppTot[29];
  static const G4double ss[29];
  static const G4double npTot[29];
  static const G4int tableSize;

  std::map <const G4ParticleDefinition *, G4LowEXsection *,
  std::less<const G4ParticleDefinition *> > theCrossSections;
  typedef std::map <const G4ParticleDefinition *, G4LowEXsection*, std::less<const G4ParticleDefinition *> > LowEMap;

};

#endif











































