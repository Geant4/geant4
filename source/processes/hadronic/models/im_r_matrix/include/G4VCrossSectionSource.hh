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

#ifndef G4VCrossSectionSource_h
#define G4VCrossSectionSource_h

#include "globals.hh"
#include "G4HadronicException.hh"
#include "G4CrossSectionVector.hh"

class G4ParticleDefinition;
class G4KineticTrack;

class G4VCrossSectionSource 
{

public:

  // Constructors
  G4VCrossSectionSource();

  virtual ~G4VCrossSectionSource();

  G4bool operator==(const G4VCrossSectionSource &right) const;
  G4bool operator!=(const G4VCrossSectionSource &right) const;
  //  G4bool operator<(const G4VCrossSectionSource &right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const = 0;
 
  virtual const G4CrossSectionVector* GetComponents() const = 0;
  
  virtual G4String Name() const = 0;

  virtual void Print() const;
  virtual void PrintAll(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;

  virtual G4bool IsValid(G4double e) const;

  virtual G4double HighLimit() const;
  virtual G4double LowLimit() const;

protected:

  G4bool InLimits(G4double e, G4double eLow, G4double eHigh) const;

  // Determine lighter particle
  const G4ParticleDefinition* FindLightParticle(const G4KineticTrack& trk1, 
						const G4KineticTrack& trk2) const;

  // Parameterisation
  G4double FcrossX(G4double e, G4double e0, G4double sigma, 
                   G4double eParam, G4double power) const;

  const G4ParticleDefinition * FindKeyParticle(const G4KineticTrack& trk1,const G4KineticTrack& trk2) const;
  
private:  

  G4VCrossSectionSource(const G4VCrossSectionSource &right);
  G4VCrossSectionSource& operator=(const G4VCrossSectionSource &right);
  
};

#endif














