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
  const G4double FcrossX(G4double e, G4double e0, G4double sigma, 
			 G4double eParam, G4double power) const;

  G4String FindKeyParticle(const G4KineticTrack& trk1,const G4KineticTrack& trk2) const;
  
  // Transverse pion mass
  const G4double GetTransversePionMass() const;

  // Min string mass
  const G4double GetMinStringMass() const;

private:  

  G4VCrossSectionSource(const G4VCrossSectionSource &right);
  G4VCrossSectionSource& operator=(const G4VCrossSectionSource &right);
  
};

#endif














