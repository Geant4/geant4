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

#ifndef G4XAnnihilationChannel_h
#define G4XAnnihilationChannel_h

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include "G4CrossSectionVector.hh"
#include "G4Clebsch.hh"
#include "G4ResonanceNames.hh"

#include <map>

class G4ParticleDefinition;
class G4PhysicsVector;
class G4KineticTrack;
class G4ResonanceWidth;
class G4ResonancePartialWidth;
class G4PartialWidthTable;

class G4XAnnihilationChannel : public G4VCrossSectionSource
{
public:

  G4XAnnihilationChannel();

  G4XAnnihilationChannel(const G4ParticleDefinition* resDefinition,
			 const G4ResonanceWidth& resWidths,
			 const G4ResonancePartialWidth& resPartWidths,
			 const G4String& partWidthLabel);

  virtual ~G4XAnnihilationChannel();

  G4bool operator==(const G4XAnnihilationChannel &right) const;
  G4bool operator!=(const G4XAnnihilationChannel &right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
 
  virtual const G4CrossSectionVector* GetComponents() const { return 0; }

  virtual G4bool IsValid(G4double e) const;

  virtual G4String Name() const;


protected:

private:  

  G4XAnnihilationChannel(const G4XAnnihilationChannel &right);
  const G4XAnnihilationChannel& operator=(const G4XAnnihilationChannel &right);
  
  const G4double Branch(const G4KineticTrack& trk1, 
			const G4KineticTrack& trk2) const;

  const G4double VariableWidth(const G4KineticTrack& trk1, 
			       const G4KineticTrack& trk2) const;
 
  const G4double VariablePartialWidth(const G4KineticTrack& trk1, 
				     const G4KineticTrack& trk2) const;

  const G4double NormalizedClebsch(const G4KineticTrack& trk1, 
				   const G4KineticTrack& trk2) const;

  G4double lowLimit;
  G4double highLimit;

  G4Clebsch clebsch;
  G4ResonanceNames theNames;
  
  // Owned pointers
  G4PhysicsVector* widthTable;
  G4PhysicsVector* partWidthTable;

  // Unowned pointer
  const G4ParticleDefinition* resonance;
};

#endif


















