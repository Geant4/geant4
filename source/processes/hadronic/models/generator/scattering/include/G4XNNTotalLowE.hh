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

#ifndef G4XNNTotalLowE_h
#define G4XNNTotalLowE_h

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include "G4CrossSectionVector.hh"
#include "G4LowEXsection.hh"
#include "g4std/map"

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
  virtual G4double HighLimit() const { return 3.*GeV; }


protected:


private:  

  G4XNNTotalLowE(const G4XNNTotalLowE &right);
  const G4XNNTotalLowE& operator=(const G4XNNTotalLowE &right);
  
  static const G4double ppTot[29];
  static const G4double ss[29];
  static const G4double npTot[29];
  static const G4int tableSize;

  G4std::map <G4String, G4LowEXsection *, G4std::less<G4String> > theCrossSections;
  typedef G4std::map <G4String, G4LowEXsection*, G4std::less<G4String> > LowEMap;

};

#endif











































