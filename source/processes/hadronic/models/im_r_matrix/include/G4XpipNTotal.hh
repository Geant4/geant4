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
#ifndef G4XpipNTotal_h
#define G4XpipNTotal_h

#include "G4XPDGTotal.hh"
#include "G4VCrossSectionSource.hh"
#include "G4KineticTrack.hh"
#include "G4Pair.hh"
#include <vector>

class G4XpipNTotal : public G4VCrossSectionSource
{
public:
  G4XpipNTotal();
  virtual ~G4XpipNTotal() {}
  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
  virtual const G4CrossSectionVector* GetComponents() const { return 0; }
  virtual G4String Name() const {return "G4XpipNTotal";}
private:
  G4XPDGTotal thePDGData;
  std::vector<G4Pair<double,double> > theLowEData;
};
#endif
