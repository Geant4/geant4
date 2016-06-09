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
#ifndef G4XpimNTotal_h
#define G4XpimNTotal_h

#include "G4XPDGTotal.hh"
#include "G4VCrossSectionSource.hh"
#include "G4KineticTrack.hh"
#include <utility>
#include <vector>

class G4XpimNTotal : public G4VCrossSectionSource
{
public:
  G4XpimNTotal();
  virtual ~G4XpimNTotal() {}
  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
  virtual const G4CrossSectionVector* GetComponents() const { return 0; }
  virtual G4String Name() const {return "G4XpimNTotal";}
private:
  G4XPDGTotal thePDGData;
  std::vector<std::pair<double,double> > theLowEData;
};
#endif
