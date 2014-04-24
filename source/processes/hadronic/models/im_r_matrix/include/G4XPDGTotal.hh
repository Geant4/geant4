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
//      File name:     G4XPDGTotal
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4XPDGTOTAL_HH
#define G4XPDGTOTAL_HH

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include <algorithm>
#include <map>
#include <vector>

class G4KineticTrack;
class G4ParticleDefinition;

typedef std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> G4pDefPair;

class G4XPDGTotal : public G4VCrossSectionSource
{

public:

  G4XPDGTotal();

  virtual ~G4XPDGTotal();

  G4bool operator==(const G4XPDGTotal &right) const;
  G4bool operator!=(const G4XPDGTotal &right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
 
  virtual const G4CrossSectionVector* GetComponents() const { return 0; }

  virtual G4String Name() const;

  virtual G4bool IsValid(G4double e) const;

  virtual G4double LowLimit() const { return _lowLimit; }


protected:

private:  

  G4double PDGTotal(G4double coeff, G4double ecm) const;

  G4XPDGTotal(const G4XPDGTotal &right);
  const G4XPDGTotal& operator=(const G4XPDGTotal &right);
  
  static const G4double _lowLimit;
  static const G4double _highLimit;

  static const G4int nFit;
  static const G4double ppPDGFit[5];
  static const G4double npPDGFit[5];
  static const G4double pipPDGFit[5];
  static const G4double KpPDGFit[5];
  static const G4double KnPDGFit[5];
  static const G4double gammapPDGFit[5];
  static const G4double gammagammaPDGFit[5];

  std::map <G4pDefPair, std::vector<G4double>, std::less<G4pDefPair> > xMap;
  typedef std::map <G4pDefPair, std::vector<G4double>, std::less<G4pDefPair> > PairDoubleMap;

  //  G4double eMinFit;
  //  G4double eMaxFit;
  //  G4double xFit;
  //  G4double y1Fit;
  //  G4double y2Fit;

};

#endif
