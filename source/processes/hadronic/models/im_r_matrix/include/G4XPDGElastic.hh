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
// $Id: G4XPDGElastic.hh,v 1.3 2010-03-12 15:45:18 gunter Exp $ //
// ---------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4XPDGElastic
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4XPDGELASTIC_HH
#define G4XPDGELASTIC_HH

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include <algorithm>
#include <map>

class G4KineticTrack;
class G4ParticleDefinition;

typedef std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> G4pDefPair;

class G4XPDGElastic : public G4VCrossSectionSource
{

public:

  G4XPDGElastic();

  virtual ~G4XPDGElastic();

  G4bool operator==(const G4XPDGElastic &right) const;
  G4bool operator!=(const G4XPDGElastic &right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
 
  virtual const G4CrossSectionVector* GetComponents() const { return 0; }

  virtual G4String Name() const;

  virtual G4bool IsValid(G4double e) const;

  virtual G4double LowLimit() const { return _lowLimit; }


protected:


private:  

  G4XPDGElastic(const G4XPDGElastic &right);
  const G4XPDGElastic& operator=(const G4XPDGElastic &right);
  
  static const G4double _lowLimit;
  static const G4double _highLimit;

  static const G4int nPar;
  static const G4double pPiPlusPDGFit[7];
  static const G4double pPiMinusPDGFit[7];
  static const G4double pKPlusPDGFit[7];
  static const G4double pKMinusPDGFit[7];
  static const G4double ppPDGFit[7];
  static const G4double ppbarPDGFit[7];
  static const G4double npbarPDGFit[7];

  std::map <G4pDefPair, std::vector<G4double>, std::less<G4pDefPair> > xMap;

  typedef std::map <G4pDefPair, std::vector<G4double>, std::less<G4pDefPair> > PairDoubleMap;

  //  G4double pMinFit;
  //  G4double pMaxFit;
  //  G4double aFit;
  //  G4double bFit;
  //  G4double nFit;
  //  G4double cFit;
  //  G4double dFit;

};

#endif
