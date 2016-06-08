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
// $Id: G4XPDGElastic.hh,v 1.2 2002/12/12 19:17:44 gunter Exp $ //
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
#include "g4std/algorithm"
#include "g4std/map"

typedef G4std::pair<G4String,G4String> G4StringPair;

class G4KineticTrack;
class G4ParticleDefinition;

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

  G4std::map <G4StringPair, G4std::vector<G4double>, G4std::less<G4StringPair> > xMap;

  typedef G4std::map <G4StringPair, G4std::vector<G4double>, G4std::less<G4StringPair> > PairDoubleMap;

  //  G4double pMinFit;
  //  G4double pMaxFit;
  //  G4double aFit;
  //  G4double bFit;
  //  G4double nFit;
  //  G4double cFit;
  //  G4double dFit;

};

#endif
