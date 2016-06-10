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
// $Id: G4IonProtonCrossSection.hh 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4IonProtonCrossSection
//
//
// Cross-sections for ion proton scattering up to 20 GeV, getting the low
// energy threshold behaviour right.
// H.P. Wellisch (TRIUMF), D. Axen (British Columbia U.). 1996. 
// Published in Phys.Rev.C54:1329-1332,1996 
//
// Original by H.P.Wellisch 28 June 2001
//
// Modifications:
// 30-07-2010 V.Ivanchenko move virtual methods to source,
//                         add constructor and destructor,
//                         make G4ProtonInelasticCrossSection class member,
//                         fix bug in kinematics
// 18.08.2011 V.Ivanchenko migration to new design
//

#ifndef G4IonProtonCrossSection_h
#define G4IonProtonCrossSection_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"

class G4ProtonInelasticCrossSection;

class G4IonProtonCrossSection : public G4VCrossSectionDataSet
{
public:

  G4IonProtonCrossSection();

  virtual ~G4IonProtonCrossSection();

  G4double GetProtonCrossSection(G4double kinEnergy, G4int Z);

  virtual
  G4bool IsElementApplicable(const G4DynamicParticle* aPart, G4int Z,
			     const G4Material*);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle* aPart, 
				  G4int Z, const G4Material*);

  virtual void CrossSectionDescription(std::ostream&) const;

private: // Without Description

  G4IonProtonCrossSection & operator=(const G4IonProtonCrossSection &right);
  G4IonProtonCrossSection(const G4IonProtonCrossSection&);

  G4ProtonInelasticCrossSection* theForward;

};

#endif
