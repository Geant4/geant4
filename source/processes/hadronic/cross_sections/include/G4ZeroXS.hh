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
// $Id: G4ZeroXS.hh 76889 2013-11-18 13:01:55Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ZeroXS
//
// Author:        Tatsumi Koi
//
// Creation date: 26.10.2015
// Modifications:
//
//
// Class Description:
//
// An artificial cross section data set which always replys zero
//
// -------------------------------------------------------------------
//

#ifndef G4ZeroXS_h
#define G4ZeroXS_h 1

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

class G4ZeroXS : public G4VCrossSectionDataSet
{

public:

  G4ZeroXS ();

  virtual ~G4ZeroXS();
   
  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int , 
			     const G4Material* mat = 0);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, G4int ,
				  const G4Material* mat = 0);

  virtual void CrossSectionDescription(std::ostream&) const;

private:

  G4ZeroXS & operator=(const G4ZeroXS &right);
  G4ZeroXS(const G4ZeroXS&);

};

#endif
