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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4StoppingHadronBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard and add mu-
//----------------------------------------------------------------------------
//

#ifndef G4StoppingHadronBuilder_h
#define G4StoppingHadronBuilder_h 1

#include "globals.hh"

class G4MuonMinusCaptureAtRest;
class G4AntiProtonAnnihilationAtRest;
class G4AntiNeutronAnnihilationAtRest;
class G4PionMinusAbsorptionAtRest;
class G4KaonMinusAbsorption;


class G4StoppingHadronBuilder
{
public:
  G4StoppingHadronBuilder();
  virtual ~G4StoppingHadronBuilder();

  virtual void Build();

private:

  G4MuonMinusCaptureAtRest*        theMuonMinusAbsorption;
  G4PionMinusAbsorptionAtRest*     thePionMinusAbsorption;
  G4KaonMinusAbsorption*           theKaonMinusAbsorption;
  G4AntiProtonAnnihilationAtRest*  theAntiProtonAnnihilation;
  G4AntiNeutronAnnihilationAtRest* theAntiNeutronAnnihilation;

  G4bool wasActivated;
};

#endif

