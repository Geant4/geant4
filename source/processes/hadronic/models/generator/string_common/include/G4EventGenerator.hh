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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4EventGenerator.hh,v 1.4 2001/08/01 17:08:57 hpw Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
#ifndef G4EventGenerator_h
#define G4EventGenerator_h 1

#include "G4HadronicInteraction.hh"
#include "G4VertexCode.hh"
#include "G4InteractionCode.hh"

class G4EventGenerator : public G4HadronicInteraction
{
  public:
      G4EventGenerator();
      ~G4EventGenerator();

  private:
      G4EventGenerator(const G4EventGenerator &right);
      const G4EventGenerator & operator=(const G4EventGenerator &right);
      int operator==(const G4EventGenerator &right) const;
      int operator!=(const G4EventGenerator &right) const;

  public:
      virtual G4double GetWidth(G4VertexCode &theCode) = 0;
      virtual G4double GetCrossSection(G4InteractionCode &theCode) = 0;
      virtual G4ReactionProduct * GetFinalState(G4InteractionCode &theCode) = 0;
     // please also define ApplyYourself for scattering off Hydrogen
};


#endif


