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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli,University of Wollongong, Australia
// Contributions by F. Ambroglini INFN Perugia, Italy
// This class manages colours
//
#ifndef G4HumanPhantomColour_H
#define G4HumanPhantomColour_H 1

#include "globals.hh"

class G4Colour;
class G4HumanPhantomColour
{ 
public:
  G4HumanPhantomColour();
  ~G4HumanPhantomColour()=default;

public:
  //void  DefineColour();
  G4Colour GetColour(const G4String&); //returns the colour

private:
  G4Colour fWhite;
  G4Colour fPink;
  G4Colour fGrey;
  G4Colour fYellow;
  G4Colour fBlue;
  G4Colour fLightBlue;
  G4Colour fGreen;
  G4Colour fBrown;
  G4Colour fPurple;
  G4Colour fRed;
  G4Colour fOrange;
  G4Colour fBlack;
};
#endif
