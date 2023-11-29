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
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia

#include "G4VisAttributes.hh"
#include "G4HumanPhantomColour.hh"
#include "G4Colour.hh"

G4HumanPhantomColour::G4HumanPhantomColour()
{
  fWhite = G4Colour(1.0, 1.0, 1.0);
  fPink = G4Colour(0.94, 0.5, 0.5);
  fGrey = G4Colour(0.46, 0.53, 0.6);
  fYellow = G4Colour(1.0, 1.0, 0.);
  fBlue = G4Colour(0.25,0.41, 0.88 );
  fLightBlue = G4Colour(0.28, 0.82, 0.8);
  fGreen = G4Colour(0., 1., 0.);
  fBrown = G4Colour(0.5, 0.5, 0.);
  fPurple = G4Colour(0.85,0.44,0.84);
  fRed = G4Colour(1.0, 0.0, 0.0);
  fOrange = G4Colour(1.,0.5,0.); 
  fBlack =  G4Colour(0.,0.,0.); 
}

G4Colour G4HumanPhantomColour::GetColour(const G4String& colourName)
{
 // Returns the colour
  if (colourName == "pink") return fPink;
  else if(colourName == "white") return fWhite;
  else if (colourName == "grey") return fGrey;
  else if (colourName == "yellow") return fYellow;
  else if (colourName == "blue") return fBlue;
  else if (colourName == "lightBlue") return fLightBlue;
  else if (colourName == "green") return fGreen;
  else if (colourName == "brown") return fBrown;
  else if (colourName == "purple") return fPurple;
  else if (colourName == "red") return fRed;
  else if (colourName == "orange") return fOrange;
  else if  (colourName == "black") return fBlack; 
else {G4cout<< colourName << "does not exist !!!"<< G4endl; return fWhite;}
}
