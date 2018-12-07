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
//
//
#ifndef G4MaleBuilder_h
#define G4MaleBuilder_h 1

//#include "G4VPhysicalVolume.hh"

#include "G4PhantomBuilder.hh"

class G4PhantomBuilder;
class G4MaleBuilder: public G4PhantomBuilder
{
public:
  G4MaleBuilder();
  ~G4MaleBuilder();

  void BuildMaleGenitalia(const G4String& colourName, G4bool solidVis, G4bool sensitivity);
  void BuildLeftTeste(const G4String& colourName, G4bool solidVis, G4bool sensitivity);
  void BuildRightTeste(const G4String& colourName, G4bool solidVis, G4bool sensitivity);
  // void SetModel(G4String);
};
#endif
