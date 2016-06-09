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
//      File name:     G4XnpTotal
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4XNPTOTAL_HH
#define G4XNPTOTAL_HH

#include "globals.hh"
#include "G4CrossSectionPatch.hh"
#include "G4CrossSectionVector.hh"

class G4KineticTrack;

class G4XnpTotal : public G4CrossSectionPatch
{

public:

  G4XnpTotal();

  virtual ~G4XnpTotal();

  G4bool operator==(const G4XnpTotal &right) const;
  G4bool operator!=(const G4XnpTotal &right) const;

  virtual const G4CrossSectionVector* GetComponents() const { return components; } 
 
  virtual G4String Name() const;


protected:


private:  

  G4XnpTotal(const G4XnpTotal &right);
  const G4XnpTotal& operator=(const G4XnpTotal &right);

  G4CrossSectionVector* components;

};

#endif


















