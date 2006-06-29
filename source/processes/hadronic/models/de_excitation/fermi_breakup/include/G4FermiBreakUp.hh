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
// $Id: G4FermiBreakUp.hh,v 1.3 2006-06-29 20:11:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4FermiBreakUp_h
#define G4FermiBreakUp_h 1

#include "G4VFermiBreakUp.hh"
#include "G4FermiConfiguration.hh"
#include "G4FermiConfigurationList.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

class G4FermiBreakUp : public G4VFermiBreakUp 
{
public:
  G4FermiBreakUp();
  ~G4FermiBreakUp();
  
private:
  G4FermiBreakUp(const G4FermiBreakUp &right);
  
  const G4FermiBreakUp & operator=(const G4FermiBreakUp &right);
  G4bool operator==(const G4FermiBreakUp &right) const;
  G4bool operator!=(const G4FermiBreakUp &right) const;
  
public:
  G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);
};


#endif


