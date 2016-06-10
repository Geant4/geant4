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
// $Id: G4GenericIon.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      4-th Dec 1998, H.Kurashige
// ****************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ****************************************************************
// This class is used only by G4IonTable and not for tracking
// G4IonTable creates various ions other than alpha,deuteron,triton,
// and He3. Processes for these ions will be same as ones for 
// this "GenericIon". So, user should register processes for ions
// to this class in his/her UserPhysicsList
// ----------------------------------------------------------------------

#ifndef G4GenericIon_h
#define G4GenericIon_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4Ions.hh"

// ######################################################################
// ###                          GenericIon                            ###
// ######################################################################

class G4GenericIon : public G4Ions
{
 private:
   static G4GenericIon* theInstance;
   G4GenericIon(){}
   ~G4GenericIon(){}

 public:
   static G4GenericIon* Definition();
   static G4GenericIon* GenericIonDefinition();
   static G4GenericIon* GenericIon();
};

#endif
