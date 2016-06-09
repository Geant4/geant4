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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4GenericIon.hh,v 1.9 2005/01/14 03:49:13 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

// ######################################################################
// ###                          GenericIon                            ###
// ######################################################################

class G4GenericIon : public G4ParticleDefinition
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
