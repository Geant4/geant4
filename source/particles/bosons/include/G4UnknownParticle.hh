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
// $Id: G4UnknownParticle.hh,v 1.2 2004-09-02 01:52:26 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//  first implementation : M.Asai Jul 07, 2004
// ----------------------------------------------------------------
//  New impelemenataion as an utility class  H.Kurashige, 14 July 2004
// ----------------------------------------------------------------


#ifndef G4UnknownParticle_h
#define G4UnknownParticle_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                         UNKNOWN                                ###
// ######################################################################

class G4UnknownParticle
{
 private:
   static G4ParticleDefinition* theInstance;

 private:
  G4UnknownParticle(){}

 public:
   ~G4UnknownParticle(){}
 
   static G4ParticleDefinition* Definition();
   static G4ParticleDefinition* UnknownParticleDefinition();
   static G4ParticleDefinition* UnknownParticle();
};

#endif











