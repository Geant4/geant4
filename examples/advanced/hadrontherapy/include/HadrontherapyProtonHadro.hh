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
//    **************************************
//    *                                    *
//    *    HadrontherapyProtonHadro.hh        *
//    *                                    *
//    **************************************
//
// $Id: HadrontherapyProtonHadro.hh,v 1.3 2005-05-18 07:53:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it
// 

#ifndef HadrontherapyProtonHadro_h
#define HadrontherapyProtonHadro_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4QGSParticipants.hh"
#include "G4QGSModel.hh"

class HadrontherapyProtonHadro: public G4VPhysicsConstructor 
{
  public:
    HadrontherapyProtonHadro(const G4String& name = "proton-hadronic");
    virtual ~HadrontherapyProtonHadro();

  protected:
    // Construct particle and physics
    void ConstructParticle(){};
    void ConstructProcess();
};
#endif








