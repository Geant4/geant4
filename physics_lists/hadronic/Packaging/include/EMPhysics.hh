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
// $Id: EMPhysics.hh,v 1.2 2005-11-25 15:38:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   EMPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 09.11.2005 G.Folger: standard EM is now seperate
//
//----------------------------------------------------------------------------
//
#ifndef EMPhysics_h
#define EMPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4EmStandardBuilder.hh"
#include "G4EMBuilder.hh"


class EMPhysics : public G4VPhysicsConstructor
{
  public: 
    EMPhysics(const G4String& name ="EM");
    virtual ~EMPhysics();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4EmStandardBuilder theEMStandardPhysics;
    G4EMBuilder theEMPhysics;
};

// 2002 by J.P. Wellisch

#endif





