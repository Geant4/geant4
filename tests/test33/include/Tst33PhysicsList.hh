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
// $Id: Tst33PhysicsList.hh,v 1.2 2002-11-20 13:09:16 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class Tst33PhysicsList
//
// Class description:
//
// taken from Tst12PhysicsList

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33PhysicsList_hh
#define Tst33PhysicsList_hh Tst33PhysicsList_hh

#include "G4VUserPhysicsList.hh"
#include "globals.hh"




class Tst33PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst33PhysicsList();
    virtual ~Tst33PhysicsList();

  protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // 
    virtual void SetCuts();
    
  protected:
  // these methods Construct physics processes and register them
    virtual void ConstructGeneral();
    virtual void ConstructEM();
    virtual void ConstructHad();
    virtual void ConstructLeptHad();
 //
    void  ConstructAllBosons();
    void  ConstructAllLeptons();
    void  ConstructAllMesons();
    void  ConstructAllBaryons();
    void  ConstructAllIons();
    void  ConstructAllShortLiveds();

};

#endif
