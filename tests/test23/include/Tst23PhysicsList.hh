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
// $Id: Tst23PhysicsList.hh,v 1.3 2004-03-05 15:23:18 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is an derived class of G4VPhysicsConstructor
//
// ------------------------------------------- 
//	History
//        first version                   Dec. 2003 by M.Kossov 
// ------------------------------------------------------------
#ifndef Tst23PhysicsList_h
#define Tst23PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class Tst23PhysicsList: public G4VModularPhysicsList
{
public:
  Tst23PhysicsList();
  virtual ~Tst23PhysicsList();
  
public:
  // SetCuts() 
  virtual void SetCuts();


};


#endif



