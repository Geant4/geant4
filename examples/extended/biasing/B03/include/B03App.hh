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
// $Id: B03App.hh,v 1.3 2002-12-13 11:26:33 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
#ifndef B03App_hh
#define B03App_hh B03App_hh

class G4RunManager;
class B03PrimaryGeneratorAction;
class B03PhysicsList;
class B03DetectorConstruction;

class B03AppBase {
public:
  ~B03AppBase();
  static B03AppBase &GetB03AppBase();
  G4RunManager &GetRunManager(){return *frunMgr;}
  
private:
  B03AppBase();
  static B03AppBase *fB03AppBase;
  G4RunManager *frunMgr;

  B03DetectorConstruction *fDetector;
  B03PrimaryGeneratorAction *fPrimary;
  B03PhysicsList *fPhysics;
};

// B03AppBase *GetB03AppBase();



#endif










