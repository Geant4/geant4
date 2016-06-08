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
// $Id: G4PersistentHitMan.cc,v 1.6 2001/07/11 10:02:15 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// class G4PersistentHitMan 
//
// History:
// 99.11.16 Y.Morita  Initial version

#include "G4PersistentHitMan.hh"


G4PersistentHitMan* G4PersistentHitMan::f_PersistentHitMan = 0;

G4PersistentHitMan* G4PersistentHitMan::GetPersistentHitMan()
{
  if(!f_PersistentHitMan)
  {
    f_PersistentHitMan = new G4PersistentHitMan;
  }
  return f_PersistentHitMan;
}

G4PersistentHitMan* G4PersistentHitMan::get_PersistentHitManIfExist()
{ return f_PersistentHitMan; }

G4PersistentHitMan::G4PersistentHitMan()
 : f_CurrentPHCofThisEvent(0)
{;}

G4PersistentHitMan::~G4PersistentHitMan()
{;}

