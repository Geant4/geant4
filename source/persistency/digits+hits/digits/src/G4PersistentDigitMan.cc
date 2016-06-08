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
// $Id: G4PersistentDigitMan.cc,v 1.6 2001/07/11 10:02:14 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// class G4PersistentDigitMan 
//
// History:
// 99.11.22 Y.Morita  Initial version

#include "G4PersistentDigitMan.hh"


G4PersistentDigitMan* G4PersistentDigitMan::f_PersistentDigitMan = 0;

G4PersistentDigitMan* G4PersistentDigitMan::GetPersistentDigitMan()
{
  if(!f_PersistentDigitMan)
  {
    f_PersistentDigitMan = new G4PersistentDigitMan;
  }
  return f_PersistentDigitMan;
}

G4PersistentDigitMan* G4PersistentDigitMan::get_PersistentDigitManIfExist()
{ return f_PersistentDigitMan; }

G4PersistentDigitMan::G4PersistentDigitMan()
 : f_CurrentPDCofThisEvent(0)
{;}

G4PersistentDigitMan::~G4PersistentDigitMan()
{;}

