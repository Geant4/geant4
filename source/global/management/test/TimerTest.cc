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
// $Id: TimerTest.cc,v 1.1 2003-04-07 13:00:54 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
#include "G4Timer.hh"
#include "globals.hh"
#include "G4ios.hh"

int main()
{
  G4Timer timer;
  G4double j=1.0;

  timer.Start();
  for (size_t i=0; i<100000000; i++) j=(i+1)*(j+1)/(5+j*i);
  G4cout << j << G4endl;
  timer.Stop();  

  G4cout << "System time: " << timer.GetSystemElapsed() << G4endl
         << "User time:   " << timer.GetUserElapsed() << G4endl;

  return 0;
}
