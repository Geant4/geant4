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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "globals.hh"
#include "Randomize.hh"
#include <iostream>

main()
{
  for(int i=0; i<10000; i++)
  {
    double x1 = RandGauss::shoot();
    double x2 = RandGauss::shoot();
    double x3 = RandGauss::shoot();
    double x4 = RandGauss::shoot();
    double x5 = RandGauss::shoot();
    double x6 = RandGauss::shoot();
    double pt_1 = x1;
    double pt_2 = x1+x2;
    double pt_3 = x1+x2+x3;
    double pt_4 = x1+x2+x3+x4;
    double pt_5 = x1+x2+x3+x4+x5;
    double pt_6 = x1+x2+x3+x4+x5+x6;
    cout <<pt_1<<" "<<pt_2<<" "<<pt_3<<" "<<pt_4<<" "<<pt_5<<" "<<pt_6<<" "<<G4endl;
  }
}
