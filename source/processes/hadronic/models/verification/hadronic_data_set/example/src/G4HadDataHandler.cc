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

#include "G4HadDataHandler.hh"
#include "G4DataVector.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include <fstream>
#include <strstream>
#include "globals.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>








G4HadDataHandler::G4HadDataHandler()
{
  G4int n = 0;
  G4int pos = 0;

  for (pos = 0 ; pos < 100; pos++) numberOfSecondaries.push_back(n);
      
}

G4HadDataHandler::~G4HadDataHandler()
{ 

}

//////////////////////////////////////////////////////////////////
//
//






//
//
/////////////////////////////////////////////////////////////
