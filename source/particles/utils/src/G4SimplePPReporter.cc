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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SimplePPReporter.cc,v 1.1 2003/09/21 19:38:50 kurasige Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 
// ---------------------------------------------------------------
#include "G4SimplePPReporter.hh"
#include "G4ios.hh"
#include "globals.hh"

//////////////////////////////
G4SimplePPReporter::G4SimplePPReporter():G4VParticlePropertyReporter()
{ 
}

////////////////////////////
G4SimplePPReporter::~G4SimplePPReporter()
{
}    

/////////////////////
void G4SimplePPReporter::Print(const G4String& )
{
  for (size_t i=0; i<pList.size(); i++){
    G4ParticlePropertyData*    ptr = (pList)[i];
    ptr->Print();  
    G4cout << G4endl;
  }
}    
