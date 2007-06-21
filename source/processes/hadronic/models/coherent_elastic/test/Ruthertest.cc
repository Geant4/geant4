//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// J.P. Wellisch, X-mas 2002.
#include "G4Rutherford.hh"
#include "G4ParticleTable.hh"

int main()
{
  G4Rutherford theR;
  G4ParticleDefinition * theP 
     = G4ParticleTable::GetParticleTable()->FindIon(4,8,0,4);
  G4Nucleus theN(28, 14);
  int i=0;
  while (i++<1000000) 
  {
    if(i==1000*(i/1000)) G4cerr << "Event # "<<i<<G4endl;
    G4cout << theR.Apply(theP, theN)<<G4endl;
  }
}
