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
#ifndef DumpFrame_h
#define DumpFrame_h

#include "G4FastVector.hh"
#include "G4ReactionProduct.hh"

class DumpFrames
{
public:
static void DumpFrame(G4FastVector<G4ReactionProduct,128> &vec, G4int vecLen)
{
//  cout << vecLen<<G4endl;
//  for(G4int i=0; i<vecLen; i++)
//  {
//    cout << vec[i]->GetDefinition()->GetPDGEncoding()<<" ";
//    cout << vec[i]->GetPositionInNucleus()<<" ";
//    G4double x,y,z;
//    do
//    {
//      x = G4UniformRand()*10-5;
//      y = G4UniformRand()*10-5;
//      z = G4UniformRand()*10-5;
//    }
//    while(std::sqrt(x*x+y*y+z*z)>5);  /* Loop checking, 30-Oct-2015, G.Folger :)*/
//    cout << x<<" "<<y<<" "<<z;
//    cout << vec[i]->GetMomentum()<<" ";
//    cout << vec[i]->GetTotalEnergy();
//    cout << G4endl;
//  }
}
};

#endif
