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
#ifndef DumpFrame_h
#define DumpFrame_h

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
//    while(sqrt(x*x+y*y+z*z)>5);
//    cout << x<<" "<<y<<" "<<z;
//    cout << vec[i]->GetMomentum()<<" ";
//    cout << vec[i]->GetTotalEnergy();
//    cout << G4endl;
//  }
}
};

#endif
