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
//
/// \file radiobiology/src/IonLet.cc
/// \brief Implementation of the RadioBio::IonLet class

#include "IonLet.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonLet::IonLet(G4int trackID, G4int PDG, G4String fullname, G4String name, G4int Z, G4int A,
               G4int voxNumber)
  : fIsPrimary(trackID == 1), fPDGencoding(PDG), fFullName(fullname), fName(name), fZ(Z), fA(A)
{
  fLETDN = array_type(0.0, voxNumber);
  fLETDD = array_type(0.0, voxNumber);
  fLETTN = array_type(0.0, voxNumber);
  fLETTD = array_type(0.0, voxNumber);

  fLETD = array_type(0.0, voxNumber);
  fLETT = array_type(0.0, voxNumber);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IonLet::Update(G4int voxel, G4double DE, G4double DEELETrons, G4double Lsn, G4double DX)
{
  // ions dose LET Numerator, including secondary electrons energy deposit
  fLETDN[voxel] += (DE + DEELETrons) * Lsn;
  // ions dose LET Denominator, including secondary electrons energy deposit
  fLETDD[voxel] += DE + DEELETrons;
  // ions track LET Numerator
  fLETTN[voxel] += DX * Lsn;
  // ions track LET Denominator
  fLETTD[voxel] += DX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IonLet::Merge(const IonLet* rhs)
{
  // If programmed correctly, this exception should never appear
  if (rhs->GetPDGencoding() != fPDGencoding || rhs->IsPrimary() != fIsPrimary)
    G4Exception("IonLet::merge", "mergingdifferentions", FatalException,
                "Cannotmerge ions, probably merging of data from different ions");

  fLETDN += rhs->GetLETDN();
  fLETDD += rhs->GetLETDD();
  fLETTN += rhs->GetLETTN();
  fLETTD += rhs->GetLETTD();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IonLet::Calculate()
{
  for (unsigned int v = 0; v < fLETD.size(); v++) {
    if (fLETDD[v] > 0.) fLETD[v] = fLETDN[v] / fLETDD[v];
    if (fLETTD[v] > 0.) fLETT[v] = fLETTN[v] / fLETTD[v];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio