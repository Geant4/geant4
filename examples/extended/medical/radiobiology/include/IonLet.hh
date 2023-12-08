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
/// \file radiobiology/include/IonLet.hh
/// \brief Definition of the RadioBio::IonLet class

#ifndef RadiobiologyIonLet_HH
#define RadiobiologyIonLet_HH

#include "globals.hh"

#include <valarray>

namespace RadioBio
{

/// class to save and hold data for LET of different ions
class IonLet
{
  public:
    // Constructor wants ion data, trackID and total voxel number
    // trackID used only to see if particle is primary
    IonLet(G4int trackID, G4int PDG, G4String fullname, G4String name, G4int Z, G4int A,
           G4int voxNumber);
    ~IonLet() = default;

    // Alias for matrix type
    using array_type = std::valarray<G4double>;

    G4bool IsPrimary() const { return fIsPrimary; }
    G4int GetPDGencoding() const { return fPDGencoding; }
    G4String GetFullName() const { return fFullName; }
    G4String GetName() const { return fName; }
    G4int GetZ() const { return fZ; }
    G4int GetA() const { return fA; }

    // Track and Dose LET numerator and denominator
    array_type GetLETDN() const { return fLETDN; }
    array_type GetLETDD() const { return fLETDD; }
    array_type GetLETTN() const { return fLETTN; }
    array_type GetLETTD() const { return fLETTD; }

    // Final Dose and Track LET for the ion
    array_type GetLETD() const { return fLETD; }
    array_type GetLETT() const { return fLETT; }

    void SetLETDN(array_type LETDN) { fLETDN = LETDN; }
    void SetLETDD(array_type LETDD) { fLETDD = LETDD; }
    void SetLETTN(array_type LETTN) { fLETTN = LETTN; }
    void SetLETTD(array_type LETTD) { fLETTD = LETTD; }

    // To update data inside this IonLet
    void Update(G4int voxel, G4double DE, G4double DEELETrons, G4double Lsn, G4double DX);

    // To merge data from another IonLet object inside this one
    void Merge(const IonLet* lhs);

    // To calculate LET given the numerator and denominator
    void Calculate();

    // To sort by the mass number, else sort by the atomic one.
    G4bool operator<(const IonLet& a) const
    {
      return (this->fZ == a.fZ) ? this->fA < a.fA : this->fZ < a.fZ;
    }

  private:
    G4bool fIsPrimary = true;  // True if particle is primary
    G4int fPDGencoding = -1;  // Particle data group id for the particle
    G4String fFullName ;  // AZ[excitation energy]: like He3[1277.4], ...
    G4String fName ;  // simple name no excitation energy: He3, He4, ...
    G4int fZ = -1;  // atomic number
    G4int fA = -1;  // mass number

    // Track averaged LET and Dose averaged LET
    // Numerator, denominator, actual value.
    array_type fLETDN = {};
    array_type fLETDD = {};
    array_type fLETTN = {};
    array_type fLETTD = {};
    array_type fLETD = {};
    array_type fLETT = {};
};

}  // namespace RadioBio

#endif  // IonLet_HH
