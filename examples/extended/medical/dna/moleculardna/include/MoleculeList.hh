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
/// file:  MolecularMoleculeList.hh
/// brief: Enumerators for DNA molecules

#ifndef MOLECULAR_MOLECULE_LIST_HH
#define MOLECULAR_MOLECULE_LIST_HH

#include "globals.hh"
#include "DNAHashing.hh"

typedef enum molecule
{
  UNSPECIFIED,
  SUGAR,
  PHOSPHATE,
  GUANINE,
  ADENINE,
  CYTOSINE,
  THYMINE,
  Count
} molecule;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace utility
{
  using namespace G4::hashing;
  // Larson: works for the given set of words and is the fastest
#define CT_HASHER(x) larson::Cthash(x)
#define DEFINE_HASHER
#define RT_HASHER(x) larson::Hash(x)

  // Default
  //#define CT_HASHER(x) "x"_hash
  //#define DEFINE_HASHER G4ThreadLocalStatic G4::hashing::hasher<std::string>
  //_hasher; #define RT_HASHER(x) _hasher(x)

  // CRC32
  //#define CT_HASHER(x) crc32::hash(x)
  //#define DEFINE_HASHER
  //#define RT_HASHER(x) crc32::hash(x)

  // STL
  //  static constexpr size_t sugar_h = 18188831749337764501llu;
  //  static constexpr size_t p_h = 18387783588576939323llu;
  //  static constexpr size_t a_h = 13783237927007415;
  //  static constexpr size_t g_h = 8059069078009542073;
  //  static constexpr size_t t_h = 11553912749450711402llu;
  //  static constexpr size_t c_h = 17880634351340259280llu;
  //#define CT_HASHER(x) crc32::hash(x)
  //#define DEFINE_HASHER G4ThreadLocalStatic std::hash<std::string> stdhasher;
  //#define RT_HASHER(x) stdhasher(x)

  static constexpr size_t fSugar_h = CT_HASHER("SUGAR");
  static constexpr size_t fP_h     = CT_HASHER("PHOSPHATE");
  static constexpr size_t fA_h     = CT_HASHER("ADENINE");
  static constexpr size_t fG_h     = CT_HASHER("GUANINE");
  static constexpr size_t fT_h     = CT_HASHER("THYMINE");
  static constexpr size_t fC_h     = CT_HASHER("CYTOSINE");

  inline ::molecule GetMoleculeEnum(const G4String& mol)
  {
    G4String _mol(mol);
    std::transform(_mol.begin(), _mol.end(), _mol.begin(), ::toupper);

    DEFINE_HASHER
    size_t mol_h = RT_HASHER(std::string(_mol));

    switch(mol_h)
    {
      case fSugar_h:
        //          G4cout << "RETURN SUGAR FOR : " << mol << G4endl;
        return SUGAR;
      case fP_h:
        //          G4cout << "RETURN P FOR : " << mol << G4endl;
        return PHOSPHATE;
      case fA_h:
        //          G4cout << "RETURN A FOR : " << mol << G4endl;
        return ADENINE;
      case fG_h:
        //          G4cout << "RETURN G FOR : " << mol << G4endl;
        return GUANINE;
      case fC_h:
        //          G4cout << "RETURN C FOR : " << mol << G4endl;
        return CYTOSINE;
      case fT_h:
        //          G4cout << "RETURN T FOR : " << mol << G4endl;
        return THYMINE;
      default:
        //          G4cout << "RETURN UNSPECIFIED FOR : " << mol << G4endl;
        //          G4cout << "p_h " << p_h << G4endl;
        G4Exception("MolecularMoleculeList::GetMoleculeEnum",
                    "ERR_UNKNOWN_MOLECULE", JustWarning, "Unknown molecule");
        return UNSPECIFIED;
    }
  }

  inline G4String GetMoleculeEnumString(::molecule mol)
  {
    switch(mol)
    {
      case PHOSPHATE:
        return "Phosphate";
      case SUGAR:
        return "Sugar";
      case GUANINE:
        return "Guanine";
      case ADENINE:
        return "Adenine";
      case CYTOSINE:
        return "Cytosine";
      case THYMINE:
        return "Thymine";
      default:
        return "Unspecified";
    }
  }

  // TODO : move to test part
  inline int TestMoleculeEnum()
  {
    G4cout << GetMoleculeEnumString(
                utility::GetMoleculeEnum(G4String("PHOSPHATE")))
           << G4endl;

    G4cout << GetMoleculeEnumString(utility::GetMoleculeEnum(G4String("SUGAR")))
           << G4endl;

    G4cout << GetMoleculeEnumString(
                utility::GetMoleculeEnum(G4String("GUANINE")))
           << G4endl;

    G4cout << GetMoleculeEnumString(
                utility::GetMoleculeEnum(G4String("ADENINE")))
           << G4endl;

    G4cout << GetMoleculeEnumString(
                utility::GetMoleculeEnum(G4String("THYMINE")))
           << G4endl;

    G4cout << GetMoleculeEnumString(
                utility::GetMoleculeEnum(G4String("CYTOSINE")))
           << G4endl;

    G4cout << GetMoleculeEnumString(
                utility::GetMoleculeEnum(G4String("CyToSiNe")))
           << G4endl;

    return 0;
  }
}  // namespace utility

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_MOLECULE_LIST_HH
