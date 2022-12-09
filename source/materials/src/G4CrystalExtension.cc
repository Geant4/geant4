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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 21-04-16, created by E.Bagli

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4CrystalExtension.hh"
#include "G4AtomicFormFactor.hh"

G4CrystalExtension::G4CrystalExtension(G4Material* mat, const G4String& name)
  : G4VMaterialExtension(name)
  , fMaterial(mat)
  , theUnitCell(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrystalExtension::~G4CrystalExtension()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4complex G4CrystalExtension::
ComputeStructureFactor(G4double kScatteringVector,
                       G4int h,
                       G4int k,
                       G4int l){
    //SF == Structure Factor
    //AFF == Atomic Form Factor
    //GFS == Geometrical Structure Factor
    G4complex SF = G4complex(0.,0.);

    for(auto & anElement: *(fMaterial->GetElementVector())){
        G4double AFF = G4AtomicFormFactor::GetManager()->Get(kScatteringVector,anElement->GetZ());
        
        G4complex GFS = G4complex(0.,0.);

        for(const auto& anAtomPos : GetAtomBase(anElement)->GetPos())
        {
            G4double aDouble = h * anAtomPos.x()
            + k * anAtomPos.y()
            + l * anAtomPos.z();
            GFS += G4complex(std::cos(CLHEP::twopi * aDouble),
                             std::sin(CLHEP::twopi * aDouble));
        }

        
        SF += G4complex(AFF * GFS.real(),AFF * GFS.imag());
    }
    return SF;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4complex G4CrystalExtension::
ComputeStructureFactorGeometrical(G4int h,
                                  G4int k,
                                  G4int l){
    //GFS == Geometrical Structure Form Factor
    G4complex GFS = G4complex(0.,0.);
    
    for(auto & anElement: *(fMaterial->GetElementVector())){
      for(const auto& anAtomPos : GetAtomBase(anElement)->GetPos())
      {
        G4double aDouble =
          h * anAtomPos.x() + k * anAtomPos.y() + l * anAtomPos.z();
        GFS += G4complex(std::cos(CLHEP::twopi * aDouble),
                         std::sin(CLHEP::twopi * aDouble));
        }
    }
    return GFS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CrystalExtension::SetElReduced(const ReducedElasticity& mat) {
  for (size_t i=0; i<6; ++i) {
    for (size_t j=0; j<6; ++j) {
      fElReduced[i][j] = mat[i][j];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CrystalExtension::SetCpq(G4int p, G4int q, G4double value) {
  if(p > 0 && p < 7 && q > 0 && q < 7)
  {
    fElReduced[p - 1][q - 1] = value;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrystalAtomBase* G4CrystalExtension::GetAtomBase(const G4Element* anElement){
    if((theCrystalAtomBaseMap.count(anElement)<1)){
        G4String astring =  "Atom base for element " + anElement->GetName()
        + " is not registered." ;
        G4Exception ("G4CrystalExtension::GetAtomBase()", "cry001", JustWarning,astring);
        
        AddAtomBase(anElement, new G4CrystalAtomBase());
    }
    return theCrystalAtomBaseMap[anElement];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalExtension::GetAtomPos(const G4Element* anEl, std::vector<G4ThreeVector>& vecout){
    std::vector<G4ThreeVector> pos;
    for(auto & asinglepos: GetAtomBase(anEl)->GetPos()){
        pos.clear();
        theUnitCell->FillAtomicPos(asinglepos,pos);
        vecout.insert(std::end(vecout), std::begin(pos), std::end(pos));
    }
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalExtension::GetAtomPos(std::vector<G4ThreeVector>& vecout){
    std::vector<G4ThreeVector> pos;
    vecout.clear();
    for(auto & anElement: *(fMaterial->GetElementVector())){
        pos.clear();
        GetAtomPos(anElement,pos);
        vecout.insert(std::end(vecout), std::begin(pos), std::end(pos));
    }
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

