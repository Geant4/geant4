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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4NistMaterialBuilder
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.12.2004
//
// Modifications:
// 31-10-05 Add chemical effect and gas properties (V.Ivanchenko)
// 27.02.06 V.Ivanchenko add ConstructNewGasMaterial
// 11.05.06 V.Ivanchenko add warning flag to FindMaterial method
// 27.06.06 V.Ivanchenko fix graphite description
// 27.07.07 V.Ivanchenko remove dependence on NistManager
// 30.10.09 V.Ivanchenko update density of G4_GRAFITE from PDG'2008
//                       added G4_GRAPHITE_POROUS
// 03.11.09 A.Lechner changed following material names:
//                    From G4_NYLON-6/6 to G4_NYLON-6-6
//                    From G4_NYLON-6/10 to G4_NYLON-6-10
// 12.12.10 A.Ivantchenko added following materials methodes:
//                    BioChemicalMaterials() and SpaceMaterials(),
//                    where new materials are introduced
// 14.06.11 A.Ivantchenko updated body materials (G4_....ICRP)
//                    according ICRU Report 46 (1992) instead of 1975 
//                    data from ICRU Report 37 used previously
// 26.10.11 new scheme for G4Exception  (mma)
// 09.02.12 P.Gumplinger add ConstructNewIdealGasMaterial
// 30.04.13 M.Trocme & S.Seltzer:
//        - Replace AddElementByWeightFraction() by AddElementByAtomCount() 
//          as much as possible
//        - Comment out ill-defined material GLUCOSE 
//        - Fixed density and atom composition of  POLYCHLOROSTYRENE, 
//          POLYVINYL_BUTYRAL, TERPHENYL
// -------------------------------------------------------------------
//
// Class Description:
//
// Element data from the NIST DB on Atomic Weights and Isotope Compositions
// http://physics.nist.gov/PhysRefData/Compositions/index.html
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4NistMaterialBuilder.hh"
#include "G4NistElementBuilder.hh"
#include "G4Element.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4AutoLock.hh"
#include <iomanip>

namespace
{
  G4Mutex nistMaterialMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NistMaterialBuilder::G4NistMaterialBuilder(G4NistElementBuilder* eb, G4int vb)
: elmBuilder(eb),
  verbose(vb),
  nMaterials(0),
  nComponents(0),
  nCurrent(0)
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::FindOrBuildMaterial(const G4String& matname,
                                                       G4bool warning)
{
  if(verbose > 1) {
    G4cout << "G4NistMaterialBuilder::FindOrBuildMaterial " 
	   << matname << G4endl;
  }
  G4String name = matname;
  if("G4_NYLON-6/6" == matname)  { name = "G4_NYLON-6-6"; }
  else if(name == "G4_NYLON-6/10") { name = "G4_NYLON-6-10"; }

  G4Material* mat = FindMaterial(name);
  if(mat != nullptr) { return mat; }

  mat = BuildNistMaterial(name, warning);
  return mat; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::BuildNistMaterial(const G4String& name, 
						     G4bool warning)
{
  G4Material* mat = nullptr;
  // Check if name inside DB
  for (G4int i=0; i<nMaterials; ++i) {

    if (name == names[i]) {
      if(matIndex[i] == -1) { 
	// Build new Nist material 
	mat = BuildMaterial(i); 
      } else { 
	// Nist material was already built
	const G4MaterialTable* theMaterialTable = 
	  G4Material::GetMaterialTable();
	mat = (*theMaterialTable)[matIndex[i]]; 
      }
      return mat;
    }
  }

  if( (verbose == 1 && warning) || verbose > 1) {
    G4cout << "G4NistMaterialBuilder::FindOrBuildMaterial WARNING:"
	   << " material <" << name
	   << "> is not found." << G4endl;
  }
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* 
G4NistMaterialBuilder::FindOrBuildSimpleMaterial(G4int Z, G4bool warn)
{
  G4Material* mat = FindSimpleMaterial(Z);
  if(mat == nullptr) {
    mat = BuildNistMaterial(names[Z], warn); 
  }
  return mat;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::BuildMaterial(G4int i)
{
  G4Material* mat = nullptr;
  if(i >= nMaterials) { return mat; }

  G4AutoLock l(&nistMaterialMutex);
  if(matIndex[i] >= 0) {
    // Nist material was already built
    const G4MaterialTable* theMaterialTable =
      G4Material::GetMaterialTable();
    mat = (*theMaterialTable)[matIndex[i]];

  } else {
    if (verbose > 1) {
      G4cout << "G4NistMaterialBuilder: BuildMaterial #" << i
	     << G4endl;
    }

    G4int nc = components[i];

    // Check gas parameters:
    // defaults may be changed via AddGas() method
    G4double t = NTP_Temperature;
    G4double p = STP_Pressure;
    if(kStateGas == states[i]) {
      size_t nn = idxGas.size();
      if(nn > 0) {
	for(size_t j=0; j<nn; ++j) {
	  if(i == idxGas[j]) {
	    t = gasTemperature[j];
	    p = gasPressure[j];
	    break;
	  }
	}
      }
    }
    mat = new G4Material(names[i],densities[i],nc,states[i],t,p);

    if (verbose>1) { G4cout << "New material nComponents= " << nc << G4endl; }
    if (nc > 0) {
      G4int idx = indexes[i];
      for (G4int j=0; j<nc; ++j) {
	G4int Z = elements[idx+j];
	G4Element* el = elmBuilder->FindOrBuildElement(Z);
	if(el == nullptr)
	  {
	    G4cout << "G4NistMaterialBuilder::BuildMaterial:"
		   << "  ERROR: elements Z= " << Z << " is not found"
		   << " for material " << names[i] << G4endl;
	    G4Exception("G4NistMaterialBuilder::BuildMaterial()", "mat103",
			FatalException, "Failed to construct material");
	    return nullptr;
	  }
	if(atomCount[i]) {
	  mat->AddElement(el,G4lrint(fractions[idx+j]));
	} else {
	  mat->AddElement(el,fractions[idx+j]);
	}
      }
    }

    // Ionisation potential can be defined via NIST DB or
    // Chemical Formula (ICRU37 Report data)
    G4IonisParamMat* ion = mat->GetIonisation();
    G4double exc0 = ion->GetMeanExcitationEnergy();
    G4double exc1 = exc0;
    if(!chFormulas[i].empty()) {
      mat->SetChemicalFormula(chFormulas[i]);
      exc1 = ion->FindMeanExcitationEnergy(mat);
    }
    // If exists, NIST DB data always overwrites other data
    if(ionPotentials[i] > 0.0) { exc1 = ionPotentials[i]; }
    if(exc0 != exc1) { ion->SetMeanExcitationEnergy(exc1); }

    // Index in Material Table
    matIndex[i] = (G4int)mat->GetIndex();
  }
  l.unlock();
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4int>& nbAtoms,
				      G4double dens,
				      G4State state,     
				      G4double temp,  
				      G4double pres)
{
  // Material is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat != nullptr)
  {
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: the material <" << name
	   << "> already exists." << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return mat;
  }

  // Material not in DB
  G4int els = (G4int)elm.size();
  if(els == 0) { 
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: empty list of elements for " << name
	   << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return nullptr;
  } 

  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential is not defined
  G4bool stp = true;
  if(state == kStateGas && 
     (temp != NTP_Temperature || pres != CLHEP::STP_Pressure))
    { stp = false; }

  AddMaterial(name,dens*CLHEP::cm3/CLHEP::g,0,0.,els,state,stp);
  if(!stp) { AddGas(name,temp,pres); }

  for (G4int i=0; i<els; ++i) {
    AddElementByAtomCount(elmBuilder->GetZ(elm[i]), nbAtoms[i]);
  }

  return BuildMaterial(nMaterials-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4double>& w,
				      G4double dens,
				      G4State state,     
				      G4double temp,  
				      G4double pres)
{
  // Material is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat != nullptr)
  {
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: the material <" << name
	   << "> already exists." << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return mat;
  }

  // Material not in DB
  G4int els = (G4int)elm.size();
  if(els == 0) { 
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: empty list of elements for " << name
	   << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return nullptr;
  } 

  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential is not defined
  G4bool stp = true;
  if(state == kStateGas && 
     (temp != NTP_Temperature || pres != CLHEP::STP_Pressure))
    { stp = false; }
  AddMaterial(name,dens*CLHEP::cm3/CLHEP::g,0,0.,els,state,stp);
  if(!stp) { AddGas(name,temp,pres); }

  for (G4int i=0; i<els; ++i) {
    AddElementByWeightFraction(elmBuilder->GetZ(elm[i]), w[i]);
  }
  
  return BuildMaterial(nMaterials-1); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::ConstructNewGasMaterial(
				      const G4String& name,
				      const G4String& nameDB,
				      G4double temp, 
				      G4double pres)
{
  // Material name is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat != nullptr)
  {
    G4cout << "G4NistMaterialBuilder::ConstructNewGasMaterial:"
           << "  WARNING: the material <" << name
	   << "> already exists." << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return mat;
  }

  G4Material* bmat = FindOrBuildMaterial(nameDB);
  if(bmat == nullptr)
  {
    G4cout << "G4NistMaterialBuilder::ConstructNewGasMaterial:"
	   << "  WARNING: the Name <" << nameDB 
	   << "> is NOT in the database: no new gas will be constructed."
	   << G4endl;
    return nullptr;
  }
  if(bmat->GetState() != kStateGas) {
    G4cout << "G4NistMaterialBuilder::ConstructNewGasMaterial:"
	   << "  WARNING:  <" << nameDB 
	   << "> is NOT a gas -  no new gas will be constructed."
	   << G4endl;
    return nullptr;
  }

  G4double dens = bmat->GetDensity()*pres*bmat->GetTemperature()
    /(temp*bmat->GetPressure());
  mat = new G4Material(name,dens,bmat,kStateGas,temp,pres);

  if (verbose>1) {
    G4cout << "G4NistMaterialBuilder::ConstructNewGasMaterial: done" << G4endl;
    G4cout << &mat << G4endl; 
  }	
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::ConstructNewIdealGasMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4int>& nbAtoms,
                                      G4double temp,
                                      G4double pres)
{
  G4State state = kStateGas;

  // Material is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat != nullptr)
  {
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: the material <" << name
           << "> already exists." << G4endl;
    G4cout << "      New material will NOT be built!"
           << G4endl;
    return mat;
  }

  // Material not in DB
  G4int els = (G4int)elm.size();
  if(els == 0) {
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: empty list of elements for " << name
           << G4endl;
    G4cout << "      New material will NOT be built!"
           << G4endl;
    return nullptr;
  }

  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential is not defined
  G4bool stp = true;
  if(temp != NTP_Temperature || pres != CLHEP::STP_Pressure)
    { stp = false; }

  G4double massPerMole = 0.;

  G4int Z = 0;
  for (G4int i=0; i<els; ++i) {
    Z = elmBuilder->GetZ(elm[i]);
    massPerMole += nbAtoms[i]*elmBuilder->GetAtomicMassAmu(Z)*CLHEP::amu_c2;
  }

  G4double dens = massPerMole/(CLHEP::Avogadro*CLHEP::k_Boltzmann*temp/pres);

  if (els == 1) { AddMaterial(name,dens,Z,0.,els,state,stp); }
  else {
    AddMaterial(name,dens,0,0.,els,state,stp);
    for (G4int i=0; i<els; ++i) {
      AddElementByAtomCount(elmBuilder->GetZ(elm[i]), nbAtoms[i]);
    }
  }

  if(!stp) { AddGas(name,temp,pres); }

  return BuildMaterial(nMaterials-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::AddMaterial(const G4String& nameMat, G4double dens,
					G4int Z, G4double pot, 
					G4int ncomp, G4State state, 
					G4bool stp)
{
  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential in eV

  // if ncomp == 1 then Z should be defined and 
  // AddElement should not be applied

  if (nCurrent != 0) {
    G4cout << "G4NistMaterialBuilder::AddMaterial WARNING: previous "
	   << "mixture " << nMaterials << " " << names[nMaterials]
	   << " is not yet complete!"
	   << G4endl;
    G4cout << "         New material " << nameMat << " will not be added."
	   << G4endl;
    return;
  }

  // density in g/cm3, mean ionisation potential in eV

  names.push_back(nameMat);
  chFormulas.emplace_back("");
  densities.push_back(dens*CLHEP::g/CLHEP::cm3);
  ionPotentials.push_back(pot*CLHEP::eV);
  states.push_back(state);
  components.push_back(ncomp);
  indexes.push_back(nComponents);
  STP.push_back(stp);
  matIndex.push_back(-1);
  atomCount.push_back(false);

  if (1 == ncomp && Z > 0) {
    elements.push_back(Z);
    fractions.push_back(1.0);
    atomCount[nMaterials] = true;
    ++nComponents;
    nCurrent = 0;
  } else {
    nCurrent = ncomp;
  }

  ++nMaterials;

  if(verbose > 1) {
    G4cout << "New material " << nameMat << " is prepared; "
           << " nMaterials= " << nMaterials
           << " nComponents= " << nComponents
           << " nCurrent= " << nCurrent
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::SetVerbose(G4int val)
{
  verbose = val;
  elmBuilder->SetVerbose(verbose);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListMaterials(const G4String& mnam) const
{
  if (mnam == "simple")           { ListNistSimpleMaterials(); }
  else if (mnam == "compound")    { ListNistCompoundMaterials(); }
  else if (mnam == "hep")         { ListHepMaterials(); }
  else if (mnam == "space")       { ListSpaceMaterials(); }
  else if (mnam == "bio")         { ListBioChemicalMaterials(); }

  else if (mnam == "all") {
    ListNistSimpleMaterials();
    ListNistCompoundMaterials();
    ListHepMaterials();
    ListSpaceMaterials();
    ListBioChemicalMaterials();

  } else {
    G4cout << "### G4NistMaterialBuilder::ListMaterials: Warning " 
	   << mnam << " list is not known." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListNistSimpleMaterials() const
{
  G4cout << "=======================================================" << G4endl;
  G4cout << "###   Simple Materials from the NIST Data Base      ###" << G4endl;
  G4cout << "=======================================================" << G4endl;
  G4cout << " Z   Name   density(g/cm^3)  I(eV)                     " << G4endl;
  G4cout << "=======================================================" << G4endl;
  for (G4int i=1; i<nElementary; ++i) {DumpElm(i);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListNistCompoundMaterials() const
{
  G4cout << "=============================================================" << G4endl;
  G4cout << "###    Compound Materials from the NIST Data Base          ##" << G4endl;
  G4cout << "=============================================================" << G4endl;
  G4cout << " Ncomp             Name      density(g/cm^3)  I(eV) ChFormula" << G4endl;
  G4cout << "=============================================================" << G4endl;
  for (G4int i=nElementary; i<nNIST; ++i) {DumpMix(i);}
  DumpMix(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListHepMaterials() const
{
  G4cout << "=============================================================" << G4endl;
  G4cout << "###           HEP & Nuclear Materials                      ##" << G4endl;
  G4cout << "=============================================================" << G4endl;
  G4cout << " Ncomp             Name      density(g/cm^3)  I(eV) ChFormula" << G4endl;
  G4cout << "=============================================================" << G4endl;
  for (G4int i=nNIST; i<nHEP; ++i) {DumpMix(i);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListSpaceMaterials() const
{
  G4cout << "=============================================================" << G4endl;
  G4cout << "###           Space ISS Materials                          ##" << G4endl;
  G4cout << "=============================================================" << G4endl;
  G4cout << " Ncomp             Name      density(g/cm^3)  I(eV) ChFormula" << G4endl;
  G4cout << "=============================================================" << G4endl;
  for (G4int i=nHEP; i<nSpace; ++i) {DumpMix(i);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListBioChemicalMaterials() const
{
  G4cout << "=============================================================" << G4endl;
  G4cout << "###          Bio-Chemical Materials                        ##" << G4endl;
  G4cout << "=============================================================" << G4endl;
  G4cout << " Ncomp             Name      density(g/cm^3)  I(eV) ChFormula" << G4endl;
  G4cout << "=============================================================" << G4endl;
  for (G4int i=nSpace; i<nMaterials; ++i) {DumpMix(i);}
  G4cout << "=============================================================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::DumpElm(G4int i) const
{
  G4cout << std::setw(2)  << i << " " 
         << std::setw(6)  << names[i] 
         << std::setw(14) << densities[i]*cm3/g 
         << std::setw(11) << ionPotentials[i]/eV
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::DumpMix(G4int i) const
{
  G4int nc = components[i];
  G4cout << std::setw(2)  << nc << " "
	 << std::setw(26) << names[i] << " " 
         << std::setw(10) << densities[i]*cm3/g 
	 << std::setw(10) << ionPotentials[i]/eV
	 << "   " << chFormulas[i]
	 << G4endl;
  if (nc > 1) {
    G4int imin = indexes[i];
    G4int imax = imin + nc;
    for (G4int j=imin; j<imax; ++j) {
      G4cout << std::setw(10) << elements[j] << std::setw(14) << fractions[j] 
             << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4NistMaterialBuilder::AddGas(const G4String& nameMat, G4double t, G4double p)
{
  for(G4int i=0; i<nMaterials; ++i) {
    if(nameMat == names[i]) {
      idxGas.push_back(i);
      gasTemperature.push_back(t);
      gasPressure.push_back(p);
      return;
    }
  }
  G4cout << "WARNING: G4NistMaterialBuilder::AddGas problem: there is no "
	 << nameMat << " in the list of materials."
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::AddElementByWeightFraction(G4int Z, G4double w)
{
  elements.push_back(Z);
  fractions.push_back(w);
  --nCurrent;
  ++nComponents;
  if (nCurrent == 0) {
    G4int n = nMaterials - 1;
    G4double sum = 0.0;
    G4int imin = indexes[n];
    G4int imax = imin + components[n];

    if(!atomCount[n]) {
      for(G4int i=imin; i<imax; ++i) {sum += fractions[i];}
      if(sum > 0.0)
      {
        for(G4int i = imin; i < imax; ++i)
        {
          fractions[i] /= sum;
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::AddElementByWeightFraction(const G4String& name,
                                                       G4double w)
{
  G4int Z = elmBuilder->GetZ(name);
  AddElementByWeightFraction(Z, w);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::AddElementByAtomCount(G4int Z, G4int nb)
{
  atomCount[nMaterials-1] = true;
  G4double w = (G4double)nb;
  AddElementByWeightFraction(Z, w);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::AddElementByAtomCount(const G4String& name,
						  G4int nb)
{
  atomCount[nMaterials-1] = true;
  G4int Z = elmBuilder->GetZ(name);
  G4double w = (G4double)nb;
  AddElementByWeightFraction(Z, w);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::Initialise()
{
  if (verbose > 0) {
    G4cout << "### G4NistMaterialBuilder::Initialise()" << G4endl;
  }
  NistSimpleMaterials();
  NistCompoundMaterials();
  NistCompoundMaterials2();
  HepAndNuclearMaterials();
  SpaceMaterials();
  BioChemicalMaterials();

  if (verbose > 1) { ListMaterials("all"); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::NistSimpleMaterials()
{
  // density in g/cm3, mean ionisation potential in eV

  AddMaterial("G4_WATER", 1.0,0, 78., 2);
  AddElementByAtomCount("H" ,  2);
  AddElementByAtomCount("O" ,  1);
  chFormulas[nMaterials-1] = "H_2O";

  AddMaterial("G4_H" ,  8.37480e-5,  1,  19.2, 1, kStateGas);
  AddMaterial("G4_He",  1.66322e-4,  2,  41.8, 1, kStateGas);
  AddMaterial("G4_Li",  0.534     ,  3,  40. );
  AddMaterial("G4_Be",  1.848     ,  4,  63.7);
  AddMaterial("G4_B" ,  2.37      ,  5,  76. );
  AddMaterial("G4_C" ,  2.        ,  6,  78. );
  AddMaterial("G4_N" ,  1.16520e-3,  7,  82. , 1, kStateGas);
  AddMaterial("G4_O" ,  1.33151e-3,  8,  95. , 1, kStateGas);
  AddMaterial("G4_F" ,  1.58029e-3,  9, 115. , 1, kStateGas);
  AddMaterial("G4_Ne",  8.38505e-4, 10, 137. , 1, kStateGas);
  AddMaterial("G4_Na",  0.971     , 11, 149. );
  AddMaterial("G4_Mg",  1.74      , 12, 156. );
  AddMaterial("G4_Al",  2.699     , 13, 166. );
  AddMaterial("G4_Si",  2.33      , 14, 173. );
  AddMaterial("G4_P" ,  2.2       , 15, 173. );
  AddMaterial("G4_S" ,  2.0       , 16, 180. );
  AddMaterial("G4_Cl",  2.99473e-3, 17, 174. , 1, kStateGas);
  AddMaterial("G4_Ar",  1.66201e-3, 18, 188.0, 1, kStateGas);
  AddMaterial("G4_K" ,  0.862     , 19, 190. );
  AddMaterial("G4_Ca",  1.55      , 20, 191. );
  AddMaterial("G4_Sc",  2.989     , 21, 216. );
  AddMaterial("G4_Ti",  4.54      , 22, 233. );
  AddMaterial("G4_V" ,  6.11      , 23, 245. );
  AddMaterial("G4_Cr",  7.18      , 24, 257. );
  AddMaterial("G4_Mn",  7.44      , 25, 272. );
  AddMaterial("G4_Fe",  7.874     , 26, 286. );
  AddMaterial("G4_Co",  8.9       , 27, 297. );
  AddMaterial("G4_Ni",  8.902     , 28, 311. );
  AddMaterial("G4_Cu",  8.96      , 29, 322. );
  AddMaterial("G4_Zn",  7.133     , 30, 330. );
  AddMaterial("G4_Ga",  5.904     , 31, 334. );
  AddMaterial("G4_Ge",  5.323     , 32, 350. );
  AddMaterial("G4_As",  5.73      , 33, 347. );
  AddMaterial("G4_Se",  4.5       , 34, 348. );
  AddMaterial("G4_Br",  7.07210e-3, 35, 343. , 1, kStateGas);
  AddMaterial("G4_Kr",  3.47832e-3, 36, 352. , 1, kStateGas);
  AddMaterial("G4_Rb",  1.532     , 37, 363. );
  AddMaterial("G4_Sr",  2.54      , 38, 366. );
  AddMaterial("G4_Y" ,  4.469     , 39, 379. );
  AddMaterial("G4_Zr",  6.506     , 40, 393. );
  AddMaterial("G4_Nb",  8.57      , 41, 417. );
  AddMaterial("G4_Mo", 10.22      , 42, 424. );
  AddMaterial("G4_Tc", 11.50      , 43, 428. );
  AddMaterial("G4_Ru", 12.41      , 44, 441. );
  AddMaterial("G4_Rh", 12.41      , 45, 449. );
  AddMaterial("G4_Pd", 12.02      , 46, 470. );
  AddMaterial("G4_Ag", 10.5       , 47, 470. );
  AddMaterial("G4_Cd",  8.65      , 48, 469. );
  AddMaterial("G4_In",  7.31      , 49, 488. );
  AddMaterial("G4_Sn",  7.31      , 50, 488. );
  AddMaterial("G4_Sb",  6.691     , 51, 487. );
  AddMaterial("G4_Te",  6.24      , 52, 485. );
  AddMaterial("G4_I" ,  4.93      , 53, 491. );
  AddMaterial("G4_Xe",  5.48536e-3, 54, 482. , 1, kStateGas);
  AddMaterial("G4_Cs",  1.873     , 55, 488. );
  AddMaterial("G4_Ba",  3.5       , 56, 491. );
  AddMaterial("G4_La",  6.154     , 57, 501. );
  AddMaterial("G4_Ce",  6.657     , 58, 523. );
  AddMaterial("G4_Pr",  6.71      , 59, 535. );
  AddMaterial("G4_Nd",  6.9       , 60, 546. );
  AddMaterial("G4_Pm",  7.22      , 61, 560. );
  AddMaterial("G4_Sm",  7.46      , 62, 574. );
  AddMaterial("G4_Eu",  5.243     , 63, 580. );
  AddMaterial("G4_Gd",  7.9004    , 64, 591. );
  AddMaterial("G4_Tb",  8.229     , 65, 614. );
  AddMaterial("G4_Dy",  8.55      , 66, 628. );
  AddMaterial("G4_Ho",  8.795     , 67, 650. );
  AddMaterial("G4_Er",  9.066     , 68, 658. );
  AddMaterial("G4_Tm",  9.321     , 69, 674. );
  AddMaterial("G4_Yb",  6.73      , 70, 684. );
  AddMaterial("G4_Lu",  9.84      , 71, 694. );
  AddMaterial("G4_Hf", 13.31      , 72, 705. );
  AddMaterial("G4_Ta", 16.654     , 73, 718. );
  AddMaterial("G4_W" , 19.30      , 74, 727. );
  AddMaterial("G4_Re", 21.02      , 75, 736. );
  AddMaterial("G4_Os", 22.57      , 76, 746. );
  AddMaterial("G4_Ir", 22.42      , 77, 757. );
  AddMaterial("G4_Pt", 21.45      , 78, 790. );
  AddMaterial("G4_Au", 19.32      , 79, 790. );
  AddMaterial("G4_Hg", 13.546     , 80, 800. );
  AddMaterial("G4_Tl", 11.72      , 81, 810. );
  AddMaterial("G4_Pb", 11.35      , 82, 823. );
  AddMaterial("G4_Bi",  9.747     , 83, 823. );
  AddMaterial("G4_Po",  9.32      , 84, 830. );
  AddMaterial("G4_At",  9.32      , 85, 825. );
  AddMaterial("G4_Rn",  9.00662e-3, 86, 794. , 1, kStateGas);
  AddMaterial("G4_Fr",  1.00      , 87, 827. );
  AddMaterial("G4_Ra",  5.00      , 88, 826. );
  AddMaterial("G4_Ac", 10.07      , 89, 841. );
  AddMaterial("G4_Th", 11.72      , 90, 847. );
  AddMaterial("G4_Pa", 15.37      , 91, 878. );
  AddMaterial("G4_U" , 18.95      , 92, 890. );
  AddMaterial("G4_Np", 20.25      , 93, 902. );
  AddMaterial("G4_Pu", 19.84      , 94, 921. );
  AddMaterial("G4_Am", 13.67      , 95, 934. );
  AddMaterial("G4_Cm", 13.51      , 96, 939. );
  AddMaterial("G4_Bk", 14.00      , 97, 952. );
  AddMaterial("G4_Cf", 10.00      , 98, 966. );

  nElementary = nMaterials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::NistCompoundMaterials()
{
  AddMaterial("G4_A-150_TISSUE", 1.127, 0, 65.1, 6);
  AddElementByWeightFraction( 1, 0.101327);
  AddElementByWeightFraction( 6, 0.775501);
  AddElementByWeightFraction( 7, 0.035057);
  AddElementByWeightFraction( 8, 0.052316);
  AddElementByWeightFraction( 9, 0.017422);
  AddElementByWeightFraction(20, 0.018378);

  AddMaterial("G4_ACETONE", 0.7899, 0, 64.2, 3);
  AddElementByAtomCount("C" ,  3);
  AddElementByAtomCount("H" ,  6);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_ACETYLENE", 0.0010967, 0, 58.2, 2, kStateGas);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  2);

  //Tan, Z., et al. NIMB,2006(248)
  AddMaterial("G4_ADENINE", 1.35, 0, 71.4, 3);
  AddElementByAtomCount("C" ,  5);
  AddElementByAtomCount("H" ,  5);
  AddElementByAtomCount("N" ,  5);

  AddMaterial("G4_ADIPOSE_TISSUE_ICRP", 0.95, 0, 63.2, 7);
  AddElementByWeightFraction( 1, 0.114);
  AddElementByWeightFraction( 6, 0.598);
  AddElementByWeightFraction( 7, 0.007);
  AddElementByWeightFraction( 8, 0.278);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(16, 0.001);
  AddElementByWeightFraction(17, 0.001);

  AddMaterial("G4_AIR", 0.00120479, 0, 85.7, 4, kStateGas);
  AddElementByWeightFraction( 6, 0.000124);
  AddElementByWeightFraction( 7, 0.755267);
  AddElementByWeightFraction( 8, 0.231781);
  AddElementByWeightFraction(18, 0.012827);

  AddMaterial("G4_ALANINE", 1.42, 0, 71.9, 4);
  AddElementByAtomCount("C" ,  3);
  AddElementByAtomCount("H" ,  7);
  AddElementByAtomCount("N" ,  1);
  AddElementByAtomCount("O" ,  2);

  AddMaterial("G4_ALUMINUM_OXIDE", 3.97, 0, 145.2, 2);
  AddElementByAtomCount("Al",  2);
  AddElementByAtomCount("O" ,  3);
  chFormulas[nMaterials-1] = "Al_2O_3";

  AddMaterial("G4_AMBER", 1.1, 0, 63.2, 3);
  AddElementByWeightFraction( 1, 0.10593 );
  AddElementByWeightFraction( 6, 0.788973);
  AddElementByWeightFraction( 8, 0.105096);

  AddMaterial("G4_AMMONIA", 0.000826019, 0, 53.7, 2, kStateGas);
  AddElementByAtomCount("N" ,  1);
  AddElementByAtomCount("H" ,  3);

  AddMaterial("G4_ANILINE", 1.0235, 0, 66.2, 3);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" ,  7);
  AddElementByAtomCount("N" ,  1);

  AddMaterial("G4_ANTHRACENE", 1.283, 0, 69.5, 2);
  AddElementByAtomCount("C" , 14);
  AddElementByAtomCount("H" , 10);

  AddMaterial("G4_B-100_BONE", 1.45, 0, 85.9, 6);
  AddElementByWeightFraction( 1, 0.065471);
  AddElementByWeightFraction( 6, 0.536945);
  AddElementByWeightFraction( 7, 0.0215  );
  AddElementByWeightFraction( 8, 0.032085);
  AddElementByWeightFraction( 9, 0.167411);
  AddElementByWeightFraction(20, 0.176589);

  AddMaterial("G4_BAKELITE", 1.25, 0, 72.4, 3);
  AddElementByWeightFraction( 1, 0.057441);
  AddElementByWeightFraction( 6, 0.774591);
  AddElementByWeightFraction( 8, 0.167968);

  AddMaterial("G4_BARIUM_FLUORIDE", 4.89 ,0, 375.9, 2);
  AddElementByAtomCount("Ba",  1);
  AddElementByAtomCount("F" ,  2);

  AddMaterial("G4_BARIUM_SULFATE", 4.5, 0, 285.7, 3);
  AddElementByAtomCount("Ba",  1);
  AddElementByAtomCount("S" ,  1);
  AddElementByAtomCount("O" ,  4);

  AddMaterial("G4_BENZENE", 0.87865, 0, 63.4, 2);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" ,  6);

  AddMaterial("G4_BERYLLIUM_OXIDE", 3.01, 0, 93.2, 2);
  AddElementByAtomCount("Be",  1);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_BGO", 7.13, 0, 534.1, 3);
  AddElementByAtomCount("Bi",  4);
  AddElementByAtomCount("Ge",  3);
  AddElementByAtomCount("O" , 12);

  AddMaterial("G4_BLOOD_ICRP", 1.06, 0, 75.2, 10);
  AddElementByWeightFraction( 1, 0.102);
  AddElementByWeightFraction( 6, 0.110);
  AddElementByWeightFraction( 7, 0.033);
  AddElementByWeightFraction( 8, 0.745);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.003);
  AddElementByWeightFraction(19, 0.002);
  AddElementByWeightFraction(26, 0.001);

  AddMaterial("G4_BONE_COMPACT_ICRU", 1.85, 0, 91.9, 8);
  AddElementByWeightFraction( 1, 0.064);
  AddElementByWeightFraction( 6, 0.278);
  AddElementByWeightFraction( 7, 0.027);
  AddElementByWeightFraction( 8, 0.410);
  AddElementByWeightFraction(12, 0.002);
  AddElementByWeightFraction(15, 0.07 );
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(20, 0.147);

  // Sceleton Cortical bone for Adult ICRU 46
  AddMaterial("G4_BONE_CORTICAL_ICRP", 1.92, 0, 110, 9);
  AddElementByWeightFraction( 1, 0.034);
  AddElementByWeightFraction( 6, 0.155);
  AddElementByWeightFraction( 7, 0.042);
  AddElementByWeightFraction( 8, 0.435);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(12, 0.002);
  AddElementByWeightFraction(15, 0.103);
  AddElementByWeightFraction(16, 0.003);
  AddElementByWeightFraction(20, 0.225);

  AddMaterial("G4_BORON_CARBIDE", 2.52, 0, 84.7, 2);
  AddElementByAtomCount("B" ,  4);
  AddElementByAtomCount("C" ,  1);

  AddMaterial("G4_BORON_OXIDE", 1.812, 0, 99.6, 2);
  AddElementByAtomCount("B" ,  2);
  AddElementByAtomCount("O" ,  3);

  AddMaterial("G4_BRAIN_ICRP", 1.04, 0, 73.3, 9);
  AddElementByWeightFraction( 1, 0.107);
  AddElementByWeightFraction( 6, 0.145);
  AddElementByWeightFraction( 7, 0.022);
  AddElementByWeightFraction( 8, 0.712);
  AddElementByWeightFraction(11, 0.002);
  AddElementByWeightFraction(15, 0.004);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.003);
  AddElementByWeightFraction(19, 0.003);

  AddMaterial("G4_BUTANE", 0.00249343, 0, 48.3, 2, kStateGas);
  AddElementByAtomCount("C" ,  4);
  AddElementByAtomCount("H" , 10);

  AddMaterial("G4_N-BUTYL_ALCOHOL", 0.8098, 0, 59.9, 3);
  AddElementByAtomCount("C" ,  4);
  AddElementByAtomCount("H" , 10);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_C-552", 1.76, 0, 86.8, 5);
  AddElementByWeightFraction( 1, 0.02468 );
  AddElementByWeightFraction( 6, 0.50161 );
  AddElementByWeightFraction( 8, 0.004527);
  AddElementByWeightFraction( 9, 0.465209);
  AddElementByWeightFraction(14, 0.003973);

  AddMaterial("G4_CADMIUM_TELLURIDE", 6.2, 0, 539.3, 2);
  AddElementByAtomCount("Cd",  1);
  AddElementByAtomCount("Te",  1);

  AddMaterial("G4_CADMIUM_TUNGSTATE", 7.9, 0, 468.3, 3);
  AddElementByAtomCount("Cd",  1);
  AddElementByAtomCount("W" ,  1);
  AddElementByAtomCount("O" ,  4);

  AddMaterial("G4_CALCIUM_CARBONATE", 2.8, 0, 136.4, 3);
  AddElementByAtomCount("Ca",  1);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("O" ,  3);

  AddMaterial("G4_CALCIUM_FLUORIDE", 3.18, 0, 166., 2);
  AddElementByAtomCount("Ca",  1);
  AddElementByAtomCount("F" ,  2);

  AddMaterial("G4_CALCIUM_OXIDE", 3.3, 0, 176.1, 2);
  AddElementByAtomCount("Ca",  1);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_CALCIUM_SULFATE", 2.96, 0, 152.3, 3);
  AddElementByAtomCount("Ca",  1);
  AddElementByAtomCount("S" ,  1);
  AddElementByAtomCount("O" ,  4);

  AddMaterial("G4_CALCIUM_TUNGSTATE", 6.062, 0, 395., 3);
  AddElementByAtomCount("Ca",  1);
  AddElementByAtomCount("W" ,  1);
  AddElementByAtomCount("O" ,  4);

  AddMaterial("G4_CARBON_DIOXIDE", 0.00184212, 0, 85., 2, kStateGas);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("O" ,  2);
  chFormulas[nMaterials-1] = "CO_2";

  AddMaterial("G4_CARBON_TETRACHLORIDE", 1.594, 0, 166.3, 2);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("Cl",  4);

  AddMaterial("G4_CELLULOSE_CELLOPHANE", 1.42, 0, 77.6, 3);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" , 10);
  AddElementByAtomCount("O" ,  5);

  AddMaterial("G4_CELLULOSE_BUTYRATE", 1.2, 0, 74.6, 3);
  AddElementByWeightFraction( 1, 0.067125);
  AddElementByWeightFraction( 6, 0.545403);
  AddElementByWeightFraction( 8, 0.387472);

  AddMaterial("G4_CELLULOSE_NITRATE", 1.49, 0, 87., 4);
  AddElementByWeightFraction( 1, 0.029216);
  AddElementByWeightFraction( 6, 0.271296);
  AddElementByWeightFraction( 7, 0.121276);
  AddElementByWeightFraction( 8, 0.578212);

  AddMaterial("G4_CERIC_SULFATE", 1.03, 0, 76.7, 5);
  AddElementByWeightFraction( 1, 0.107596);
  AddElementByWeightFraction( 7, 0.0008  );
  AddElementByWeightFraction( 8, 0.874976);
  AddElementByWeightFraction(16, 0.014627);
  AddElementByWeightFraction(58, 0.002001);

  AddMaterial("G4_CESIUM_FLUORIDE", 4.115, 0, 440.7, 2);
  AddElementByAtomCount("Cs",  1);
  AddElementByAtomCount("F" ,  1);

  AddMaterial("G4_CESIUM_IODIDE", 4.51, 0, 553.1, 2);
  AddElementByAtomCount("Cs",  1);
  AddElementByAtomCount("I" ,  1);

  AddMaterial("G4_CHLOROBENZENE", 1.1058, 0, 89.1, 3);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" ,  5);
  AddElementByAtomCount("Cl",  1);

  AddMaterial("G4_CHLOROFORM", 1.4832, 0, 156., 3);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("H" ,  1);
  AddElementByAtomCount("Cl",  3);

  AddMaterial("G4_CONCRETE", 2.3, 0, 135.2, 10);
  AddElementByWeightFraction( 1, 0.01    );
  AddElementByWeightFraction( 6, 0.001   );
  AddElementByWeightFraction( 8, 0.529107);
  AddElementByWeightFraction(11, 0.016   );
  AddElementByWeightFraction(12, 0.002   );
  AddElementByWeightFraction(13, 0.033872);
  AddElementByWeightFraction(14, 0.337021);
  AddElementByWeightFraction(19, 0.013   );
  AddElementByWeightFraction(20, 0.044   );
  AddElementByWeightFraction(26, 0.014   );

  AddMaterial("G4_CYCLOHEXANE", 0.779, 0, 56.4, 2);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" , 12);

  AddMaterial("G4_1,2-DICHLOROBENZENE", 1.3048, 0, 106.5, 3);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" ,  4);
  AddElementByAtomCount("Cl",  2);

  AddMaterial("G4_DICHLORODIETHYL_ETHER", 1.2199, 0, 103.3, 4);
  AddElementByAtomCount("C" ,  4);
  AddElementByAtomCount("H" ,  8);
  AddElementByAtomCount("O" ,  1);
  AddElementByAtomCount("Cl",  2);

  AddMaterial("G4_1,2-DICHLOROETHANE", 1.2351, 0, 111.9, 3);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  4);
  AddElementByAtomCount("Cl",  2);

  AddMaterial("G4_DIETHYL_ETHER", 0.71378, 0, 60., 3);
  AddElementByAtomCount("C" ,  4);
  AddElementByAtomCount("H" , 10);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_N,N-DIMETHYL_FORMAMIDE", 0.9487, 0, 66.6, 4);
  AddElementByAtomCount("C" ,  3);
  AddElementByAtomCount("H" ,  7);
  AddElementByAtomCount("N" ,  1);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_DIMETHYL_SULFOXIDE", 1.1014, 0, 98.6, 4);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  6);
  AddElementByAtomCount("O" ,  1);
  AddElementByAtomCount("S" ,  1);

  AddMaterial("G4_ETHANE", 0.00125324, 0, 45.4, 2, kStateGas);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  6);

  AddMaterial("G4_ETHYL_ALCOHOL", 0.7893, 0, 62.9, 3);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  6);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_ETHYL_CELLULOSE", 1.13, 0, 69.3, 3);
  AddElementByWeightFraction( 1, 0.090027);
  AddElementByWeightFraction( 6, 0.585182);
  AddElementByWeightFraction( 8, 0.324791);

  AddMaterial("G4_ETHYLENE", 0.00117497, 0, 50.7, 2, kStateGas);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  4);

  AddMaterial("G4_EYE_LENS_ICRP", 1.07, 0, 73.3, 8);
  AddElementByWeightFraction( 1, 0.096);
  AddElementByWeightFraction( 6, 0.195);
  AddElementByWeightFraction( 7, 0.057);
  AddElementByWeightFraction( 8, 0.646);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.003);
  AddElementByWeightFraction(17, 0.001);

  AddMaterial("G4_FERRIC_OXIDE", 5.2, 0, 227.3, 2);
  AddElementByAtomCount("Fe",  2);
  AddElementByAtomCount("O" ,  3);

  AddMaterial("G4_FERROBORIDE", 7.15, 0, 261., 2);
  AddElementByAtomCount("Fe",  1);
  AddElementByAtomCount("B" ,  1);

  AddMaterial("G4_FERROUS_OXIDE", 5.7, 0, 248.6, 2);
  AddElementByAtomCount("Fe",  1);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_FERROUS_SULFATE", 1.024, 0, 76.4, 7);
  AddElementByWeightFraction( 1, 0.108259);
  AddElementByWeightFraction( 7, 2.7e-05 );
  AddElementByWeightFraction( 8, 0.878636);
  AddElementByWeightFraction(11, 2.2e-05 );
  AddElementByWeightFraction(16, 0.012968);
  AddElementByWeightFraction(17, 3.4e-05 );
  AddElementByWeightFraction(26, 5.4e-05 );

  AddMaterial("G4_FREON-12", 1.12, 0, 143., 3);
  AddElementByWeightFraction( 6, 0.099335);
  AddElementByWeightFraction( 9, 0.314247);
  AddElementByWeightFraction(17, 0.586418);

  AddMaterial("G4_FREON-12B2", 1.8, 0, 284.9, 3);
  AddElementByWeightFraction( 6, 0.057245);
  AddElementByWeightFraction( 9, 0.181096);
  AddElementByWeightFraction(35, 0.761659);

  AddMaterial("G4_FREON-13", 0.95, 0, 126.6, 3);
  AddElementByWeightFraction( 6, 0.114983);
  AddElementByWeightFraction( 9, 0.545622);
  AddElementByWeightFraction(17, 0.339396);

  AddMaterial("G4_FREON-13B1", 1.5, 0, 210.5, 3);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("F" ,  3);
  AddElementByAtomCount("Br",  1);

  AddMaterial("G4_FREON-13I1", 1.8, 0, 293.5, 3);
  AddElementByWeightFraction( 6, 0.061309);
  AddElementByWeightFraction( 9, 0.290924);
  AddElementByWeightFraction(53, 0.647767);

  AddMaterial("G4_GADOLINIUM_OXYSULFIDE", 7.44, 0, 493.3, 3);
  AddElementByAtomCount("Gd",  2);
  AddElementByAtomCount("O" ,  2);
  AddElementByAtomCount("S" ,  1);

  AddMaterial("G4_GALLIUM_ARSENIDE", 5.31, 0, 384.9, 2);
  AddElementByAtomCount("Ga",  1);
  AddElementByAtomCount("As",  1);

  AddMaterial("G4_GEL_PHOTO_EMULSION", 1.2914, 0, 74.8, 5);
  AddElementByWeightFraction( 1, 0.08118);
  AddElementByWeightFraction( 6, 0.41606);
  AddElementByWeightFraction( 7, 0.11124);
  AddElementByWeightFraction( 8, 0.38064);
  AddElementByWeightFraction(16, 0.01088);

  AddMaterial("G4_Pyrex_Glass", 2.23, 0, 134., 6);
  AddElementByWeightFraction( 5, 0.040064);
  AddElementByWeightFraction( 8, 0.539562);
  AddElementByWeightFraction(11, 0.028191);
  AddElementByWeightFraction(13, 0.011644);
  AddElementByWeightFraction(14, 0.37722 );
  AddElementByWeightFraction(19, 0.003321);

  AddMaterial("G4_GLASS_LEAD", 6.22, 0, 526.4, 5);
  AddElementByWeightFraction( 8, 0.156453);
  AddElementByWeightFraction(14, 0.080866);
  AddElementByWeightFraction(22, 0.008092);
  AddElementByWeightFraction(33, 0.002651);
  AddElementByWeightFraction(82, 0.751938);

  AddMaterial("G4_GLASS_PLATE", 2.4, 0, 145.4, 4);
  AddElementByWeightFraction( 8, 0.4598  );
  AddElementByWeightFraction(11, 0.096441);
  AddElementByWeightFraction(14, 0.336553);
  AddElementByWeightFraction(20, 0.107205);

  //  AddMaterial("G4_GLUCOSE", 1.54, 0, 77.2, 3);
  //  AddElementByWeightFraction( 1, 0.071204);
  //  AddElementByWeightFraction( 6, 0.363652);
  //  AddElementByWeightFraction( 8, 0.565144);
  //  These weight fractions correspond to 0.90909 glucose + 0.09091 water 
  //  (in mass), not pure glucose. The density used is that of pure glucose!
  //  M.Trocme & S.Seltzer

  AddMaterial("G4_GLUTAMINE", 1.46, 0, 73.3, 4);
  AddElementByAtomCount("C" ,  5);
  AddElementByAtomCount("H" , 10);
  AddElementByAtomCount("N" ,  2);
  AddElementByAtomCount("O" ,  3);

  AddMaterial("G4_GLYCEROL", 1.2613, 0, 72.6, 3);
  AddElementByAtomCount("C" ,  3);
  AddElementByAtomCount("H" ,  8);
  AddElementByAtomCount("O" ,  3);

  //Tan, Z., et al. NIMB,2006(248)
  AddMaterial("G4_GUANINE", 1.58, 0, 75. ,4);
  AddElementByAtomCount("C" ,  5);
  AddElementByAtomCount("H" ,  5);
  AddElementByAtomCount("N" ,  5);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_GYPSUM", 2.32, 0, 129.7, 4);
  AddElementByAtomCount("Ca",  1);
  AddElementByAtomCount("S" ,  1);
  AddElementByAtomCount("O" ,  6);
  AddElementByAtomCount("H" ,  4);

  AddMaterial("G4_N-HEPTANE", 0.68376, 0, 54.4, 2);
  AddElementByAtomCount("C" ,  7);
  AddElementByAtomCount("H" , 16);

  AddMaterial("G4_N-HEXANE", 0.6603, 0, 54., 2);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" , 14);

  AddMaterial("G4_KAPTON", 1.42, 0, 79.6, 4);
  AddElementByAtomCount("C" , 22);
  AddElementByAtomCount("H" , 10);
  AddElementByAtomCount("N" ,  2);
  AddElementByAtomCount("O" ,  5);

  AddMaterial("G4_LANTHANUM_OXYBROMIDE", 6.28, 0, 439.7, 3);
  AddElementByAtomCount("La",  1);
  AddElementByAtomCount("Br",  1);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_LANTHANUM_OXYSULFIDE", 5.86, 0, 421.2, 3);
  AddElementByAtomCount("La",  2);
  AddElementByAtomCount("O" ,  2);
  AddElementByAtomCount("S" ,  1);

  AddMaterial("G4_LEAD_OXIDE", 9.53, 0, 766.7, 2);
  AddElementByWeightFraction( 8, 0.071682);
  AddElementByWeightFraction(82, 0.928318);

  AddMaterial("G4_LITHIUM_AMIDE", 1.178, 0, 55.5, 3);
  AddElementByAtomCount("Li",  1);
  AddElementByAtomCount("N" ,  1);
  AddElementByAtomCount("H" ,  2);

  AddMaterial("G4_LITHIUM_CARBONATE", 2.11, 0, 87.9, 3);
  AddElementByAtomCount("Li",  2);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("O" ,  3);

  AddMaterial("G4_LITHIUM_FLUORIDE", 2.635, 0, 94., 2);
  AddElementByAtomCount("Li",  1);
  AddElementByAtomCount("F" ,  1);

  AddMaterial("G4_LITHIUM_HYDRIDE", 0.82, 0, 36.5, 2);
  AddElementByAtomCount("Li",  1);
  AddElementByAtomCount("H" ,  1);

  AddMaterial("G4_LITHIUM_IODIDE", 3.494, 0, 485.1, 2);
  AddElementByAtomCount("Li",  1);
  AddElementByAtomCount("I" ,  1);

  AddMaterial("G4_LITHIUM_OXIDE", 2.013, 0, 73.6, 2);
  AddElementByAtomCount("Li",  2);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_LITHIUM_TETRABORATE", 2.44, 0, 94.6, 3);
  AddElementByAtomCount("Li",  2);
  AddElementByAtomCount("B" ,  4);
  AddElementByAtomCount("O" ,  7);
}

void G4NistMaterialBuilder::NistCompoundMaterials2()
{
  //Adult Lung congested
  AddMaterial("G4_LUNG_ICRP", 1.04, 0, 75.3, 9);
  AddElementByWeightFraction( 1, 0.105);
  AddElementByWeightFraction( 6, 0.083);
  AddElementByWeightFraction( 7, 0.023);
  AddElementByWeightFraction( 8, 0.779);
  AddElementByWeightFraction(11, 0.002);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.003);
  AddElementByWeightFraction(19, 0.002);

  AddMaterial("G4_M3_WAX", 1.05, 0, 67.9, 5);
  AddElementByWeightFraction( 1, 0.114318);
  AddElementByWeightFraction( 6, 0.655823);
  AddElementByWeightFraction( 8, 0.092183);
  AddElementByWeightFraction(12, 0.134792);
  AddElementByWeightFraction(20, 0.002883);

  AddMaterial("G4_MAGNESIUM_CARBONATE", 2.958, 0, 118., 3);
  AddElementByAtomCount("Mg",  1);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("O" ,  3);

  AddMaterial("G4_MAGNESIUM_FLUORIDE", 3.0, 0, 134.3, 2);
  AddElementByAtomCount("Mg",  1);
  AddElementByAtomCount("F" ,  2);

  AddMaterial("G4_MAGNESIUM_OXIDE", 3.58, 0, 143.8, 2);
  AddElementByAtomCount("Mg",  1);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_MAGNESIUM_TETRABORATE", 2.53, 0, 108.3, 3);
  AddElementByAtomCount("Mg",  1);
  AddElementByAtomCount("B" ,  4);
  AddElementByAtomCount("O" ,  7);

  AddMaterial("G4_MERCURIC_IODIDE", 6.36, 0, 684.5, 2);
  AddElementByAtomCount("Hg",  1);
  AddElementByAtomCount("I" ,  2);

  AddMaterial("G4_METHANE", 0.000667151, 0, 41.7, 2, kStateGas);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("H" ,  4);

  AddMaterial("G4_METHANOL", 0.7914, 0, 67.6, 3);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("H" ,  4);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_MIX_D_WAX", 0.99, 0, 60.9, 5);
  AddElementByWeightFraction( 1, 0.13404 );
  AddElementByWeightFraction( 6, 0.77796 );
  AddElementByWeightFraction( 8, 0.03502 );
  AddElementByWeightFraction(12, 0.038594);
  AddElementByWeightFraction(22, 0.014386);

  AddMaterial("G4_MS20_TISSUE", 1.0, 0, 75.1, 6);
  AddElementByWeightFraction( 1, 0.081192);
  AddElementByWeightFraction( 6, 0.583442);
  AddElementByWeightFraction( 7, 0.017798);
  AddElementByWeightFraction( 8, 0.186381);
  AddElementByWeightFraction(12, 0.130287);
  AddElementByWeightFraction(17, 0.0009  );

  AddMaterial("G4_MUSCLE_SKELETAL_ICRP", 1.05, 0, 75.3, 9);
  AddElementByWeightFraction( 1, 0.102);
  AddElementByWeightFraction( 6, 0.143);
  AddElementByWeightFraction( 7, 0.034);
  AddElementByWeightFraction( 8, 0.710);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.002);
  AddElementByWeightFraction(16, 0.003);
  AddElementByWeightFraction(17, 0.001);
  AddElementByWeightFraction(19, 0.004);

  // from old ICRU report
  AddMaterial("G4_MUSCLE_STRIATED_ICRU", 1.04, 0, 74.7, 8);
  AddElementByWeightFraction( 1, 0.102);
  AddElementByWeightFraction( 6, 0.123);
  AddElementByWeightFraction( 7, 0.035);
  AddElementByWeightFraction( 8, 0.729);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.002);
  AddElementByWeightFraction(16, 0.004);
  AddElementByWeightFraction(19, 0.003);

  AddMaterial("G4_MUSCLE_WITH_SUCROSE", 1.11, 0, 74.3, 4);
  AddElementByWeightFraction( 1, 0.098234);
  AddElementByWeightFraction( 6, 0.156214);
  AddElementByWeightFraction( 7, 0.035451);
  AddElementByWeightFraction( 8, 0.7101  );
  
  AddMaterial("G4_MUSCLE_WITHOUT_SUCROSE", 1.07, 0, 74.2, 4);
  AddElementByWeightFraction( 1, 0.101969);
  AddElementByWeightFraction( 6, 0.120058);
  AddElementByWeightFraction( 7, 0.035451);
  AddElementByWeightFraction( 8, 0.742522);

  AddMaterial("G4_NAPHTHALENE", 1.145, 0, 68.4, 2);
  AddElementByAtomCount("C" , 10);
  AddElementByAtomCount("H" ,  8);

  AddMaterial("G4_NITROBENZENE", 1.19867, 0, 75.8, 4);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" ,  5);
  AddElementByAtomCount("N" ,  1);
  AddElementByAtomCount("O" ,  2);

  AddMaterial("G4_NITROUS_OXIDE", 0.00183094, 0, 84.9, 2, kStateGas);
  AddElementByAtomCount("N" ,  2);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_NYLON-8062", 1.08, 0, 64.3, 4);
  AddElementByWeightFraction( 1, 0.103509);
  AddElementByWeightFraction( 6, 0.648415);
  AddElementByWeightFraction( 7, 0.099536);
  AddElementByWeightFraction( 8, 0.148539);

  AddMaterial("G4_NYLON-6-6", 1.14, 0, 63.9, 4);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" , 11);
  AddElementByAtomCount("N" ,  1);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_NYLON-6-10", 1.14, 0, 63.2, 4);
  AddElementByWeightFraction( 1, 0.107062);
  AddElementByWeightFraction( 6, 0.680449);
  AddElementByWeightFraction( 7, 0.099189);
  AddElementByWeightFraction( 8, 0.1133  );

  AddMaterial("G4_NYLON-11_RILSAN", 1.425, 0, 61.6, 4);
  AddElementByWeightFraction( 1, 0.115476);
  AddElementByWeightFraction( 6, 0.720819);
  AddElementByWeightFraction( 7, 0.076417);
  AddElementByWeightFraction( 8, 0.087289);

  AddMaterial("G4_OCTANE", 0.7026, 0, 54.7, 2);
  AddElementByAtomCount("C" ,  8);
  AddElementByAtomCount("H" , 18);

  AddMaterial("G4_PARAFFIN", 0.93, 0, 55.9, 2);
  AddElementByAtomCount("C" , 25);
  AddElementByAtomCount("H" , 52);

  AddMaterial("G4_N-PENTANE", 0.6262, 0, 53.6, 2);
  AddElementByAtomCount("C" ,  5);
  AddElementByAtomCount("H" , 12);

  AddMaterial("G4_PHOTO_EMULSION", 3.815, 0, 331., 8);
  AddElementByWeightFraction( 1, 0.0141  );
  AddElementByWeightFraction( 6, 0.072261);
  AddElementByWeightFraction( 7, 0.01932 );
  AddElementByWeightFraction( 8, 0.066101);
  AddElementByWeightFraction(16, 0.00189 );
  AddElementByWeightFraction(35, 0.349103);
  AddElementByWeightFraction(47, 0.474105);
  AddElementByWeightFraction(53, 0.00312 );

  AddMaterial("G4_PLASTIC_SC_VINYLTOLUENE", 1.032, 0, 64.7, 2);
  // AddElementByWeightFraction( 1, 0.085);
  // AddElementByWeightFraction( 6, 0.915);
  // Watch out! These weight fractions do not correspond to pure PVT 
  // (PolyVinylToluene, C_9H_10) but to an unknown mixture...
  // M.Trocme & S.Seltzer
  AddElementByAtomCount("C" ,  9);
  AddElementByAtomCount("H" , 10);

  AddMaterial("G4_PLUTONIUM_DIOXIDE", 11.46, 0, 746.5, 2);
  AddElementByAtomCount("Pu",  1);
  AddElementByAtomCount("O" ,  2);

  AddMaterial("G4_POLYACRYLONITRILE", 1.17, 0, 69.6, 3);
  AddElementByAtomCount("C" ,  3);
  AddElementByAtomCount("H" ,  3);
  AddElementByAtomCount("N" ,  1);

  AddMaterial("G4_POLYCARBONATE", 1.2, 0, 73.1, 3);
  AddElementByAtomCount("C" , 16);
  AddElementByAtomCount("H" , 14);
  AddElementByAtomCount("O" ,  3);

  AddMaterial("G4_POLYCHLOROSTYRENE", 1.3, 0, 81.7, 3);
  //  AddElementByWeightFraction( 1, 0.061869);
  //  AddElementByWeightFraction( 6, 0.696325);
  //  AddElementByWeightFraction(17, 0.241806);
  //  These weight fractions correspond to C_17H_18Cl_2 which is not 
  //  POLYCHLOROSTYRENE. POLYCHLOROSTYRENE is C_8H_7Cl!
  //  M.Trocme & S.Seltzer
  AddElementByAtomCount("C" ,  8);
  AddElementByAtomCount("H" ,  7);
  AddElementByAtomCount("Cl",  1);

  AddMaterial("G4_POLYETHYLENE", 0.94, 0, 57.4, 2);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("H" ,  2);
  chFormulas[nMaterials-1] = "(C_2H_4)_N-Polyethylene";

  AddMaterial("G4_MYLAR", 1.4, 0, 78.7, 3);
  AddElementByAtomCount("C" , 10);
  AddElementByAtomCount("H" ,  8);
  AddElementByAtomCount("O" ,  4);

  AddMaterial("G4_PLEXIGLASS", 1.19, 0, 74., 3);
  AddElementByAtomCount("C" ,  5);
  AddElementByAtomCount("H" ,  8);
  AddElementByAtomCount("O" ,  2);

  AddMaterial("G4_POLYOXYMETHYLENE", 1.425 ,0, 77.4, 3);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("H" ,  2);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_POLYPROPYLENE", 0.9, 0, 56.5, 2);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  4);
  chFormulas[nMaterials-1] = "(C_2H_4)_N-Polypropylene";

  AddMaterial("G4_POLYSTYRENE", 1.06, 0, 68.7, 2);
  AddElementByAtomCount("C" ,  8);
  AddElementByAtomCount("H" ,  8);

  AddMaterial("G4_TEFLON", 2.2, 0, 99.1, 2);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("F" ,  4);

  AddMaterial("G4_POLYTRIFLUOROCHLOROETHYLENE", 2.1, 0, 120.7, 3);
  // correct chemical name Polychlorotrifluoroethylene [CF2CClF]n, IvantchenkoA.
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("F" ,  3);
  AddElementByAtomCount("Cl",  1);

  AddMaterial("G4_POLYVINYL_ACETATE", 1.19, 0, 73.7, 3);
  AddElementByAtomCount("C" ,  4);
  AddElementByAtomCount("H" ,  6);
  AddElementByAtomCount("O" ,  2);

  AddMaterial("G4_POLYVINYL_ALCOHOL", 1.3, 0, 69.7, 3);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  4);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_POLYVINYL_BUTYRAL", 1.12, 0, 67.2, 3);
  //  AddElementByWeightFraction( 1, 0.092802);
  //  AddElementByWeightFraction( 6, 0.680561);
  //  AddElementByWeightFraction( 8, 0.226637);
  //  These weight fractions correspond to C_8H_13O_2 which is not 
  //  POLYVINYL_BUTYRAL. POLYVINYL_BUTYRAL is C_8H_14O_2!
  //  M.Trocme & S.Seltzer
  AddElementByAtomCount("C" ,  8);
  AddElementByAtomCount("H" , 14);
  AddElementByAtomCount("O" ,  2);

  AddMaterial("G4_POLYVINYL_CHLORIDE", 1.3, 0, 108.2, 3);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  3);
  AddElementByAtomCount("Cl",  1);

  AddMaterial("G4_POLYVINYLIDENE_CHLORIDE", 1.7, 0, 134.3, 3);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  2);
  AddElementByAtomCount("Cl",  2);

  AddMaterial("G4_POLYVINYLIDENE_FLUORIDE", 1.76, 0, 88.8, 3);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  2);
  AddElementByAtomCount("F" ,  2);

  AddMaterial("G4_POLYVINYL_PYRROLIDONE", 1.25, 0, 67.7, 4);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" ,  9);
  AddElementByAtomCount("N" ,  1);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_POTASSIUM_IODIDE", 3.13, 0, 431.9, 2);
  AddElementByAtomCount("K" ,  1);
  AddElementByAtomCount("I" ,  1);

  AddMaterial("G4_POTASSIUM_OXIDE", 2.32, 0, 189.9, 2);
  AddElementByAtomCount("K" ,  2);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_PROPANE", 0.00187939, 0, 47.1, 2, kStateGas);
  AddElementByAtomCount("C" ,  3);
  AddElementByAtomCount("H" ,  8);

  AddMaterial("G4_lPROPANE", 0.43, 0, 52., 2);
  AddElementByAtomCount("C" ,  3);
  AddElementByAtomCount("H" ,  8);

  AddMaterial("G4_N-PROPYL_ALCOHOL", 0.8035, 0, 61.1, 3);
  AddElementByAtomCount("C" ,  3);
  AddElementByAtomCount("H" ,  8);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_PYRIDINE", 0.9819, 0, 66.2, 3);
  AddElementByAtomCount("C" ,  5);
  AddElementByAtomCount("H" ,  5);
  AddElementByAtomCount("N" ,  1);

  AddMaterial("G4_RUBBER_BUTYL", 0.92, 0, 56.5, 2);
  AddElementByWeightFraction( 1, 0.143711);
  AddElementByWeightFraction( 6, 0.856289);

  AddMaterial("G4_RUBBER_NATURAL", 0.92, 0, 59.8, 2);
  AddElementByWeightFraction( 1, 0.118371);
  AddElementByWeightFraction( 6, 0.881629);

  AddMaterial("G4_RUBBER_NEOPRENE", 1.23, 0, 93., 3);
  AddElementByWeightFraction( 1, 0.05692 );
  AddElementByWeightFraction( 6, 0.542646);
  AddElementByWeightFraction(17, 0.400434);

  AddMaterial("G4_SILICON_DIOXIDE", 2.32, 0, 139.2, 2);
  AddElementByAtomCount("Si",  1);
  AddElementByAtomCount("O" ,  2);
  chFormulas[nMaterials-1] = "SiO_2";

  AddMaterial("G4_SILVER_BROMIDE", 6.473, 0, 486.6, 2);
  AddElementByAtomCount("Ag",  1);
  AddElementByAtomCount("Br",  1);

  AddMaterial("G4_SILVER_CHLORIDE", 5.56, 0, 398.4, 2);
  AddElementByAtomCount("Ag",  1);
  AddElementByAtomCount("Cl",  1);

  AddMaterial("G4_SILVER_HALIDES", 6.47, 0, 487.1, 3);
  AddElementByWeightFraction(35, 0.422895);
  AddElementByWeightFraction(47, 0.573748);
  AddElementByWeightFraction(53, 0.003357);

  AddMaterial("G4_SILVER_IODIDE", 6.01, 0, 543.5, 2);
  AddElementByAtomCount("Ag",  1);
  AddElementByAtomCount("I" ,  1);

  AddMaterial("G4_SKIN_ICRP", 1.09, 0, 72.7, 9);
  AddElementByWeightFraction( 1, 0.100);
  AddElementByWeightFraction( 6, 0.204);
  AddElementByWeightFraction( 7, 0.042);
  AddElementByWeightFraction( 8, 0.645);
  AddElementByWeightFraction(11, 0.002);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.003);
  AddElementByWeightFraction(19, 0.001);

  AddMaterial("G4_SODIUM_CARBONATE", 2.532, 0, 125., 3);
  AddElementByAtomCount("Na",  2);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("O" ,  3);

  AddMaterial("G4_SODIUM_IODIDE", 3.667, 0, 452., 2);
  AddElementByAtomCount("Na",  1);
  AddElementByAtomCount("I" ,  1);

  AddMaterial("G4_SODIUM_MONOXIDE", 2.27, 0, 148.8, 2);
  AddElementByAtomCount("Na",  2);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_SODIUM_NITRATE", 2.261, 0, 114.6, 3);
  AddElementByAtomCount("Na",  1);
  AddElementByAtomCount("N" ,  1);
  AddElementByAtomCount("O" ,  3);

  AddMaterial("G4_STILBENE", 0.9707, 0, 67.7, 2);
  AddElementByAtomCount("C" , 14);
  AddElementByAtomCount("H" , 12);

  AddMaterial("G4_SUCROSE", 1.5805, 0, 77.5, 3);
  AddElementByAtomCount("C" , 12);
  AddElementByAtomCount("H" , 22);
  AddElementByAtomCount("O" , 11);

  AddMaterial("G4_TERPHENYL", 1.24 /*1.234*/, 0, 71.7, 2);
  //  AddElementByWeightFraction( 1, 0.044543);
  //  AddElementByWeightFraction( 6, 0.955457);
  //  These weight fractions correspond to C_18H_10 which is not TERPHENYL. 
  //  TERPHENYL is C_18H_14! The current density is 1.24 g/cm3
  //  M.Trocme & S.Seltzer
  AddElementByAtomCount("C" , 18);
  AddElementByAtomCount("H" , 14);

  AddMaterial("G4_TESTIS_ICRP", 1.04, 0, 75., 9);
  AddElementByWeightFraction( 1, 0.106);
  AddElementByWeightFraction( 6, 0.099);
  AddElementByWeightFraction( 7, 0.020);
  AddElementByWeightFraction( 8, 0.766);
  AddElementByWeightFraction(11, 0.002);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.002);
  AddElementByWeightFraction(19, 0.002);

  AddMaterial("G4_TETRACHLOROETHYLENE", 1.625, 0, 159.2, 2);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("Cl",  4);

  AddMaterial("G4_THALLIUM_CHLORIDE", 7.004, 0, 690.3, 2);
  AddElementByAtomCount("Tl",  1);
  AddElementByAtomCount("Cl",  1);

  // TISSUE_SOFT_MALE ICRU-44/46 (1989)
  AddMaterial("G4_TISSUE_SOFT_ICRP", 1.03, 0, 72.3, 9);
  AddElementByWeightFraction( 1, 0.105);
  AddElementByWeightFraction( 6, 0.256);
  AddElementByWeightFraction( 7, 0.027);
  AddElementByWeightFraction( 8, 0.602);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.002);
  AddElementByWeightFraction(16, 0.003);
  AddElementByWeightFraction(17, 0.002);
  AddElementByWeightFraction(19, 0.002);

  // Tissue soft adult ICRU-33 (1980)
  AddMaterial("G4_TISSUE_SOFT_ICRU-4", 1.0, 0, 74.9, 4);
  AddElementByWeightFraction( 1, 0.101);
  AddElementByWeightFraction( 6, 0.111);
  AddElementByWeightFraction( 7, 0.026);
  AddElementByWeightFraction( 8, 0.762);

  AddMaterial("G4_TISSUE-METHANE", 0.00106409, 0, 61.2, 4, kStateGas);
  AddElementByWeightFraction( 1, 0.101869);
  AddElementByWeightFraction( 6, 0.456179);
  AddElementByWeightFraction( 7, 0.035172);
  AddElementByWeightFraction( 8, 0.40678 );

  AddMaterial("G4_TISSUE-PROPANE", 0.00182628, 0, 59.5, 4, kStateGas);
  AddElementByWeightFraction( 1, 0.102672);
  AddElementByWeightFraction( 6, 0.56894 );
  AddElementByWeightFraction( 7, 0.035022);
  AddElementByWeightFraction( 8, 0.293366);

  AddMaterial("G4_TITANIUM_DIOXIDE", 4.26, 0, 179.5, 2);
  AddElementByAtomCount("Ti",  1);
  AddElementByAtomCount("O" ,  2);

  AddMaterial("G4_TOLUENE", 0.8669, 0, 62.5, 2);
  AddElementByAtomCount("C" ,  7);
  AddElementByAtomCount("H" ,  8);

  AddMaterial("G4_TRICHLOROETHYLENE", 1.46, 0, 148.1, 3);
  AddElementByAtomCount("C" ,  2);
  AddElementByAtomCount("H" ,  1);
  AddElementByAtomCount("Cl",  3);

  AddMaterial("G4_TRIETHYL_PHOSPHATE", 1.07, 0, 81.2, 4);
  AddElementByAtomCount("C" ,  6);
  AddElementByAtomCount("H" , 15);
  AddElementByAtomCount("O" ,  4);
  AddElementByAtomCount("P" ,  1);

  AddMaterial("G4_TUNGSTEN_HEXAFLUORIDE", 2.4, 0, 354.4, 2);
  AddElementByAtomCount("W" ,  1);
  AddElementByAtomCount("F" ,  6);

  AddMaterial("G4_URANIUM_DICARBIDE", 11.28, 0, 752., 2);
  AddElementByAtomCount("U" ,  1);
  AddElementByAtomCount("C" ,  2);

  AddMaterial("G4_URANIUM_MONOCARBIDE", 13.63, 0, 862., 2);
  AddElementByAtomCount("U" ,  1);
  AddElementByAtomCount("C" ,  1);

  AddMaterial("G4_URANIUM_OXIDE", 10.96, 0, 720.6, 2);
  AddElementByAtomCount("U" ,  1);
  AddElementByAtomCount("O" ,  2);

  AddMaterial("G4_UREA", 1.323, 0, 72.8, 4);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("H" ,  4);
  AddElementByAtomCount("N" ,  2);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_VALINE", 1.23, 0, 67.7, 4);
  AddElementByAtomCount("C" ,  5);
  AddElementByAtomCount("H" , 11);
  AddElementByAtomCount("N" ,  1);
  AddElementByAtomCount("O" ,  2);

  AddMaterial("G4_VITON", 1.8, 0, 98.6, 3);
  AddElementByWeightFraction( 1, 0.009417);
  AddElementByWeightFraction( 6, 0.280555);
  AddElementByWeightFraction( 9, 0.710028);

  AddMaterial("G4_WATER_VAPOR", 0.000756182, 0, 71.6, 2, kStateGas);
  AddElementByAtomCount("H" ,  2);
  AddElementByAtomCount("O" ,  1);
  chFormulas[nMaterials-1] = "H_2O-Gas";

  AddMaterial("G4_XYLENE", 0.87, 0, 61.8, 2);
  AddElementByAtomCount("C" ,  8);
  AddElementByAtomCount("H" , 10);

  AddMaterial("G4_GRAPHITE", 2.21, 6, 78.);
  chFormulas[nMaterials-1] = "Graphite";

  nNIST = nMaterials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::HepAndNuclearMaterials()
{
  AddMaterial("G4_lH2", 0.0708,  1,  21.8, 1, kStateLiquid, false);
  AddMaterial("G4_lN2", 0.807,   7,  82.,  1, kStateLiquid, false);
  AddMaterial("G4_lO2", 1.141,   8,  95.,  1, kStateLiquid, false);
  AddMaterial("G4_lAr", 1.396 , 18, 188. , 1, kStateLiquid, false);
  AddMaterial("G4_lBr", 3.1028, 35, 343. , 1, kStateLiquid, false);
  AddMaterial("G4_lKr", 2.418 , 36, 352. , 1, kStateLiquid, false);
  AddMaterial("G4_lXe", 2.953 , 54, 482. , 1, kStateLiquid, false);

  AddMaterial("G4_PbWO4", 8.28, 0, 0.0, 3);
  AddElementByAtomCount("O" , 4);
  AddElementByAtomCount("Pb", 1);
  AddElementByAtomCount("W" , 1);

  G4double density = universe_mean_density*cm3/g;
  AddMaterial("G4_Galactic", density, 1, 21.8, 1, kStateGas);
  AddGas("G4_Galactic",2.73*kelvin, 3.e-18*hep_pascal);

  AddMaterial("G4_GRAPHITE_POROUS", 1.7, 6, 78.);
  chFormulas[nMaterials-1] = "Graphite";

  // LUCITE is equal to plustiglass
  AddMaterial("G4_LUCITE", 1.19, 0, 74., 3);
  AddElementByWeightFraction( 1, 0.080538);
  AddElementByWeightFraction( 6, 0.599848);
  AddElementByWeightFraction( 8, 0.319614);

  // SRIM-2008 materials
  AddMaterial("G4_BRASS", 8.52, 0, 0.0, 3);
  AddElementByAtomCount("Cu", 62);
  AddElementByAtomCount("Zn", 35);
  AddElementByAtomCount("Pb" , 3);

  AddMaterial("G4_BRONZE", 8.82, 0, 0.0, 3);
  AddElementByAtomCount("Cu", 89);
  AddElementByAtomCount("Zn",  9);
  AddElementByAtomCount("Pb" , 2);

  // parameters are taken from
  //  http://www.azom.com/article.aspx?ArticleID=965
  AddMaterial("G4_STAINLESS-STEEL", 8.00, 0, 0.0, 3);
  AddElementByAtomCount("Fe", 74);
  AddElementByAtomCount("Cr", 18);
  AddElementByAtomCount("Ni" , 8);

  AddMaterial("G4_CR39", 1.32, 0, 0.0, 3);
  AddElementByAtomCount("H", 18);
  AddElementByAtomCount("C", 12);
  AddElementByAtomCount("O", 7);

  AddMaterial("G4_OCTADECANOL", 0.812, 0, 0.0, 3);
  AddElementByAtomCount("H", 38);
  AddElementByAtomCount("C", 18);
  AddElementByAtomCount("O", 1);

  nHEP = nMaterials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::SpaceMaterials()
{
  // density in g/cm3
  AddMaterial("G4_KEVLAR" , 1.44, 0, 0.0, 4);
  AddElementByAtomCount("C", 14);
  AddElementByAtomCount("H", 10);
  AddElementByAtomCount("O", 2);
  AddElementByAtomCount("N", 2);

  AddMaterial("G4_DACRON" , 1.40, 0, 0.0, 3);   // G4_POLYETHYLENE_TEREPHTALATE
  AddElementByAtomCount("C", 10);
  AddElementByAtomCount("H", 8);
  AddElementByAtomCount("O", 4);

  AddMaterial("G4_NEOPRENE" , 1.23, 0, 0.0, 3);   // POLYCLOROPRENE
  AddElementByAtomCount("C", 4);
  AddElementByAtomCount("H", 5);
  AddElementByAtomCount("Cl",1);

  nSpace = nMaterials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
  Build biochemical materials used in G4DNA Applications.
  Materials are defined in bonded and unbonded forms according to the
  following schema:
  G4_MATERIAL: Molecule in its free state
  G4_DNA_MATERIAL: Molecule, considering atoms lost in bonding
*/

void G4NistMaterialBuilder::BioChemicalMaterials()
{
  // BEGIN MATERIALS IN THEIR UNBONDED FORM

  // G4_ADENINE, G4_GUANINE are defined in
  // G4NistMaterialBuilder::NistCompoundMaterials()

  //Tan, Z., et al. NIMB,2006(248)
  AddMaterial("G4_CYTOSINE", 1.3, 0, 72., 4);
  AddElementByAtomCount("H", 5);
  AddElementByAtomCount("C", 4);
  AddElementByAtomCount("N", 3);
  AddElementByAtomCount("O", 1);

  AddMaterial("G4_THYMINE", 1.48, 0, 72., 4);
  AddElementByAtomCount("H", 6);
  AddElementByAtomCount("C", 5);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 2);

  AddMaterial("G4_URACIL", 1.32, 0, 72., 4);
  AddElementByAtomCount("H", 4);
  AddElementByAtomCount("C", 4);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 2);

  //ACD Labs Percepta Plateform - PhysChem Module
  //(https://www.acdlabs.com/products/percepta/index.php)
  AddMaterial("G4_DEOXYRIBOSE", 1.5, 0, 72, 3);
  AddElementByAtomCount("H", 10);
  AddElementByAtomCount("C", 5);
  AddElementByAtomCount("O", 4);

  //Egan, E.P. and B.B. Luff,
  //Industrial & Engineering Chemistry, 1955. 47(6): p. 1280-1281
  AddMaterial("G4_PHOSPHORIC_ACID", 1.87, 0, 72, 3);
  AddElementByAtomCount("H", 3);
  AddElementByAtomCount("P", 1);
  AddElementByAtomCount("O", 4);
    
  // END UNBONDED MATERIALS / BEGIN BONDED MATERIALS

  // Deoxyribose loses 3 OH groups in bonding to bond with PO4 and a base pair
  AddMaterial("G4_DNA_DEOXYRIBOSE", 1, 0, 72., 3);
  AddElementByAtomCount("H", 7);
  AddElementByAtomCount("C", 5);
  AddElementByAtomCount("O", 1);

  // Typically there are no H atoms considered in the Phosphate group
  AddMaterial("G4_DNA_PHOSPHATE", 1, 0, 72., 2);
  AddElementByAtomCount("P", 1);
  AddElementByAtomCount("O", 4);

  // GATCU bases bonded to a deoxyribose (they drop one H)
  AddMaterial("G4_DNA_ADENINE", 1, 0, 72., 3);
  AddElementByAtomCount("H",4 );
  AddElementByAtomCount("C",5 );
  AddElementByAtomCount("N",5 );

  AddMaterial("G4_DNA_GUANINE", 1, 0, 72., 4);
  AddElementByAtomCount("H",4 );
  AddElementByAtomCount("C",5 );
  AddElementByAtomCount("N",5 );
  AddElementByAtomCount("O",1 );

  AddMaterial("G4_DNA_CYTOSINE", 1, 0, 72., 4);
  AddElementByAtomCount("H", 4);
  AddElementByAtomCount("C", 4);
  AddElementByAtomCount("N", 3);
  AddElementByAtomCount("O", 1);

  AddMaterial("G4_DNA_THYMINE", 1, 0, 72., 4);
  AddElementByAtomCount("H", 5);
  AddElementByAtomCount("C", 5);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 2);

  AddMaterial("G4_DNA_URACIL", 1, 0, 72., 4);
  AddElementByAtomCount("H", 3);
  AddElementByAtomCount("C", 4);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 2);

  // END BONDED MATERIALS

  /*
  // Complete 70 kg body of adult men from en.wikipedia.org/ see References there
  AddMaterial("G4_BODY", 1.8, 0, 78, 12);
  AddElementByWeightFraction( 8, 0.650);
  AddElementByWeightFraction( 6, 0.180);
  AddElementByWeightFraction( 1, 0.100);
  AddElementByWeightFraction( 7, 0.030);
  AddElementByWeightFraction(20, 0.015);
  AddElementByWeightFraction(15, 0.010);
  AddElementByWeightFraction(19, 0.0025);
  AddElementByWeightFraction(16, 0.0025);
  AddElementByWeightFraction(11, 0.0015);
  AddElementByWeightFraction(17, 0.0015);
  AddElementByWeightFraction(12, 0.0005);
  AddElementByWeightFraction(26, 0.00006);
  */
}




