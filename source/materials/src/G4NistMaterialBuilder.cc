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
// $Id: G4NistMaterialBuilder.cc 67044 2013-01-30 08:50:06Z gcosmo $
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
// 27.02.06 V.Ivanchneko add ConstructNewGasMaterial
// 11.05.06 V.Ivanchneko add warning flag to FindMaterial method
// 27.06.06 V.Ivanchneko fix graphite description
// 27.07.07 V.Ivanchneko remove dependence on NistManager
// 30.10.09 V.Ivanchneko update density of G4_GRAFITE from PDG'2008
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
//
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NistMaterialBuilder::G4NistMaterialBuilder(G4NistElementBuilder* eb, G4int vb)
: elmBuilder(eb),
  verbose(vb),
  nMaterials(0),
  nComponents(0),
  nCurrent(0),
  first(true)
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NistMaterialBuilder::~G4NistMaterialBuilder()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::FindOrBuildMaterial(const G4String& matname,
                                                       G4bool isotopes,
						       G4bool warning)
{
  if(first) {
    if(verbose > 0) {
      G4cout << "### NIST DataBase for Materials is used" << G4endl;
    }
    first = false;
  }

  G4String name = matname;
  if("G4_NYLON-6/6" == matname)  { name = "G4_NYLON-6-6"; }
  if("G4_NYLON-6/10" == matname) { name = "G4_NYLON-6-10";}

  if (verbose > 1) {
    G4cout << "G4NistMaterialBuilder::FindOrBuildMaterial " << name << G4endl;
  }
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int nmat = theMaterialTable->size();

  // Check if name inside DB
  G4Material* mat = 0;

  for (G4int i=0; i<nMaterials; ++i) {

    if (name == names[i]) {
      // Build new Nist material 
      if(matIndex[i] == -1) { mat = BuildMaterial(i, isotopes); }
      // Nist material was already built
      else                  { mat = (*theMaterialTable)[matIndex[i]]; }
      return mat;
    }
  }

  // Check the list of all materials
  if (nmat > 0) {
    for (G4int i=0; i<nmat; ++i) {
      if(name == ((*theMaterialTable)[i])->GetName()) {
        mat = (*theMaterialTable)[i];
	return mat;
      }
    }
  }

  if( (verbose == 1 && warning) || verbose > 1) {
    G4cout << "G4NistMaterialBuilder::FindOrBuildMaterial WARNING:"
	   << " material <" << name
	   << "> is not found out" << G4endl;
  }
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::BuildMaterial(G4int i, G4bool isotopes)
{
  if (verbose > 1) {
    G4cout << "G4NistMaterialBuilder: BuildMaterial #" << i
	   << G4endl;
  }
  G4Material* mat = 0;
  if (nMaterials == 0) { return mat; }

  G4int nc = components[i];

  // Check gas parameters
  G4double t = STP_Temperature;
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
    // liquids
  } else if( !STP[i] ) { t = 0.0; }

  mat = new G4Material(names[i],densities[i],nc,states[i],t,p);

  if (verbose>1) { G4cout << "New material nComponents= " << nc << G4endl; }
  if (nc > 0) {
    G4int idx = indexes[i];
    for (G4int j=0; j<nc; ++j) {
      G4int Z = elements[idx+j];
      G4Element* el = elmBuilder->FindOrBuildElement(Z, isotopes);
      if(!el) {
	G4cout << "G4NistMaterialBuilder::BuildMaterial:"
	       << "  ERROR: elements Z= " << Z << " is not found "
	       << " for material " << names[i]
	       << G4endl;
	G4Exception("G4NistMaterialBuilder::BuildMaterial()", "mat103",
	             FatalException, "Fail to construct material");
	return 0;
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
  if(chFormulas[i] != "") {
    mat->SetChemicalFormula(chFormulas[i]);
    exc1 = ion->FindMeanExcitationEnergy(chFormulas[i]);
  }
  if(ionPotentials[i] > 0.0) { exc1 = ionPotentials[i]; }
  if(exc0 != exc1) { ion->SetMeanExcitationEnergy(exc1); }

  // Index in Material Table
  matIndex[i] = mat->GetIndex();
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4int>& nbAtoms,
				      G4double dens, 
				      G4bool isotopes,
				      G4State state,     
				      G4double temp,  
				      G4double pres)
{
  // Material is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat) { 
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: the material <" << name
	   << "> is already exist" << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return mat; 
  }

  // Material not in DB
  G4int els = elm.size();
  if(els == 0) { 
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: empty list of elements for " << name
	   << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return 0;
  } 

  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential is not defined
  G4bool stp = true;
  if(state == kStateGas && temp != STP_Temperature && pres != STP_Pressure)
    { stp = false; }

  AddMaterial(name,dens*cm3/g,0,0.,els,state,stp);
  if(!stp) { AddGas(name,temp,pres); }

  for (G4int i=0; i<els; ++i) {
    AddElementByAtomCount(elmBuilder->GetZ(elm[i]), nbAtoms[i]);
  }

  return BuildMaterial(nMaterials-1, isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4double>& w,
				      G4double dens, 
				      G4bool isotopes,
				      G4State state,     
				      G4double temp,  
				      G4double pres)
{
  // Material is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat) { 
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: the material <" << name
	   << "> is already exist" << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return mat; 
  }

  // Material not in DB
  G4int els = elm.size();
  if(els == 0) { 
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: empty list of elements for " << name
	   << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return 0;
  } 

  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential is not defined
  G4bool stp = true;
  if(state == kStateGas && temp != STP_Temperature && pres != STP_Pressure)
    { stp = false; }
  AddMaterial(name,dens*cm3/g,0,0.,els,state,stp);
  if(!stp) { AddGas(name,temp,pres); }

  for (G4int i=0; i<els; ++i) {
    AddElementByWeightFraction(elmBuilder->GetZ(elm[i]), w[i]);
  }
  
  return BuildMaterial(nMaterials-1, isotopes);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialBuilder::ConstructNewGasMaterial(
				      const G4String& name,
				      const G4String& nameDB,
				      G4double temp, 
				      G4double pres, 
				      G4bool)
{
  // Material name is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat) { 
    G4cout << "G4NistMaterialBuilder::ConstructNewGasMaterial:"
           << "  WARNING: the material <" << name
	   << "> is already exist" << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return mat; 
  }

  G4Material* bmat = FindOrBuildMaterial(nameDB);
  if(!bmat) {
    G4cout << "G4NistMaterialBuilder::ConstructNewGasMaterial:"
	   << "  WARNING: the Name <" << nameDB 
	   << "> is NOT in the DB: no new gas will be constructed"
	   << G4endl;
    return 0;
  }
  if(bmat->GetState() != kStateGas) {
    G4cout << "G4NistMaterialBuilder::ConstructNewGasMaterial:"
	   << "  WARNING:  <" << nameDB 
	   << "> is NOT a gas -  no new gas will be constructed"
	   << G4endl;
    return 0;
  }

  G4double dens = bmat->GetDensity()*pres*STP_Temperature/(temp*STP_Pressure);
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
                                      G4bool isotopes,
                                      G4double temp,
                                      G4double pres)
{
  G4State state = kStateGas;

  // Material is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat) {
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: the material <" << name
           << "> is already exist" << G4endl;
    G4cout << "      New material will NOT be built!"
           << G4endl;
    return mat;
  }

  // Material not in DB
  G4int els = elm.size();
  if(els == 0) {
    G4cout << "G4NistMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: empty list of elements for " << name
           << G4endl;
    G4cout << "      New material will NOT be built!"
           << G4endl;
    return 0;
  }

  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential is not defined
  G4bool stp = true;
  if(temp != STP_Temperature && pres != STP_Pressure)
    { stp = false; }

  G4double massPerMole = 0;

  G4int Z = 0;
  for (G4int i=0; i<els; ++i) {
    Z = elmBuilder->GetZ(elm[i]);
    massPerMole += nbAtoms[i] * elmBuilder->GetAtomicMassAmu(Z) * amu_c2;
  }

  G4double dens = massPerMole / (Avogadro*k_Boltzmann*temp/pres);

  if (els == 1) { AddMaterial(name,dens,Z,0.,els,state,stp); }
  else {
    AddMaterial(name,dens,0,0.,els,state,stp);
    for (G4int i=0; i<els; ++i) {
      AddElementByAtomCount(elmBuilder->GetZ(elm[i]), nbAtoms[i]);
    }
  }

  if(!stp) { AddGas(name,temp,pres); }

  return BuildMaterial(nMaterials-1, isotopes);
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
    G4cout << "         New material " << nameMat << " will not be added" 
	   << G4endl;
    return;
  }

  // density in g/cm3, mean ionisation potential in eV

  names.push_back(nameMat);
  chFormulas.push_back("");
  densities.push_back(dens*g/cm3);
  ionPotentials.push_back(pot*eV);
  states.push_back(state);
  components.push_back(ncomp);
  indexes.push_back(nComponents);
  STP.push_back(stp);
  matIndex.push_back(-1);
  atomCount.push_back(false);

  if (ncomp == 1 && Z > 0) {
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
    G4cout << "New material " << nameMat << " is prepeared; "
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

void G4NistMaterialBuilder::ListMaterials(const G4String& mnam)
{
  if (mnam == "simple")           { ListNistSimpleMaterials(); }
  else if (mnam == "compound")    { ListNistCompoundMaterials(); }
  else if (mnam == "hep")         { ListHepMaterials(); }
  else if (mnam == "space")       { ListSpaceMaterials(); }
  else if (mnam == "biochemical") { ListBioChemicalMaterials(); }

  else if (mnam == "all") {
    ListNistSimpleMaterials();
    ListNistCompoundMaterials();
    ListHepMaterials();
    ListSpaceMaterials();
    ListBioChemicalMaterials();

  } else {
    G4cout << "### G4NistMaterialBuilder::ListMaterials: Warning " 
	   << mnam << " list is not known" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListNistSimpleMaterials()
{
  G4cout << "=======================================================" << G4endl;
  G4cout << "###   Simple Materials from the NIST Data Base   ###" << G4endl;
  G4cout << "=======================================================" << G4endl;
  G4cout << " Z Name  ChFormula        density(g/cm^3)  I(eV)       " << G4endl;
  G4cout << "=======================================================" << G4endl;
  for (G4int i=0; i<nElementary; ++i) {DumpElm(i);}
  G4cout << "=======================================================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListNistCompoundMaterials()
{
  G4cout << "###    Compound Materials from the NIST Data Base    ##" << G4endl;
  G4cout << "=======================================================" << G4endl;
  G4cout << " Ncomp Name  ChFormula        density(g/cm^3)  I(eV)   " << G4endl;
  G4cout << "=======================================================" << G4endl;
  for (G4int i=nElementary; i<nNIST; ++i) {DumpMix(i);}
  G4cout << "=======================================================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListHepMaterials()
{
  G4cout << "=======================================================" << G4endl;
  G4cout << "###           HEP & Nuclear Materials                ##" << G4endl;
  G4cout << "=======================================================" << G4endl;
  G4cout << " Ncomp Name  ChFormula        density(g/cm^3)  I(eV)   " << G4endl;
  G4cout << "=======================================================" << G4endl;
  for (G4int i=nNIST; i<nHEP; ++i) {DumpMix(i);}
  G4cout << "=======================================================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListSpaceMaterials()
{
  G4cout << "=======================================================" << G4endl;
  G4cout << "###           Space ISS Materials                    ##" << G4endl;
  G4cout << "=======================================================" << G4endl;
  G4cout << " Ncomp Name  ChFormula        density(g/cm^3)  I(eV)   " << G4endl;
  G4cout << "=======================================================" << G4endl;
  for (G4int i=nHEP; i<nSpace; ++i) {DumpMix(i);}
  G4cout << "=======================================================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::ListBioChemicalMaterials()
{
  G4cout << "=======================================================" << G4endl;
  G4cout << "###          Bio-Chemical Materials                  ##" << G4endl;
  G4cout << "=======================================================" << G4endl;
  G4cout << " Ncomp Name  ChFormula        density(g/cm^3)  I(eV)   " << G4endl;
  G4cout << "=======================================================" << G4endl;
  for (G4int i=nSpace; i<nMaterials; ++i) {DumpMix(i);}
  G4cout << "=======================================================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::DumpElm(G4int i)
{
  G4cout << i+1 << "  " << names[i] << "  " << chFormulas[i]
         << densities[i]*cm3/g << "     " << ionPotentials[i]/eV
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::DumpMix(G4int i)
{
  G4int nc = components[i];
  G4cout << nc << "  " << names[i] << "  " << chFormulas[i]
         << densities[i]*cm3/g << "     " << ionPotentials[i]/eV
	 << G4endl;
  if (nc > 1) {
    G4int imin = indexes[i];
    G4int imax = imin + nc;
    for (G4int j=imin; j<imax; ++j) {
      G4cout << "              " << elements[j] << "     " << fractions[j] 
             << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4NistMaterialBuilder::AddGas(const G4String& nameMat, G4double t, G4double p)
{
  G4int idx = nMaterials-1;
  if(nameMat != names[idx]) {
    idx = -1;
    for(G4int i=0; i<nMaterials; ++i) {
      if(nameMat == names[i]) {
        idx = i; break;
      }
    }
  }
  if(idx >= 0) {
    idxGas.push_back(idx);
    gasTemperature.push_back(t);
    gasPressure.push_back(p);
  } else {
    G4cout << "WARNING: G4NistMaterialBuilder::AddGas problem: there is no "
	   << nameMat << " in the list of materials;"
	   << G4endl;
  }
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
      if (sum > 0.0) for (G4int i=imin; i<imax; ++i) {fractions[i] /= sum;}
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
  HepAndNuclearMaterials();
  SpaceMaterials();
  BioChemicalMaterials();

  if (verbose > 1) { ListMaterials("all"); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::NistSimpleMaterials()
{
  // density in g/cm3, mean ionisation potential in eV

  AddMaterial("G4_H" ,  8.37480e-5,  1,  19.2, 1, kStateGas);
  AddMaterial("G4_He",  1.66322e-4,  2,  41.8, 1, kStateGas);
  AddMaterial("G4_Li",  0.534     ,  3,  40. );
  AddMaterial("G4_Be",  1.848     ,  4,  63.7);
  AddMaterial("G4_B" ,  2.37      ,  5,  76. );
  AddMaterial("G4_C" ,  2.        ,  6,  81. );
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
  AddElementByWeightFraction( 1, 0.104122);
  AddElementByWeightFraction( 6, 0.620405);
  AddElementByWeightFraction( 8, 0.275473);

  AddMaterial("G4_ACETYLENE", 0.0010967, 0, 58.2, 2, kStateGas);
  AddElementByWeightFraction( 1, 0.077418);
  AddElementByWeightFraction( 6, 0.922582);

  AddMaterial("G4_ADENINE", 1.6/*1.35*/, 0, 71.4, 3);
  AddElementByAtomCount("H",5 );
  AddElementByAtomCount("C",5 );
  AddElementByAtomCount("N",5 );

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
  AddElementByWeightFraction( 1, 0.07919 );
  AddElementByWeightFraction( 6, 0.404439);
  AddElementByWeightFraction( 7, 0.157213);
  AddElementByWeightFraction( 8, 0.359159);

  AddMaterial("G4_ALUMINUM_OXIDE", 3.97, 0, 145.2, 2);
  AddElementByWeightFraction( 8, 0.470749);
  AddElementByWeightFraction(13, 0.529251);
  chFormulas[nMaterials-1] = "Al_2O_3";

  AddMaterial("G4_AMBER", 1.1, 0, 63.2, 3);
  AddElementByWeightFraction( 1, 0.10593 );
  AddElementByWeightFraction( 6, 0.788973);
  AddElementByWeightFraction( 8, 0.105096);

  AddMaterial("G4_AMMONIA", 0.000826019, 0, 53.7, 2, kStateGas);
  AddElementByWeightFraction( 1, 0.177547);
  AddElementByWeightFraction( 7, 0.822453);

  AddMaterial("G4_ANILINE", 1.0235, 0, 66.2, 3);
  AddElementByWeightFraction( 1, 0.075759);
  AddElementByWeightFraction( 6, 0.773838);
  AddElementByWeightFraction( 7, 0.150403);

  AddMaterial("G4_ANTHRACENE", 1.283, 0, 69.5, 2);
  AddElementByWeightFraction( 1, 0.05655);
  AddElementByWeightFraction( 6, 0.94345);

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
  AddElementByWeightFraction( 9, 0.21672);
  AddElementByWeightFraction(56, 0.78328);

  AddMaterial("G4_BARIUM_SULFATE", 4.5, 0, 285.7, 3);
  AddElementByWeightFraction( 8,0.274212);
  AddElementByWeightFraction(16,0.137368);
  AddElementByWeightFraction(56,0.58842 );

  AddMaterial("G4_BENZENE", 0.87865, 0, 63.4, 2);
  AddElementByWeightFraction( 1, 0.077418);
  AddElementByWeightFraction( 6, 0.922582);

  AddMaterial("G4_BERYLLIUM_OXIDE", 3.01, 0, 93.2, 2);
  AddElementByWeightFraction( 4, 0.36032);
  AddElementByWeightFraction( 8, 0.63968);

  AddMaterial("G4_BGO", 7.13, 0, 534.1, 3);
  AddElementByWeightFraction( 8, 0.154126);
  AddElementByWeightFraction(32, 0.17482 );
  AddElementByWeightFraction(83, 0.671054);

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
  AddElementByWeightFraction( 5, 0.78261);
  AddElementByWeightFraction( 6, 0.21739);

  AddMaterial("G4_BORON_OXIDE", 1.812, 0, 99.6, 2);
  AddElementByWeightFraction( 5, 0.310551);
  AddElementByWeightFraction( 8, 0.689449);

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
  AddElementByWeightFraction( 1, 0.173408);
  AddElementByWeightFraction( 6, 0.826592);

  AddMaterial("G4_N-BUTYL_ALCOHOL", 0.8098, 0, 59.9, 3);
  AddElementByWeightFraction( 1, 0.135978);
  AddElementByWeightFraction( 6, 0.648171);
  AddElementByWeightFraction( 8, 0.215851);

  AddMaterial("G4_C-552", 1.76, 0, 86.8, 5);
  AddElementByWeightFraction( 1, 0.02468 );
  AddElementByWeightFraction( 6, 0.50161 );
  AddElementByWeightFraction( 8, 0.004527);
  AddElementByWeightFraction( 9, 0.465209);
  AddElementByWeightFraction(14, 0.003973);

  AddMaterial("G4_CADMIUM_TELLURIDE", 6.2, 0, 539.3, 2);
  AddElementByWeightFraction(48, 0.468355);
  AddElementByWeightFraction(52, 0.531645);

  AddMaterial("G4_CADMIUM_TUNGSTATE", 7.9, 0, 468.3, 3);
  AddElementByWeightFraction( 8, 0.177644);
  AddElementByWeightFraction(48, 0.312027);
  AddElementByWeightFraction(74, 0.510329);

  AddMaterial("G4_CALCIUM_CARBONATE", 2.8, 0, 136.4, 3);
  AddElementByWeightFraction( 6, 0.120003);
  AddElementByWeightFraction( 8, 0.479554);
  AddElementByWeightFraction(20, 0.400443);

  AddMaterial("G4_CALCIUM_FLUORIDE", 3.18, 0, 166., 2);
  AddElementByWeightFraction( 9, 0.486659);
  AddElementByWeightFraction(20, 0.513341);

  AddMaterial("G4_CALCIUM_OXIDE", 3.3, 0, 176.1, 2);
  AddElementByWeightFraction( 8, 0.285299);
  AddElementByWeightFraction(20, 0.714701);

  AddMaterial("G4_CALCIUM_SULFATE", 2.96, 0, 152.3, 3);
  AddElementByWeightFraction( 8, 0.470095);
  AddElementByWeightFraction(16, 0.235497);
  AddElementByWeightFraction(20, 0.294408);

  AddMaterial("G4_CALCIUM_TUNGSTATE", 6.062, 0, 395., 3);
  AddElementByWeightFraction( 8, 0.22227 );
  AddElementByWeightFraction(20, 0.139202);
  AddElementByWeightFraction(74, 0.638529);

  AddMaterial("G4_CARBON_DIOXIDE", 0.00184212, 0, 85., 2, kStateGas);
  AddElementByWeightFraction( 6, 0.272916);
  AddElementByWeightFraction( 8, 0.727084);
  chFormulas[nMaterials-1] = "CO_2";

  AddMaterial("G4_CARBON_TETRACHLORIDE", 1.594, 0, 166.3, 2);
  AddElementByWeightFraction( 6, 0.078083);
  AddElementByWeightFraction(17, 0.921917);

  AddMaterial("G4_CELLULOSE_CELLOPHANE", 1.42, 0, 77.6, 3);
  AddElementByWeightFraction( 1, 0.062162);
  AddElementByWeightFraction( 6, 0.444462);
  AddElementByWeightFraction( 8, 0.493376);

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
  AddElementByWeightFraction( 9, 0.125069);
  AddElementByWeightFraction(55, 0.874931);

  AddMaterial("G4_CESIUM_IODIDE", 4.51, 0, 553.1, 2);
  AddElementByWeightFraction(53, 0.488451);
  AddElementByWeightFraction(55, 0.511549);

  AddMaterial("G4_CHLOROBENZENE", 1.1058, 0, 89.1, 3);
  AddElementByWeightFraction( 1, 0.044772);
  AddElementByWeightFraction( 6, 0.640254);
  AddElementByWeightFraction(17, 0.314974);

  AddMaterial("G4_CHLOROFORM", 1.4832, 0, 156., 3);
  AddElementByWeightFraction( 1, 0.008443);
  AddElementByWeightFraction( 6, 0.100613);
  AddElementByWeightFraction(17, 0.890944);

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
  AddElementByWeightFraction( 1, 0.143711);
  AddElementByWeightFraction( 6, 0.856289);

  AddMaterial("G4_1,2-DICHLOROBENZENE", 1.3048, 0, 106.5, 3);
  AddElementByWeightFraction( 1, 0.027425);
  AddElementByWeightFraction( 6, 0.490233);
  AddElementByWeightFraction(17, 0.482342);

  AddMaterial("G4_DICHLORODIETHYL_ETHER", 1.2199, 0, 103.3, 4);
  AddElementByWeightFraction( 1, 0.056381);
  AddElementByWeightFraction( 6, 0.335942);
  AddElementByWeightFraction( 8, 0.111874);
  AddElementByWeightFraction(17, 0.495802);

  AddMaterial("G4_1,2-DICHLOROETHANE", 1.2351, 0, 111.9, 3);
  AddElementByWeightFraction( 1, 0.04074 );
  AddElementByWeightFraction( 6, 0.242746);
  AddElementByWeightFraction(17, 0.716515);

  AddMaterial("G4_DIETHYL_ETHER", 0.71378, 0, 60., 3);
  AddElementByWeightFraction( 1, 0.135978);
  AddElementByWeightFraction( 6, 0.648171);
  AddElementByWeightFraction( 8, 0.215851);

  AddMaterial("G4_N,N-DIMETHYL_FORMAMIDE", 0.9487, 0, 66.6, 4);
  AddElementByWeightFraction( 1, 0.096523);
  AddElementByWeightFraction( 6, 0.492965);
  AddElementByWeightFraction( 7, 0.191625);
  AddElementByWeightFraction( 8, 0.218887);

  AddMaterial("G4_DIMETHYL_SULFOXIDE", 1.1014, 0, 98.6, 4);
  AddElementByWeightFraction( 1, 0.077403);
  AddElementByWeightFraction( 6, 0.307467);
  AddElementByWeightFraction( 8, 0.204782);
  AddElementByWeightFraction(16, 0.410348);

  AddMaterial("G4_ETHANE", 0.00125324, 0, 45.4, 2, kStateGas);
  AddElementByWeightFraction( 1, 0.201115);
  AddElementByWeightFraction( 6, 0.798885);

  AddMaterial("G4_ETHYL_ALCOHOL", 0.7893, 0, 62.9, 3);
  AddElementByWeightFraction( 1, 0.131269);
  AddElementByWeightFraction( 6, 0.521438);
  AddElementByWeightFraction( 8, 0.347294);

  AddMaterial("G4_ETHYL_CELLULOSE", 1.13, 0, 69.3, 3);
  AddElementByWeightFraction( 1, 0.090027);
  AddElementByWeightFraction( 6, 0.585182);
  AddElementByWeightFraction( 8, 0.324791);

  AddMaterial("G4_ETHYLENE", 0.00117497, 0, 50.7, 2, kStateGas);
  AddElementByWeightFraction( 1, 0.143711);
  AddElementByWeightFraction( 6, 0.856289);

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
  AddElementByWeightFraction( 8, 0.300567);
  AddElementByWeightFraction(26, 0.699433);

  AddMaterial("G4_FERROBORIDE", 7.15, 0, 261., 2);
  AddElementByWeightFraction( 5, 0.162174);
  AddElementByWeightFraction(26, 0.837826);

  AddMaterial("G4_FERROUS_OXIDE", 5.7, 0, 248.6, 2);
  AddElementByWeightFraction( 8, 0.222689);
  AddElementByWeightFraction(26, 0.777311);

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
  AddElementByWeightFraction( 6, 0.080659);
  AddElementByWeightFraction( 9, 0.382749);
  AddElementByWeightFraction(35, 0.536592);

  AddMaterial("G4_FREON-13I1", 1.8, 0, 293.5, 3);
  AddElementByWeightFraction( 6, 0.061309);
  AddElementByWeightFraction( 9, 0.290924);
  AddElementByWeightFraction(53, 0.647767);

  AddMaterial("G4_GADOLINIUM_OXYSULFIDE", 7.44, 0, 493.3, 3);
  AddElementByWeightFraction( 8, 0.084528);
  AddElementByWeightFraction(16, 0.08469 );
  AddElementByWeightFraction(64, 0.830782);

  AddMaterial("G4_GALLIUM_ARSENIDE", 5.31, 0, 384.9, 2);
  AddElementByWeightFraction(31, 0.482019);
  AddElementByWeightFraction(33, 0.517981);

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

  AddMaterial("G4_GLUCOSE", 1.54, 0, 77.2, 3);
  AddElementByWeightFraction( 1, 0.071204);
  AddElementByWeightFraction( 6, 0.363652);
  AddElementByWeightFraction( 8, 0.565144);

  AddMaterial("G4_GLUTAMINE", 1.46, 0, 73.3, 4);
  AddElementByWeightFraction( 1, 0.068965);
  AddElementByWeightFraction( 6, 0.410926);
  AddElementByWeightFraction( 7, 0.191681);
  AddElementByWeightFraction( 8, 0.328427);

  AddMaterial("G4_GLYCEROL", 1.2613, 0, 72.6, 3);
  AddElementByWeightFraction( 1, 0.087554);
  AddElementByWeightFraction( 6, 0.391262);
  AddElementByWeightFraction( 8, 0.521185);

  AddMaterial("G4_GUANINE", 2.2/*1.58*/, 0, 75. ,4);
  AddElementByAtomCount("H",5 );
  AddElementByAtomCount("C",5 );
  AddElementByAtomCount("N",5 );
  AddElementByAtomCount("O",1 );

  AddMaterial("G4_GYPSUM", 2.32, 0, 129.7, 4);
  AddElementByWeightFraction( 1, 0.023416);
  AddElementByWeightFraction( 8, 0.557572);
  AddElementByWeightFraction(16, 0.186215);
  AddElementByWeightFraction(20, 0.232797);

  AddMaterial("G4_N-HEPTANE", 0.68376, 0, 54.4, 2);
  AddElementByWeightFraction( 1, 0.160937);
  AddElementByWeightFraction( 6, 0.839063);

  AddMaterial("G4_N-HEXANE", 0.6603, 0, 54., 2);
  AddElementByWeightFraction( 1, 0.163741);
  AddElementByWeightFraction( 6, 0.836259);

  AddMaterial("G4_KAPTON", 1.42, 0, 79.6, 4);
  AddElementByWeightFraction( 1, 0.026362);
  AddElementByWeightFraction( 6, 0.691133);
  AddElementByWeightFraction( 7, 0.07327 );
  AddElementByWeightFraction( 8, 0.209235);

  AddMaterial("G4_LANTHANUM_OXYBROMIDE", 6.28, 0, 439.7, 3);
  AddElementByWeightFraction( 8, 0.068138);
  AddElementByWeightFraction(35, 0.340294);
  AddElementByWeightFraction(57, 0.591568);

  AddMaterial("G4_LANTHANUM_OXYSULFIDE", 5.86, 0, 421.2, 3);
  AddElementByWeightFraction( 8, 0.0936  );
  AddElementByWeightFraction(16, 0.093778);
  AddElementByWeightFraction(57, 0.812622);

  AddMaterial("G4_LEAD_OXIDE", 9.53, 0, 766.7, 2);
  AddElementByWeightFraction( 8, 0.071682);
  AddElementByWeightFraction(82, 0.928318);

  AddMaterial("G4_LITHIUM_AMIDE", 1.178, 0, 55.5, 3);
  AddElementByWeightFraction( 1, 0.087783);
  AddElementByWeightFraction( 3, 0.302262);
  AddElementByWeightFraction( 7, 0.609955);

  AddMaterial("G4_LITHIUM_CARBONATE", 2.11, 0, 87.9, 3);
  AddElementByWeightFraction( 3, 0.187871);
  AddElementByWeightFraction( 6, 0.16255 );
  AddElementByWeightFraction( 8, 0.649579);

  AddMaterial("G4_LITHIUM_FLUORIDE", 2.635, 0, 94., 2);
  AddElementByWeightFraction( 3, 0.267585);
  AddElementByWeightFraction( 9, 0.732415);

  AddMaterial("G4_LITHIUM_HYDRIDE", 0.82, 0, 36.5, 2);
  AddElementByWeightFraction( 1, 0.126797);
  AddElementByWeightFraction( 3, 0.873203);

  AddMaterial("G4_LITHIUM_IODIDE", 3.494, 0, 485.1, 2);
  AddElementByWeightFraction( 3, 0.051858);
  AddElementByWeightFraction(53, 0.948142);

  AddMaterial("G4_LITHIUM_OXIDE", 2.013, 0, 73.6, 2);
  AddElementByWeightFraction( 3, 0.46457);
  AddElementByWeightFraction( 8, 0.53543);

  AddMaterial("G4_LITHIUM_TETRABORATE", 2.44, 0, 94.6, 3);
  AddElementByWeightFraction( 3, 0.082085);
  AddElementByWeightFraction( 5, 0.25568 );
  AddElementByWeightFraction( 8, 0.662235);

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
  AddElementByWeightFraction( 6, 0.142455);
  AddElementByWeightFraction( 8, 0.569278);
  AddElementByWeightFraction(12, 0.288267);

  AddMaterial("G4_MAGNESIUM_FLUORIDE", 3.0, 0, 134.3, 2);
  AddElementByWeightFraction( 9, 0.609883);
  AddElementByWeightFraction(12, 0.390117);

  AddMaterial("G4_MAGNESIUM_OXIDE", 3.58, 0, 143.8, 2);
  AddElementByWeightFraction( 8, 0.396964);
  AddElementByWeightFraction(12, 0.603036);

  AddMaterial("G4_MAGNESIUM_TETRABORATE", 2.53, 0, 108.3, 3);
  AddElementByWeightFraction( 5, 0.240837);
  AddElementByWeightFraction( 8, 0.62379);
  AddElementByWeightFraction(12, 0.135373);

  AddMaterial("G4_MERCURIC_IODIDE", 6.36, 0, 684.5, 2);
  AddElementByWeightFraction(53, 0.55856);
  AddElementByWeightFraction(80, 0.44144);

  AddMaterial("G4_METHANE", 0.000667151, 0, 41.7, 2, kStateGas);
  AddElementByWeightFraction( 1, 0.251306);
  AddElementByWeightFraction( 6, 0.748694);

  AddMaterial("G4_METHANOL", 0.7914, 0, 67.6, 3);
  AddElementByWeightFraction( 1, 0.125822);
  AddElementByWeightFraction( 6, 0.374852);
  AddElementByWeightFraction( 8, 0.499326);

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
  AddElementByWeightFraction( 1, 0.062909);
  AddElementByWeightFraction( 6, 0.937091);

  AddMaterial("G4_NITROBENZENE", 1.19867, 0, 75.8, 4);
  AddElementByWeightFraction( 1, 0.040935);
  AddElementByWeightFraction( 6, 0.585374);
  AddElementByWeightFraction( 7, 0.113773);
  AddElementByWeightFraction( 8, 0.259918);

  AddMaterial("G4_NITROUS_OXIDE", 0.00183094, 0, 84.9, 2, kStateGas);
  AddElementByWeightFraction( 7, 0.636483);
  AddElementByWeightFraction( 8, 0.363517);

  AddMaterial("G4_NYLON-8062", 1.08, 0, 64.3, 4);
  AddElementByWeightFraction( 1, 0.103509);
  AddElementByWeightFraction( 6, 0.648415);
  AddElementByWeightFraction( 7, 0.099536);
  AddElementByWeightFraction( 8, 0.148539);

  AddMaterial("G4_NYLON-6-6", 1.14, 0, 63.9, 4);
  AddElementByWeightFraction( 1, 0.097976);
  AddElementByWeightFraction( 6, 0.636856);
  AddElementByWeightFraction( 7, 0.123779);
  AddElementByWeightFraction( 8, 0.141389);

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
  AddElementByWeightFraction( 1, 0.158821);
  AddElementByWeightFraction( 6, 0.841179);

  AddMaterial("G4_PARAFFIN", 0.93, 0, 55.9, 2);
  AddElementByWeightFraction( 1, 0.148605);
  AddElementByWeightFraction( 6, 0.851395);

  AddMaterial("G4_N-PENTANE", 0.6262, 0, 53.6, 2);
  AddElementByWeightFraction( 1, 0.167635);
  AddElementByWeightFraction (6, 0.832365);

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
  AddElementByWeightFraction( 1, 0.085);
  AddElementByWeightFraction( 6, 0.915);

  AddMaterial("G4_PLUTONIUM_DIOXIDE", 11.46, 0, 746.5, 2);
  AddElementByWeightFraction( 8, 0.118055);
  AddElementByWeightFraction(94, 0.881945);

  AddMaterial("G4_POLYACRYLONITRILE", 1.17, 0, 69.6, 3);
  AddElementByWeightFraction( 1, 0.056983);
  AddElementByWeightFraction( 6, 0.679056);
  AddElementByWeightFraction( 7, 0.263962);

  AddMaterial("G4_POLYCARBONATE", 1.2, 0, 73.1, 3);
  AddElementByWeightFraction( 1, 0.055491);
  AddElementByWeightFraction( 6, 0.755751);
  AddElementByWeightFraction( 8, 0.188758);

  AddMaterial("G4_POLYCHLOROSTYRENE", 1.3, 0, 81.7, 3);
  AddElementByWeightFraction( 1, 0.061869);
  AddElementByWeightFraction( 6, 0.696325);
  AddElementByWeightFraction(17, 0.241806);

  AddMaterial("G4_POLYETHYLENE", 0.94, 0, 57.4, 2);
  AddElementByWeightFraction( 1, 0.143711);
  AddElementByWeightFraction( 6, 0.856289);
  chFormulas[nMaterials-1] = "(C_2H_4)_N-Polyethylene";

  AddMaterial("G4_MYLAR", 1.4, 0, 78.7, 3);
  AddElementByWeightFraction( 1, 0.041959);
  AddElementByWeightFraction( 6, 0.625017);
  AddElementByWeightFraction( 8, 0.333025);

  AddMaterial("G4_PLEXIGLASS", 1.19, 0, 74., 3);
  AddElementByWeightFraction( 1, 0.080538);
  AddElementByWeightFraction( 6, 0.599848);
  AddElementByWeightFraction( 8, 0.319614);

  AddMaterial("G4_POLYOXYMETHYLENE", 1.425 ,0, 77.4, 3);
  AddElementByWeightFraction( 1, 0.067135);
  AddElementByWeightFraction( 6, 0.400017);
  AddElementByWeightFraction( 8, 0.532848);

  AddMaterial("G4_POLYPROPYLENE", 0.9, 0, 56.5, 2);
  AddElementByWeightFraction( 1, 0.143711);
  AddElementByWeightFraction( 6, 0.856289);
  chFormulas[nMaterials-1] = "(C_2H_4)_N-Polypropylene";

  AddMaterial("G4_POLYSTYRENE", 1.06, 0, 68.7, 2);
  AddElementByWeightFraction( 1, 0.077418);
  AddElementByWeightFraction( 6, 0.922582);

  AddMaterial("G4_TEFLON", 2.2, 0, 99.1, 2);
  AddElementByWeightFraction( 6, 0.240183);
  AddElementByWeightFraction( 9, 0.759817);

  AddMaterial("G4_POLYTRIFLUOROCHLOROETHYLENE", 2.1, 0, 120.7, 3);
  // correct chemical name Polychlorotrifluoroethylene [CF2CClF]n, IvantchenkoA.
  AddElementByWeightFraction( 6, 0.20625 );
  AddElementByWeightFraction( 9, 0.489354);
  AddElementByWeightFraction(17, 0.304395);

  AddMaterial("G4_POLYVINYL_ACETATE", 1.19, 0, 73.7, 3);
  AddElementByWeightFraction( 1, 0.070245);
  AddElementByWeightFraction( 6, 0.558066);
  AddElementByWeightFraction( 8, 0.371689);

  AddMaterial("G4_POLYVINYL_ALCOHOL", 1.3, 0, 69.7, 3);
  AddElementByWeightFraction( 1, 0.091517);
  AddElementByWeightFraction( 6, 0.545298);
  AddElementByWeightFraction( 8, 0.363185);

  AddMaterial("G4_POLYVINYL_BUTYRAL", 1.12, 0, 67.2, 3);
  AddElementByWeightFraction( 1, 0.092802);
  AddElementByWeightFraction( 6, 0.680561);
  AddElementByWeightFraction( 8, 0.226637);

  AddMaterial("G4_POLYVINYL_CHLORIDE", 1.3, 0, 108.2, 3);
  AddElementByWeightFraction( 1, 0.04838);
  AddElementByWeightFraction( 6, 0.38436);
  AddElementByWeightFraction(17, 0.56726);

  AddMaterial("G4_POLYVINYLIDENE_CHLORIDE", 1.7, 0, 134.3, 3);
  AddElementByWeightFraction( 1, 0.020793);
  AddElementByWeightFraction( 6, 0.247793);
  AddElementByWeightFraction(17, 0.731413);

  AddMaterial("G4_POLYVINYLIDENE_FLUORIDE", 1.76, 0, 88.8, 3);
  AddElementByWeightFraction( 1, 0.03148 );
  AddElementByWeightFraction( 6, 0.375141);
  AddElementByWeightFraction( 9, 0.593379);

  AddMaterial("G4_POLYVINYL_PYRROLIDONE", 1.25, 0, 67.7, 4);
  AddElementByWeightFraction( 1, 0.081616);
  AddElementByWeightFraction( 6, 0.648407);
  AddElementByWeightFraction( 7, 0.126024);
  AddElementByWeightFraction( 8, 0.143953);

  AddMaterial("G4_POTASSIUM_IODIDE", 3.13, 0, 431.9, 2);
  AddElementByWeightFraction(19, 0.235528);
  AddElementByWeightFraction(53, 0.764472);

  AddMaterial("G4_POTASSIUM_OXIDE", 2.32, 0, 189.9, 2);
  AddElementByWeightFraction( 8, 0.169852);
  AddElementByWeightFraction(19, 0.830148);

  AddMaterial("G4_PROPANE", 0.00187939, 0, 47.1, 2, kStateGas);
  AddElementByWeightFraction( 1, 0.182855);
  AddElementByWeightFraction( 6, 0.817145);

  AddMaterial("G4_lPROPANE", 0.43, 0, 52., 2);
  AddElementByWeightFraction( 1, 0.182855);
  AddElementByWeightFraction( 6, 0.817145);

  AddMaterial("G4_N-PROPYL_ALCOHOL", 0.8035, 0, 61.1, 3);
  AddElementByWeightFraction( 1, 0.134173);
  AddElementByWeightFraction( 6, 0.599595);
  AddElementByWeightFraction( 8, 0.266232);

  AddMaterial("G4_PYRIDINE", 0.9819, 0, 66.2, 3);
  AddElementByWeightFraction( 1, 0.06371 );
  AddElementByWeightFraction( 6, 0.759217);
  AddElementByWeightFraction( 7, 0.177073);

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
  AddElementByWeightFraction( 8, 0.532565);
  AddElementByWeightFraction(14, 0.467435);
  chFormulas[nMaterials-1] = "SiO_2";

  AddMaterial("G4_SILVER_BROMIDE", 6.473, 0, 486.6, 2);
  AddElementByWeightFraction(35, 0.425537);
  AddElementByWeightFraction(47, 0.574463);

  AddMaterial("G4_SILVER_CHLORIDE", 5.56, 0, 398.4, 2);
  AddElementByWeightFraction(17, 0.247368);
  AddElementByWeightFraction(47, 0.752632);

  AddMaterial("G4_SILVER_HALIDES", 6.47, 0, 487.1, 3);
  AddElementByWeightFraction(35, 0.422895);
  AddElementByWeightFraction(47, 0.573748);
  AddElementByWeightFraction(53, 0.003357);

  AddMaterial("G4_SILVER_IODIDE", 6.01, 0, 543.5, 2);
  AddElementByWeightFraction(47, 0.459458);
  AddElementByWeightFraction(53, 0.540542);

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
  AddElementByWeightFraction( 6, 0.113323);
  AddElementByWeightFraction( 8, 0.452861);
  AddElementByWeightFraction(11, 0.433815);

  AddMaterial("G4_SODIUM_IODIDE", 3.667, 0, 452., 2);
  AddElementByWeightFraction(11, 0.153373);
  AddElementByWeightFraction(53, 0.846627);

  AddMaterial("G4_SODIUM_MONOXIDE", 2.27, 0, 148.8, 2);
  AddElementByWeightFraction( 8, 0.258143);
  AddElementByWeightFraction(11, 0.741857);

  AddMaterial("G4_SODIUM_NITRATE", 2.261, 0, 114.6, 3);
  AddElementByWeightFraction( 7, 0.164795);
  AddElementByWeightFraction( 8, 0.56472 );
  AddElementByWeightFraction(11, 0.270485);

  AddMaterial("G4_STILBENE", 0.9707, 0, 67.7, 2);
  AddElementByWeightFraction( 1, 0.067101);
  AddElementByWeightFraction( 6, 0.932899);

  AddMaterial("G4_SUCROSE", 1.5805, 0, 77.5, 3);
  AddElementByWeightFraction( 1, 0.064779);
  AddElementByWeightFraction( 6, 0.42107);
  AddElementByWeightFraction( 8, 0.514151);

  AddMaterial("G4_TERPHENYL", 1.234, 0, 71.7, 2);
  AddElementByWeightFraction( 1, 0.044543);
  AddElementByWeightFraction( 6, 0.955457);

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
  AddElementByWeightFraction( 6, 0.144856);
  AddElementByWeightFraction(17, 0.855144);

  AddMaterial("G4_THALLIUM_CHLORIDE", 7.004, 0, 690.3, 2);
  AddElementByWeightFraction(17, 0.147822);
  AddElementByWeightFraction(81, 0.852178);

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
  AddElementByWeightFraction( 8, 0.400592);
  AddElementByWeightFraction(22, 0.599408);

  AddMaterial("G4_TOLUENE", 0.8669, 0, 62.5, 2);
  AddElementByWeightFraction( 1, 0.08751);
  AddElementByWeightFraction( 6, 0.91249);

  AddMaterial("G4_TRICHLOROETHYLENE", 1.46, 0, 148.1, 3);
  AddElementByWeightFraction( 1, 0.007671);
  AddElementByWeightFraction( 6, 0.182831);
  AddElementByWeightFraction(17, 0.809498);

  AddMaterial("G4_TRIETHYL_PHOSPHATE", 1.07, 0, 81.2, 4);
  AddElementByWeightFraction( 1, 0.082998);
  AddElementByWeightFraction( 6, 0.395628);
  AddElementByWeightFraction( 8, 0.351334);
  AddElementByWeightFraction(15, 0.17004 );

  AddMaterial("G4_TUNGSTEN_HEXAFLUORIDE", 2.4, 0, 354.4, 2);
  AddElementByWeightFraction( 9, 0.382723);
  AddElementByWeightFraction(74, 0.617277);

  AddMaterial("G4_URANIUM_DICARBIDE", 11.28, 0, 752., 2);
  AddElementByWeightFraction( 6, 0.091669);
  AddElementByWeightFraction(92, 0.908331);

  AddMaterial("G4_URANIUM_MONOCARBIDE", 13.63, 0, 862., 2);
  AddElementByWeightFraction( 6, 0.048036);
  AddElementByWeightFraction(92, 0.951964);

  AddMaterial("G4_URANIUM_OXIDE", 10.96, 0, 720.6, 2);
  AddElementByWeightFraction( 8, 0.118502);
  AddElementByWeightFraction(92, 0.881498);

  AddMaterial("G4_UREA", 1.323, 0, 72.8, 4);
  AddElementByWeightFraction( 1, 0.067131);
  AddElementByWeightFraction( 6, 0.199999);
  AddElementByWeightFraction( 7, 0.466459);
  AddElementByWeightFraction( 8, 0.266411);

  AddMaterial("G4_VALINE", 1.23, 0, 67.7, 4);
  AddElementByWeightFraction( 1, 0.094641);
  AddElementByWeightFraction( 6, 0.512645);
  AddElementByWeightFraction( 7, 0.119565);
  AddElementByWeightFraction( 8, 0.27315 );

  AddMaterial("G4_VITON", 1.8, 0, 98.6, 3);
  AddElementByWeightFraction( 1, 0.009417);
  AddElementByWeightFraction( 6, 0.280555);
  AddElementByWeightFraction( 9, 0.710028);

  AddMaterial("G4_WATER", 1.0,0, 78., 2);
  AddElementByAtomCount("H", 2);
  AddElementByAtomCount("O", 1);
  chFormulas[nMaterials-1] = "H_2O";

  AddMaterial("G4_WATER_VAPOR", 0.000756182, 0, 71.6, 2, kStateGas);
  AddElementByAtomCount("H", 2);
  AddElementByAtomCount("O", 1);
  chFormulas[nMaterials-1] = "H_2O-Gas";

  AddMaterial("G4_XYLENE", 0.87, 0, 61.8, 2);
  AddElementByWeightFraction( 1, 0.094935);
  AddElementByWeightFraction( 6, 0.905065);

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
  AddMaterial("G4_lKr", 2.418 , 36, 352. , 1, kStateLiquid, false);
  AddMaterial("G4_lXe", 2.953 , 54, 482. , 1, kStateLiquid, false);

  AddMaterial("G4_PbWO4", 8.28, 0, 0.0, 3);
  AddElementByAtomCount("O" , 4);
  AddElementByAtomCount("Pb", 1);
  AddElementByAtomCount("W" , 1);

  G4double density = universe_mean_density*cm3/g;
  AddMaterial("G4_Galactic", density, 1, 21.8, 1, kStateGas);
  AddGas("G4_Galactic",2.73*kelvin, 3.e-18*pascal);

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
  AddElementByAtomCount("Cl", 1);

  nSpace = nMaterials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialBuilder::BioChemicalMaterials()
{
  AddMaterial("G4_CYTOSINE", 1.55, 0, 72., 4);
  AddElementByAtomCount("H", 5);
  AddElementByAtomCount("C", 4);
  AddElementByAtomCount("N", 3);
  AddElementByAtomCount("O", 1);

  AddMaterial("G4_THYMINE", 1.23, 0, 72., 4);
  AddElementByAtomCount("H", 6);
  AddElementByAtomCount("C", 5);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 2);

  AddMaterial("G4_URACIL", 1.32, 0, 72., 4);
  AddElementByAtomCount("H", 4);
  AddElementByAtomCount("C", 4);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 2);

  // DNA_Nucleobase (Nucleobase-1H)
  AddMaterial("G4_DNA_ADENINE", 1, 0, 72., 3);
  AddElementByAtomCount("H",4 );
  AddElementByAtomCount("C",5 );
  AddElementByAtomCount("N",5 );

  AddMaterial("G4_DNA_GUANINE", 1, 0, 72. ,4);
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

  // DNA_Nucleoside (Nucleoside-3H)
  AddMaterial("G4_DNA_ADENOSINE", 1, 0, 72., 4);
  AddElementByAtomCount("H", 10);
  AddElementByAtomCount("C", 10);
  AddElementByAtomCount("N", 5);
  AddElementByAtomCount("O", 4);

  AddMaterial("G4_DNA_GUANOSINE", 1, 0, 72. ,4);
  AddElementByAtomCount("H", 10);
  AddElementByAtomCount("C", 10);
  AddElementByAtomCount("N", 5);
  AddElementByAtomCount("O", 5);

  AddMaterial("G4_DNA_CYTIDINE", 1, 0, 72., 4);
  AddElementByAtomCount("H", 10);
  AddElementByAtomCount("C", 9);
  AddElementByAtomCount("N", 3);
  AddElementByAtomCount("O", 5);

  AddMaterial("G4_DNA_URIDINE", 1, 0, 72., 4);
  AddElementByAtomCount("H", 9);
  AddElementByAtomCount("C", 9);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 6);

  AddMaterial("G4_DNA_METHYLURIDINE", 1, 0, 72., 4);
  AddElementByAtomCount("H", 11);
  AddElementByAtomCount("C", 10);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 6);

  AddMaterial("G4_DNA_MONOPHOSPHATE", 1, 0, 72., 2);
  AddElementByAtomCount("P", 1);
  AddElementByAtomCount("O", 3);

  AddMaterial("G4_DNA_A", 1, 0, 72., 5);  //Adenine base
  AddElementByAtomCount("H", 10);
  AddElementByAtomCount("C", 10);
  AddElementByAtomCount("N", 5);
  AddElementByAtomCount("O", 7);
  AddElementByAtomCount("P", 1);

  AddMaterial("G4_DNA_G", 1, 0, 72. ,5); //Guanine base 
  AddElementByAtomCount("H", 10);
  AddElementByAtomCount("C", 10);
  AddElementByAtomCount("N", 5);
  AddElementByAtomCount("O", 8);
  AddElementByAtomCount("P", 1);

  AddMaterial("G4_DNA_C", 1, 0, 72., 5); // Cytosine base
  AddElementByAtomCount("H", 10);
  AddElementByAtomCount("C", 9);
  AddElementByAtomCount("N", 3);
  AddElementByAtomCount("O", 8);
  AddElementByAtomCount("P", 1);

  AddMaterial("G4_DNA_U", 1, 0, 72., 5); // Uracil base
  AddElementByAtomCount("H", 9);
  AddElementByAtomCount("C", 9);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 9);
  AddElementByAtomCount("P", 1);

  AddMaterial("G4_DNA_MU", 1, 0, 72., 5);  // MethaUracil base
  AddElementByAtomCount("H", 11);
  AddElementByAtomCount("C", 10);
  AddElementByAtomCount("N", 2);
  AddElementByAtomCount("O", 9);
  AddElementByAtomCount("P", 1);
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




