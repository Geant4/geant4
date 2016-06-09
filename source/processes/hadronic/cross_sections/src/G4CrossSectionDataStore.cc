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
// $Id: G4CrossSectionDataStore.cc,v 1.15.2.2 2009/08/13 09:16:18 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-03 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4CrossSectionDataStore
//
//
// Modifications:
// 23.01.2009 V.Ivanchenko add destruction of data sets
//
//

#include "G4CrossSectionDataStore.hh"
#include "G4HadronicException.hh"
#include "G4StableIsotopes.hh"
#include "G4HadTmpUtil.hh"
#include "Randomize.hh"
#include "G4Nucleus.hh" 

G4CrossSectionDataStore::G4CrossSectionDataStore() :
  NDataSetList(0), verboseLevel(0)
{}

G4CrossSectionDataStore::~G4CrossSectionDataStore()
{}

G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* aParticle,
                                         const G4Element* anElement,
					 G4double aTemperature)
{
  if (NDataSetList == 0) 
  {
    throw G4HadronicException(__FILE__, __LINE__, 
     "G4CrossSectionDataStore: no data sets registered");
    return DBL_MIN;
  }
  for (G4int i = NDataSetList-1; i >= 0; i--) {
    if (DataSetList[i]->IsApplicable(aParticle, anElement))
      return DataSetList[i]->GetCrossSection(aParticle,anElement,aTemperature);
  }
  throw G4HadronicException(__FILE__, __LINE__, 
                      "G4CrossSectionDataStore: no applicable data set found "
                      "for particle/element");
  return DBL_MIN;
}

G4VCrossSectionDataSet*
G4CrossSectionDataStore::whichDataSetInCharge(const G4DynamicParticle* aParticle,
                                              const G4Element* anElement)
{
  if (NDataSetList == 0) 
  {
    throw G4HadronicException(__FILE__, __LINE__, 
     "G4CrossSectionDataStore: no data sets registered");
    return 0;
  }
  for (G4int i = NDataSetList-1; i >= 0; i--) {
    if (DataSetList[i]->IsApplicable(aParticle, anElement) )
    {
      return DataSetList[i];
    }
  }
  throw G4HadronicException(__FILE__, __LINE__, 
                      "G4CrossSectionDataStore: no applicable data set found "
                      "for particle/element");
  return 0;
}


G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* aParticle,
                                         const G4Isotope* anIsotope,
					 G4double aTemperature)
{
  if (NDataSetList == 0) 
  {
    throw G4HadronicException(__FILE__, __LINE__, 
     "G4CrossSectionDataStore: no data sets registered");
    return DBL_MIN;
  }
  for (G4int i = NDataSetList-1; i >= 0; i--) {
    if (DataSetList[i]->IsZAApplicable(aParticle, anIsotope->GetZ(), anIsotope->GetN()))
      return DataSetList[i]->GetIsoCrossSection(aParticle,anIsotope,aTemperature);
  }
  throw G4HadronicException(__FILE__, __LINE__, 
                      "G4CrossSectionDataStore: no applicable data set found "
                      "for particle/element");
  return DBL_MIN;
}


G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* aParticle,
                                         G4double Z, G4double A,
					 G4double aTemperature)
{
  if (NDataSetList == 0) 
  {
    throw G4HadronicException(__FILE__, __LINE__, 
     "G4CrossSectionDataStore: no data sets registered");
    return DBL_MIN;
  }
  for (G4int i = NDataSetList-1; i >= 0; i--) {
      if (DataSetList[i]->IsZAApplicable(aParticle, Z, A))
      return DataSetList[i]->GetIsoZACrossSection(aParticle,Z,A,aTemperature);
  }
  throw G4HadronicException(__FILE__, __LINE__, 
                      "G4CrossSectionDataStore: no applicable data set found "
                      "for particle/element");
  return DBL_MIN;
}


G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* aParticle,
                                         const G4Material* aMaterial)
{
  G4double sigma(0);
  G4double aTemp = aMaterial->GetTemperature();
  G4int nElements = aMaterial->GetNumberOfElements();
  const G4double* theAtomsPerVolumeVector = aMaterial->GetVecNbOfAtomsPerVolume();
  G4double xSection(0);

  for(G4int i=0; i<nElements; i++) {
    xSection = GetCrossSection( aParticle, (*aMaterial->GetElementVector())[i], aTemp);
    sigma += theAtomsPerVolumeVector[i] * xSection;
  }

  return sigma;
}


G4Element* G4CrossSectionDataStore::SampleZandA(const G4DynamicParticle* particle, 
						const G4Material* aMaterial,
						G4Nucleus& target)
{
  G4double aTemp = aMaterial->GetTemperature();
  const G4int nElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  G4Element* anElement = (*theElementVector)[0];
  G4VCrossSectionDataSet* inCharge;
  G4int i;

  // compounds
  if(1 < nElements) {
    G4double* xsec = new G4double [nElements];
    const G4double* theAtomsPerVolumeVector = aMaterial->GetVecNbOfAtomsPerVolume();
    G4double cross = 0.0;
    for(i=0; i<nElements; i++) {
      anElement= (*theElementVector)[i];
      inCharge = whichDataSetInCharge(particle, anElement);
      cross   += theAtomsPerVolumeVector[i]*
	inCharge->GetCrossSection(particle, anElement, aTemp);
      xsec[i]  = cross;
    }
    cross *= G4UniformRand();

    for(i=0; i<nElements; i++) {
      if( cross <=  xsec[i] ) {
	anElement = (*theElementVector)[i];
	break;
      }
    }
    delete [] xsec;
  }

  // element have been selected
  inCharge = whichDataSetInCharge(particle, anElement);
  G4double ZZ = anElement->GetZ();
  G4double AA;

  // Collect abundance weighted cross sections and A values for each isotope 
  // in each element
 
  const G4int nIsoPerElement = anElement->GetNumberOfIsotopes();

  // user-defined isotope abundances
  if (0 < nIsoPerElement) { 
    G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
    AA = G4double((*isoVector)[0]->GetN());
    if(1 < nIsoPerElement) {

      G4double* xsec  = new G4double [nIsoPerElement];
      G4double iso_xs = 0.0;
      G4double cross  = 0.0;

      G4double* abundVector = anElement->GetRelativeAbundanceVector();
      G4bool elementXS = false;
      for (i = 0; i<nIsoPerElement; i++) {
	if (inCharge->IsZAApplicable(particle, ZZ, G4double((*isoVector)[i]->GetN()))) {
	  iso_xs = inCharge->GetIsoCrossSection(particle, (*isoVector)[i], aTemp);
	} else if (elementXS == false) {
	  iso_xs = inCharge->GetCrossSection(particle, anElement, aTemp);
	  elementXS = true;
	}

	cross  += abundVector[i]*iso_xs;
	xsec[i] = cross;
      }
      cross *= G4UniformRand();
      for (i = 0; i<nIsoPerElement; i++) {
	if(cross <= xsec[i]) {
	  AA = G4double((*isoVector)[i]->GetN());
	  break;
	}
      }
      delete [] xsec;
    }
    // natural abundances
  } else { 

    G4StableIsotopes theDefaultIsotopes;  
    G4int Z = G4int(ZZ + 0.5);
    const G4int nIso = theDefaultIsotopes.GetNumberOfIsotopes(Z);
    G4int index = theDefaultIsotopes.GetFirstIsotope(Z);
    AA = G4double(theDefaultIsotopes.GetIsotopeNucleonCount(index));
    
    if(1 < nIso) {

      G4double* xsec  = new G4double [nIso];
      G4double iso_xs = 0.0;
      G4double cross  = 0.0;
      G4bool elementXS= false;

      for (i = 0; i<nIso; i++) {
        AA = G4double(theDefaultIsotopes.GetIsotopeNucleonCount(index+i));
	if (inCharge->IsZAApplicable(particle, ZZ, AA )) {
	  iso_xs = inCharge->GetIsoZACrossSection(particle, ZZ, AA, aTemp);
	} else if (elementXS == false) {
	  iso_xs = inCharge->GetCrossSection(particle, anElement, aTemp);
	  elementXS = true;
	}
	cross  += theDefaultIsotopes.GetAbundance(index+i)*iso_xs;
	xsec[i] = cross;
      }
      cross *= G4UniformRand();
      for (i = 0; i<nIso; i++) {
	if(cross <= xsec[i]) {
	  AA = G4double(theDefaultIsotopes.GetIsotopeNucleonCount(index+i));
	  break;
	}
      }
      delete [] xsec;
    }
  }
  //G4cout << "XS: " << particle->GetDefinition()->GetParticleName()
  //	 << " e(GeV)= " << particle->GetKineticEnergy()/GeV
  //	 << " in " << aMaterial->GetName()
  //	 << " ZZ= " << ZZ << " AA= " << AA << "  " << anElement->GetName() << G4endl;

  target.SetParameters(AA, ZZ);
  return anElement;
}

/*
G4Element* G4CrossSectionDataStore::SampleZandA(const G4DynamicParticle* particle, 
						const G4Material* aMaterial,
						G4Nucleus& target)
{
  static G4StableIsotopes theDefaultIsotopes;  // natural abundances and 
                                               // stable isotopes
  G4double aTemp = aMaterial->GetTemperature();
  G4int nElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  const G4double* theAtomsPerVolumeVector = aMaterial->GetVecNbOfAtomsPerVolume();
  std::vector<std::vector<G4double> > awicsPerElement;
  std::vector<std::vector<G4double> > AvaluesPerElement;
  G4Element* anElement;

  // Collect abundance weighted cross sections and A values for each isotope 
  // in each element
 
  for (G4int i = 0; i < nElements; i++) {
    anElement = (*theElementVector)[i];
    G4int nIsoPerElement = anElement->GetNumberOfIsotopes();
    std::vector<G4double> isoholder;
    std::vector<G4double> aholder;
    G4double iso_xs = DBL_MIN;

    if (nIsoPerElement) { // user-defined isotope abundances
      G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
      G4double* abundVector = anElement->GetRelativeAbundanceVector();
      G4VCrossSectionDataSet* inCharge = whichDataSetInCharge(particle, anElement);
      G4bool elementXS = false;
      for (G4int j = 0; j < nIsoPerElement; j++) {
        if (inCharge->IsZAApplicable(particle, (*isoVector)[j]->GetZ(), 
                                               (*isoVector)[j]->GetN() ) ) {
          iso_xs = inCharge->GetIsoCrossSection(particle, (*isoVector)[j], aTemp);
        } else if (elementXS == false) {
          iso_xs = inCharge->GetCrossSection(particle, anElement, aTemp);
          elementXS = true;
        }

        isoholder.push_back(abundVector[j]*iso_xs);
        aholder.push_back(G4double((*isoVector)[j]->GetN()));
      }

    } else { // natural abundances
      G4int ZZ = G4lrint(anElement->GetZ());
      nIsoPerElement = theDefaultIsotopes.GetNumberOfIsotopes(ZZ);
      G4int index = theDefaultIsotopes.GetFirstIsotope(ZZ);
      G4double AA;
      G4double abundance;

      G4VCrossSectionDataSet* inCharge = whichDataSetInCharge(particle, anElement);
      G4bool elementXS = false;

      for (G4int j = 0; j < nIsoPerElement; j++) {
        AA = G4double(theDefaultIsotopes.GetIsotopeNucleonCount(index+j));
        aholder.push_back(AA);
        if (inCharge->IsZAApplicable(particle, G4double(ZZ), AA) ) {
          iso_xs = inCharge->GetIsoZACrossSection(particle, G4double(ZZ), AA, aTemp);
        } else if (elementXS == false) {
          iso_xs = inCharge->GetCrossSection(particle, anElement, aTemp);
          elementXS = true;
        }
        abundance = theDefaultIsotopes.GetAbundance(index+j)/100.0;
        isoholder.push_back(abundance*iso_xs);
      }
    }

    awicsPerElement.push_back(isoholder);
    AvaluesPerElement.push_back(aholder);
  }

  // Calculate running sums for isotope selection

  G4double crossSectionTotal = 0;
  G4double xSectionPerElement;
  std::vector<G4double> runningSum;

  for (G4int i=0; i < nElements; i++) {
    xSectionPerElement = 0;
    for (G4int j=0; j < G4int(awicsPerElement[i].size()); j++)
                     xSectionPerElement += awicsPerElement[i][j];
    runningSum.push_back(theAtomsPerVolumeVector[i]*xSectionPerElement);
    crossSectionTotal += runningSum[i];
  }

  // Compare random number to running sum over element xc to choose Z

  // Initialize Z and A to first element and first isotope in case 
  // cross section is zero
   
  anElement = (*theElementVector)[0];
  G4double ZZ = anElement->GetZ();
  G4double AA = AvaluesPerElement[0][0];
  if (crossSectionTotal != 0.) {
    G4double random = G4UniformRand();
    for(G4int i=0; i < nElements; i++) {
      if(i!=0) runningSum[i] += runningSum[i-1];
      if(random <= runningSum[i]/crossSectionTotal) {
	anElement = (*theElementVector)[i];
        ZZ = anElement->GetZ();

        // Compare random number to running sum over isotope xc to choose A

        G4int nIso = awicsPerElement[i].size();
        G4double* running = new G4double[nIso];
        for (G4int j=0; j < nIso; j++) {
          running[j] = awicsPerElement[i][j];
          if(j!=0) running[j] += running[j-1];
        }

        G4double trial = G4UniformRand(); 
        for (G4int j=0; j < nIso; j++) {
          AA = AvaluesPerElement[i][j];
          if (trial <= running[j]/running[nIso-1]) break;
        }
        delete [] running;
        break;
      }
    }
  }

  //G4cout << "XS: " << particle->GetDefinition()->GetParticleName()
  //	 << " e(GeV)= " << particle->GetKineticEnergy()/GeV
  //	 << " in " << aMaterial->GetName()
  //	 << " ZZ= " << ZZ << " AA= " << AA << "  " << anElement->GetName() << G4endl;

  target.SetParameters(AA, ZZ);
  return anElement;
}

std::pair<G4double, G4double> 
G4CrossSectionDataStore::SelectRandomIsotope(const G4DynamicParticle* particle,
                                             const G4Material* aMaterial)
{
  static G4StableIsotopes theDefaultIsotopes;  // natural abundances and 
                                               // stable isotopes
  G4double aTemp = aMaterial->GetTemperature();
  G4int nElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  const G4double* theAtomsPerVolumeVector = aMaterial->GetVecNbOfAtomsPerVolume();
  std::vector<std::vector<G4double> > awicsPerElement;
  std::vector<std::vector<G4double> > AvaluesPerElement;
  G4Element* anElement;

  // Collect abundance weighted cross sections and A values for each isotope 
  // in each element
 
  for (G4int i = 0; i < nElements; i++) {
    anElement = (*theElementVector)[i];
    G4int nIsoPerElement = anElement->GetNumberOfIsotopes();
    std::vector<G4double> isoholder;
    std::vector<G4double> aholder;
    G4double iso_xs = DBL_MIN;

    if (nIsoPerElement) { // user-defined isotope abundances
      G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
      G4double* abundVector = anElement->GetRelativeAbundanceVector();
      G4VCrossSectionDataSet* inCharge = whichDataSetInCharge(particle, anElement);
      G4bool elementXS = false;
      for (G4int j = 0; j < nIsoPerElement; j++) {
        if (inCharge->IsZAApplicable(particle, (*isoVector)[j]->GetZ(), 
                                               (*isoVector)[j]->GetN() ) ) {
          iso_xs = inCharge->GetIsoCrossSection(particle, (*isoVector)[j], aTemp);
        } else if (elementXS == false) {
          iso_xs = inCharge->GetCrossSection(particle, anElement, aTemp);
          elementXS = true;
        }

        isoholder.push_back(abundVector[j]*iso_xs);
        aholder.push_back(G4double((*isoVector)[j]->GetN()));
      }

    } else { // natural abundances
      G4int ZZ = G4lrint(anElement->GetZ());
      nIsoPerElement = theDefaultIsotopes.GetNumberOfIsotopes(ZZ);
      G4int index = theDefaultIsotopes.GetFirstIsotope(ZZ);
      G4double AA;
      G4double abundance;

      G4VCrossSectionDataSet* inCharge = whichDataSetInCharge(particle, anElement);
      G4bool elementXS = false;

      for (G4int j = 0; j < nIsoPerElement; j++) {
        AA = G4double(theDefaultIsotopes.GetIsotopeNucleonCount(index+j));
        aholder.push_back(AA);
        if (inCharge->IsZAApplicable(particle, G4double(ZZ), AA) ) {
          iso_xs = inCharge->GetIsoZACrossSection(particle, G4double(ZZ), AA, aTemp);
        } else if (elementXS == false) {
          iso_xs = inCharge->GetCrossSection(particle, anElement, aTemp);
          elementXS = true;
        }
        abundance = theDefaultIsotopes.GetAbundance(index+j)/100.0;
        isoholder.push_back(abundance*iso_xs);
      }
    }

    awicsPerElement.push_back(isoholder);
    AvaluesPerElement.push_back(aholder);
  }

  // Calculate running sums for isotope selection

  G4double crossSectionTotal = 0;
  G4double xSectionPerElement;
  std::vector<G4double> runningSum;

  for (G4int i=0; i < nElements; i++) {
    xSectionPerElement = 0;
    for (G4int j=0; j < G4int(awicsPerElement[i].size()); j++)
                     xSectionPerElement += awicsPerElement[i][j];
    runningSum.push_back(theAtomsPerVolumeVector[i]*xSectionPerElement);
    crossSectionTotal += runningSum[i];
  }

  // Compare random number to running sum over element xc to choose Z

  // Initialize Z and A to first element and first isotope in case 
  // cross section is zero
   
  G4double ZZ = (*theElementVector)[0]->GetZ();
  G4double AA = AvaluesPerElement[0][0];
  if (crossSectionTotal != 0.) {
    G4double random = G4UniformRand();
    for(G4int i=0; i < nElements; i++) {
      if(i!=0) runningSum[i] += runningSum[i-1];
      if(random <= runningSum[i]/crossSectionTotal) {
        ZZ = ((*theElementVector)[i])->GetZ();

        // Compare random number to running sum over isotope xc to choose A

        G4int nIso = awicsPerElement[i].size();
        G4double* running = new G4double[nIso];
        for (G4int j=0; j < nIso; j++) {
          running[j] = awicsPerElement[i][j];
          if(j!=0) running[j] += running[j-1];
        }

        G4double trial = G4UniformRand(); 
        for (G4int j=0; j < nIso; j++) {
          AA = AvaluesPerElement[i][j];
          if (trial <= running[j]/running[nIso-1]) break;
        }
        delete [] running;
        break;
      }
    }
  }
  return std::pair<G4double, G4double>(ZZ, AA);
}
*/

void
G4CrossSectionDataStore::AddDataSet(G4VCrossSectionDataSet* aDataSet)
{
  DataSetList.push_back(aDataSet);
  NDataSetList++;
}

void
G4CrossSectionDataStore::
BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if(NDataSetList > 0) {
    for (G4int i=0; i<NDataSetList; i++) {
      DataSetList[i]->BuildPhysicsTable(aParticleType);
    } 
  }
}


void
G4CrossSectionDataStore::
DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if (NDataSetList == 0) {
    G4cout << "WARNING - G4CrossSectionDataStore::DumpPhysicsTable: "
	   << " no data sets registered"<<G4endl;
    return;
  }
  for (G4int i = NDataSetList-1; i >= 0; i--) {
    DataSetList[i]->DumpPhysicsTable(aParticleType);
  }
}
