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
// $Id: G4CrossSectionDataStore.cc,v 1.20 2010/11/19 08:14:44 gunter Exp $
// GEANT4 tag $Name: geant4-09-04 $
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
// 29.04.2010 G.Folger     modifictaions for integer A & Z
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
    if (DataSetList[i]->IsIsoApplicable(aParticle, anIsotope->GetZ(), anIsotope->GetN()))
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
  G4int iZ=G4lrint(Z);
  G4int iA=G4lrint(A);
  return GetCrossSection(aParticle, iZ, iA,aTemperature); 
  }


G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* aParticle,
                                         G4int Z, G4int A,
					 G4double aTemperature)
{
  if (NDataSetList == 0) 
  {
    throw G4HadronicException(__FILE__, __LINE__, 
     "G4CrossSectionDataStore: no data sets registered");
    return DBL_MIN;
  }
  for (G4int i = NDataSetList-1; i >= 0; i--) {
      if (DataSetList[i]->IsIsoApplicable(aParticle, Z, A))
      return DataSetList[i]->GetZandACrossSection(aParticle,Z,A,aTemperature);
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
  G4int ZZ = G4lrint(anElement->GetZ());
  G4int AA;

  // Collect abundance weighted cross sections and A values for each isotope 
  // in each element
 
  const G4int nIsoPerElement = anElement->GetNumberOfIsotopes();

  // user-defined isotope abundances
  if (0 < nIsoPerElement) { 
    G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
    AA =(*isoVector)[0]->GetN();
    if(1 < nIsoPerElement) {

      G4double* xsec  = new G4double [nIsoPerElement];
      G4double iso_xs = 0.0;
      G4double cross  = 0.0;

      G4double* abundVector = anElement->GetRelativeAbundanceVector();
      G4bool elementXS = false;
      for (i = 0; i<nIsoPerElement; i++) {
	if (inCharge->IsIsoApplicable(particle, ZZ,(*isoVector)[i]->GetN())) {
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
	  AA = (*isoVector)[i]->GetN();
	  break;
	}
      }
      delete [] xsec;
    }
    // natural abundances
  } else { 

    G4StableIsotopes theDefaultIsotopes;  
//-- Int ZZ    G4int Z = G4int(ZZ + 0.5);
    const G4int nIso = theDefaultIsotopes.GetNumberOfIsotopes(ZZ);
    G4int index = theDefaultIsotopes.GetFirstIsotope(ZZ);
    AA = theDefaultIsotopes.GetIsotopeNucleonCount(index);
    
    if(1 < nIso) {

      G4double* xsec  = new G4double [nIso];
      G4double iso_xs = 0.0;
      G4double cross  = 0.0;
      G4bool elementXS= false;

      for (i = 0; i<nIso; i++) {
        AA = theDefaultIsotopes.GetIsotopeNucleonCount(index+i);
	if (inCharge->IsIsoApplicable(particle, ZZ, AA )) {
	  iso_xs = inCharge->GetZandACrossSection(particle, ZZ, AA, aTemp);
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
	  AA = theDefaultIsotopes.GetIsotopeNucleonCount(index+i);
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
