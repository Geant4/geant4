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
//
// File name:     G4CrossSectionDataStore
//
// Modifications:
// 23.01.2009 V.Ivanchenko add destruction of data sets
// 29.04.2010 G.Folger     modifictaions for integer A & Z
// 14.03.2011 V.Ivanchenko fixed DumpPhysicsTable
// 15.08.2011 G.Folger, V.Ivanchenko, T.Koi, D.Wright redesign the class
// 07.03.2013 M.Maire cosmetic in DumpPhysicsTable
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "G4CrossSectionDataStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4Nucleus.hh"

#include "G4DynamicParticle.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include <algorithm>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4CrossSectionDataStore::G4CrossSectionDataStore()
  : nist(G4NistManager::Instance())
  , currentMaterial(nullptr)
  , matParticle(nullptr)
  , matKinEnergy(0.0)
  , matCrossSection(0.0)
  , nDataSetList(0)
  , verboseLevel(0)
  , fastPathFlags()
  , fastPathParams()
  , counters()
  , fastPathCache()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4CrossSectionDataStore::~G4CrossSectionDataStore()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         const G4Material* mat , G4bool requiresSlowPath)
{
	//The fast-path algorithm
	//   requiresSlowPath == true => Use slow path independently of other conditions
	//A. Dotti: modifications to this algorithm following the studies of P. Ruth and R. Fowler
	//          on speeding up the cross-sections calculations. Their algorithm is called in the
	// 	        following "fast-path" while the normal approach to calcualte cross-section is
	//          referred to as "slow-path".
	//Some details on the fast-path algorithm:
	//The idea is to use a cached approximation of the material cross-section.
	//Starting points:
	//1- We need the material cross-section for navigation purposes: e.g. to calculate the PIL.
	//2- If this interaction occurs at the end of a step we need to select the G4Element on which
	//   nucleus the interaction is actually happening. This is done calculating the element cross-section
	//   throwing a random number and selecting the appropriate nucleus (see SampleZandA function)
	//3- To calculate the material cross section needed for 1- we use the G4Element cross sections.
	//   Since the material cross-section is simply the weighted sum of the element cross-sections.
	//4- The slow path algorithm here accomplishes two use cases (not very good design IMHO): it
	//   calculates the material cross-section and it updates the xsecelem array that contains the
	//   relative cross-sections used by SampleZandA to select the element on which the interaction
	//   occurs.
	//The idea of the fast-path algorithm is to replace the request for 1- with a faster calculation
	//of the material cross-section from an approximation contained in cache.
	//There two complications to take into account:
	//A- If I use a fast parametrization for material cross-section, I still need to do the full
	//   calculations if an interaction occurs, because I need to calculate the element cross-sections.
	//   Since the function that updates xsecelem is the same (this one) I need to be sure that
	//   if I call this method for SampleAandZ the xsecelem is updated.
	//B- It exists the possibility to be even fast the the fast-path algorithm: this happens when
	//   to select the element of the interaction via SampleZandI I call again this method exactly
	//   with the same conditions as when I called this method to calculate the material cross-section.
	//   In such a case xsecelem is updated and we do not need to do much more. This happens when
	//   for example a neutron undergoes an interaction at the end of the step.
	//Dealing with B- complicates a bit the algorithm.
	//In summary:
	// If no fast-path algo is available (or user does not want that), go with the old plain algorithm
	// If a fast-path algo is avilable, use it whenever possible.
	//
	// In general we expect user to selectively decide for which processes, materials and particle combinations
	// we want to use the fast-path. If this is activated we also expect that the cross-section fast-path
	// cache is created during the run initialization via calls to this method.
	//
	//fastPathFlags contains control flags for the fast-path algorithm:
	// 	       .prevCalcUsedFastPath == true => Previous call to GetCrossSection used the fast-path
	//             it is used in the decision to assess if xsecelem is correctly set-up
	//	       .useFastPathIfAvailable == true => User requested the use of fast-path algorithm
	//	       .initializationPhase == true => If true we are in Geant4 Init phase before the event-loop

	//Check user-request, does he want fast-path? if not
	// OR
	// we want fast-path and we are in initialization phase?
	if ( !fastPathFlags.useFastPathIfAvailable
		||	(fastPathFlags.useFastPathIfAvailable&&fastPathFlags.initializationPhase) ) {
		//Traditional algorithm is requested
		requiresSlowPath=true;
	}

	//Logging for performance calculations and counter, active only in FPDEBUG mode
	counters.MethodCalled();
	//Measure number of cycles
	G4FastPathHadronicCrossSection::logStartCountCycles(timing);

	//This is the cache entry of the fast-path cross-section parametrization
	G4FastPathHadronicCrossSection::cycleCountEntry* entry = nullptr;
	//Did user request fast-path in first place and are we not in the initialization phase
	if ( fastPathFlags.useFastPathIfAvailable && !fastPathFlags.initializationPhase ) {
		//Important: if it is in initialization phase we should NOT use fast path: we are going to build it
		//G4FastPathHadronicCrossSection::G4CrossSectionDataStore_Key searchkey = {part->GetParticleDefinition(),mat};
		entry = fastPathCache[{part->GetParticleDefinition(),mat}];
	}

  //Super-fast-path: are we calling again this method for exactly the same conditions
  //of the triplet {particle,material,energy}?
  if(mat == currentMaterial && part->GetDefinition() == matParticle
     && part->GetKineticEnergy() == matKinEnergy) 
    {
	  G4FastPathHadronicCrossSection::logInvocationTriedOneLine(entry);
	  //If there is no user-request for the fast-path in first place?
	  //It means we built the xsecelem for sure, let's return immediately
	  if ( !fastPathFlags.useFastPathIfAvailable ) {
		  return matCrossSection;
	  } else {
		  //Check that the last time we called this method we used the slow
		  //path: we need the data-member xsecelem to be setup correctly for the current
		  //interaction. This is ensured only if: we will do the slow path right now or we
		  //did it exactly for the same conditions of the last call.
		  if ( !fastPathFlags.prevCalcUsedFastPath && ! requiresSlowPath ) {
			  counters.HitOneLine();
			  G4FastPathHadronicCrossSection::logInvocationOneLine(entry);
			  //Good everything is setup correctly, exit!
			  return matCrossSection;
		  } else {
			  //We need to follow the slow-path because
			  //xsecelem is not calculated correctly
			  requiresSlowPath = true;
		  }
	  }
    }
  
  //Ok, now check if we have cached for this {particle,material,energy} the cross-section
  //in this case let's return immediately, if we are not forced to take the slow path
  //(e.g. as before if the xsecelem is not up-to-date we need to take the slow-path).
  //Note that this is not equivalent to the previous ultra-fast check: we now have a map here
  //So we can have for example a different particle.
  if ( entry != nullptr && entry->energy == part->GetKineticEnergy() ) {
	  G4FastPathHadronicCrossSection::logHit(entry);
	  if ( !requiresSlowPath ) {
		  return entry->crossSection;
	  }
  }

  currentMaterial = mat;
  matParticle = part->GetDefinition();
  matKinEnergy = part->GetKineticEnergy();
  matCrossSection = 0;

  //Now check if the cache entry has a fast-path cross-section calculation available
  G4FastPathHadronicCrossSection::fastPathEntry* fast_entry = nullptr;
  if ( entry != nullptr && ! requiresSlowPath ) {
	  fast_entry = entry->fastPath;
	  assert(fast_entry!=nullptr && !fastPathFlags.initializationPhase);
  }

  //Each fast-path cross-section has a minimum value of validity, if energy is below
  //that skip fast-path algorithm
  if ( fast_entry != nullptr && part->GetKineticEnergy() < fast_entry->min_cutoff )
  {
	  assert(requiresSlowPath==false);
	  requiresSlowPath = true;
  }

  //Ready to use the fast-path calculation
  if ( !requiresSlowPath && fast_entry != nullptr ) {
	  counters.FastPath();
	  //Retrieve cross-section from fast-path cache
	  matCrossSection = fast_entry->GetCrossSection(part->GetKineticEnergy());
	  fastPathFlags.prevCalcUsedFastPath=true;
  } else {
	  counters.SlowPath();
	  //Remember that we are now doing the full calculation: xsecelem will
	  //be made valid
	  fastPathFlags.prevCalcUsedFastPath=false;

	  G4int nElements = mat->GetNumberOfElements();
	  const G4double* nAtomsPerVolume = mat->GetVecNbOfAtomsPerVolume();

	  if(G4int(xsecelm.size()) < nElements) { xsecelm.resize(nElements); }

	  for(G4int i=0; i<nElements; ++i) {
		  matCrossSection += nAtomsPerVolume[i] *
				  GetCrossSection(part, (*mat->GetElementVector())[i], mat);
		  xsecelm[i] = matCrossSection;
	  }
  }
  //Stop measurement of cpu cycles
  G4FastPathHadronicCrossSection::logStopCountCycles(timing);

  if ( entry != nullptr ) {
	  entry->energy = part->GetKineticEnergy();
	  entry->crossSection = matCrossSection;
  }
  //Some logging of timing
  G4FastPathHadronicCrossSection::logTiming(entry,fast_entry,timing);
  return matCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
G4CrossSectionDataStore::DumpFastPath(const G4ParticleDefinition* pd, const G4Material* mat,std::ostream& os)
{
	const G4FastPathHadronicCrossSection::cycleCountEntry* entry = fastPathCache[{pd,mat}];
	if ( entry != nullptr ) {
		if ( entry->fastPath != nullptr ) {
			os<<*entry->fastPath;
		} else {
			os<<"#Cache entry for {"<<(pd!=nullptr?pd->GetParticleName():"UNDEFINED")<<",";
			os<<(mat!=nullptr?mat->GetName():"UNDEFINED")<<"} found, but no fast path defined";
		}
	} else {
		os<<"#Cache entry for {"<<(pd!=nullptr?pd->GetParticleName():"UNDEFINED")<<",";
		os<<(mat!=nullptr?mat->GetName():"UNDEFINED")<<"} not found.";
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double 
G4CrossSectionDataStore::ComputeCrossSection(const G4DynamicParticle* part,
					     const G4Material* mat)
{
  if(part->GetKineticEnergy() == matKinEnergy && mat == currentMaterial &&
     part->GetDefinition() == matParticle)
    return matCrossSection;

  currentMaterial = mat;
  matParticle = part->GetDefinition();
  matKinEnergy = part->GetKineticEnergy();
  matCrossSection = 0.0;

  size_t nElements = mat->GetNumberOfElements();
  const G4double* nAtomsPerVolume = mat->GetVecNbOfAtomsPerVolume();

  if(xsecelm.size() < nElements) { xsecelm.resize(nElements); }

  for(size_t i=0; i<nElements; ++i) {
    matCrossSection += nAtomsPerVolume[i] *
      GetCrossSection(part, mat->GetElement(i), mat);
    xsecelm[i] = matCrossSection;
  }
  return matCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                                  const G4Element* elm,
                                                  const G4Material* mat)
{
  G4int i = nDataSetList - 1;
  G4int Z = elm->GetZasInt();

  if(elm->GetNaturalAbundanceFlag() &&
     dataSetList[i]->IsElementApplicable(part, Z, mat))
  {
    // element wise cross section
    return dataSetList[i]->GetElementCrossSection(part, Z, mat);
  }

  // isotope wise cross section
  size_t nIso = elm->GetNumberOfIsotopes();

  // user-defined isotope abundances
  const G4double* abundVector = elm->GetRelativeAbundanceVector();

  G4double sigma = 0.0;

  for(size_t j = 0; j < nIso; ++j)
  {
    const G4Isotope* iso = elm->GetIsotope(j);
    sigma += abundVector[j] *
             GetIsoCrossSection(part, Z, iso->GetN(), iso, elm, mat, i);
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetIsoCrossSection(const G4DynamicParticle* part,
					    G4int Z, G4int A, 
					    const G4Isotope* iso,
					    const G4Element* elm,
					    const G4Material* mat, 
					    G4int idx)
{
  // this methods is called after the check that dataSetList[idx] 
  // depend on isotopes, so for this DataSet only isotopes are checked

  // isotope-wise cross section does exist
  if(dataSetList[idx]->IsIsoApplicable(part, Z, A, elm, mat) ) {
    return dataSetList[idx]->GetIsoCrossSection(part, Z, A, iso, elm, mat);

  } else {
    // seach for other dataSet
    for (G4int j = nDataSetList-1; j >= 0; --j) { 
      if (dataSetList[j]->IsElementApplicable(part, Z, mat)) {
	return dataSetList[j]->GetElementCrossSection(part, Z, mat);
      } else if (dataSetList[j]->IsIsoApplicable(part, Z, A, elm, mat)) {
	return dataSetList[j]->GetIsoCrossSection(part, Z, A, iso, elm, mat);
      }
    }
  }
  G4ExceptionDescription ed;
  ed << "No isotope cross section found for " 
     << part->GetDefinition()->GetParticleName() 
     << " off Element " << elm->GetName()
     << "  in " << mat->GetName() << " Z= " << Z << " A= " << A
     << " E(MeV)= " << part->GetKineticEnergy()/MeV << G4endl; 
  G4Exception("G4CrossSectionDataStore::GetIsoCrossSection", "had001", 
              FatalException, ed);
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         G4int Z, G4int A,
					 const G4Isotope* iso,
                                         const G4Element* elm,
					 const G4Material* mat)
{
  for (G4int i = nDataSetList-1; i >= 0; --i) {
    if (dataSetList[i]->IsIsoApplicable(part, Z, A, elm, mat) ) {
      return dataSetList[i]->GetIsoCrossSection(part, Z, A, iso, elm, mat);
    }
  }
  G4ExceptionDescription ed;
  ed << "No isotope cross section found for " 
     << part->GetDefinition()->GetParticleName() 
     << " off Element " << elm->GetName()
     << "  in " << mat->GetName() << " Z= " << Z << " A= " << A
     << " E(MeV)= " << part->GetKineticEnergy()/MeV << G4endl; 
  G4Exception("G4CrossSectionDataStore::GetCrossSection", "had001", 
              FatalException, ed);
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const G4Element*
G4CrossSectionDataStore::SampleZandA(const G4DynamicParticle* part, 
                                     const G4Material* mat,
				     G4Nucleus& target)
{
  size_t nElements = mat->GetNumberOfElements();
  const G4Element* anElement = mat->GetElement(0);

  // select element from a compound 
  if(1 < nElements) {
    G4double cross = matCrossSection*G4UniformRand();
    for(size_t i=0; i<nElements; ++i) {
      if(cross <= xsecelm[i]) {
	anElement = mat->GetElement(i);
        break;
      }
    }
  }

  G4int Z = anElement->GetZasInt();
  const G4Isotope* iso = nullptr;

  G4int i = nDataSetList-1;
  if (dataSetList[i]->IsElementApplicable(part, Z, mat)) {

    //----------------------------------------------------------------
    // element-wise cross section
    // isotope cross section is not computed
    //----------------------------------------------------------------
    size_t nIso = anElement->GetNumberOfIsotopes();
    iso = anElement->GetIsotope(0);

    // more than 1 isotope
    if(1 < nIso) { 
      iso = dataSetList[i]->SelectIsotope(anElement, 
                                          part->GetKineticEnergy(),
					  part->GetLogKineticEnergy());
    }
  } else {

    //----------------------------------------------------------------
    // isotope-wise cross section
    // isotope cross section is computed
    //----------------------------------------------------------------
    size_t nIso = anElement->GetNumberOfIsotopes();
    iso = anElement->GetIsotope(0);

    // more than 1 isotope
    if(1 < nIso) {
      const G4double* abundVector = anElement->GetRelativeAbundanceVector();
      if(xseciso.size() < nIso) { xseciso.resize(nIso); }

      G4double cross = 0.0;
      size_t j;
      for (j = 0; j<nIso; ++j) {
	G4double xsec = 0.0;
	if(abundVector[j] > 0.0) {
	  iso = anElement->GetIsotope(j);
	  xsec = abundVector[j]*
	    GetIsoCrossSection(part, Z, iso->GetN(), iso, anElement, mat, i);
	}
	cross += xsec;
	xseciso[j] = cross;
      }
      cross *= G4UniformRand();
      for (j = 0; j<nIso; ++j) {
	if(cross <= xseciso[j]) {
	  iso = anElement->GetIsotope(j);
	  break;
	}
      }
    }
  }
  target.SetIsotope(iso);
  return anElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
G4CrossSectionDataStore::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if (nDataSetList == 0) {
    G4ExceptionDescription ed;
    ed << "No cross section is registered for " 
       << aParticleType.GetParticleName() << G4endl;
    G4Exception("G4CrossSectionDataStore::BuildPhysicsTable", "had001", 
                FatalException, ed);
    return;
  }
  for (G4int i=0; i<nDataSetList; ++i) {
    dataSetList[i]->BuildPhysicsTable(aParticleType);
  } 
  //A.Dotti: if fast-path has been requested we can now create the surrogate
  //         model for fast path.
  if ( fastPathFlags.useFastPathIfAvailable ) {
	  fastPathFlags.initializationPhase = true;
	  using my_value_type=G4FastPathHadronicCrossSection::G4CrossSectionDataStore_Requests::value_type;
	  //Loop on all requests, if particle matches create the corresponding fsat-path
	  std::for_each( requests.begin() , requests.end() ,
			  [&aParticleType,this](const my_value_type& req) {
  	  	  	  	  if ( aParticleType == *req.part_mat.first ) {
  	  	  	  		  G4FastPathHadronicCrossSection::cycleCountEntry* entry =
  	  	  	  				  new G4FastPathHadronicCrossSection::cycleCountEntry(aParticleType.GetParticleName(),req.part_mat.second);
  	  	  	  		  entry->fastPath =
  	  	  	  				  new G4FastPathHadronicCrossSection::fastPathEntry(&aParticleType,req.part_mat.second,req.min_cutoff);
  	  	  	  		  entry->fastPath->Initialize(this);
  	  	  	  		  fastPathCache[req.part_mat] = entry;
  	  	  	  	  }
	  	  	  }
	  	  );
	  fastPathFlags.initializationPhase = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::ActivateFastPath( const G4ParticleDefinition* pdef, 
     const G4Material* mat, G4double min_cutoff)
{
  assert(pdef!=nullptr&&mat!=nullptr);
  G4FastPathHadronicCrossSection::G4CrossSectionDataStore_Key key={pdef,mat};
  if ( requests.insert( { key , min_cutoff } ).second ) {
    G4ExceptionDescription ed;
    ed << "Attempting to request FastPath for couple: <"
       << pdef->GetParticleName() << ", " <<mat->GetName() 
       << "> but combination already exists" << G4endl;
    G4Exception("G4CrossSectionDataStore::ActivateFastPath", "had001", 
                FatalException, ed);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4CrossSectionDataStore::DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  // Print out all cross section data sets used and the energies at
  // which they apply

  if (nDataSetList == 0) {
    G4cout << "WARNING - G4CrossSectionDataStore::DumpPhysicsTable: "
	   << " no data sets registered" << G4endl;
    return;
  }

  for (G4int i = nDataSetList-1; i >= 0; --i) {
    G4double e1 = dataSetList[i]->GetMinKinEnergy();
    G4double e2 = dataSetList[i]->GetMaxKinEnergy();
     G4cout 
      << "     Cr_sctns: " << std::setw(25) << dataSetList[i]->GetName() << ": "
      <<  G4BestUnit(e1, "Energy")
      << " ---> "
      <<  G4BestUnit(e2, "Energy") << "\n";
      if (dataSetList[i]->GetName() == "G4CrossSectionPairGG") {
        dataSetList[i]->DumpPhysicsTable(aParticleType);
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
#include <typeinfo>
void G4CrossSectionDataStore::DumpHtml(const G4ParticleDefinition& /* pD */,
                                       std::ofstream& outFile) const
{
  // Write cross section data set info to html physics list
  // documentation page

  G4double ehi = 0;
  G4double elo = 0;
  G4String physListName(std::getenv("G4PhysListName"));
  for (G4int i = nDataSetList-1; i > 0; i--) {
    elo = dataSetList[i]->GetMinKinEnergy()/GeV;
    ehi = dataSetList[i]->GetMaxKinEnergy()/GeV;
    outFile << "      <li><b><a href=\"" << physListName << "_"
	         << dataSetList[i]->GetName() << ".html\"> "
            << dataSetList[i]->GetName() << "</a> from "
            << elo << " GeV to " << ehi << " GeV </b></li>\n";
	 //G4cerr << i << ": XS for " << pD.GetParticleName() << " : " << dataSetList[i]->GetName() 
	 //       << " typeid : " << typeid(dataSetList[i]).name()<< G4endl;			
	 PrintCrossSectionHtml(dataSetList[i]);			
  }

  G4double defaultHi = dataSetList[0]->GetMaxKinEnergy()/GeV;
  if (ehi < defaultHi) {
    outFile << "      <li><b><a href=\"" << dataSetList[0]->GetName() << ".html\"> "
            << dataSetList[0]->GetName() << "</a> from "
            << ehi << " GeV to " << defaultHi << " GeV </b></li>\n";
	 PrintCrossSectionHtml(dataSetList[0]);			
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::PrintCrossSectionHtml(const G4VCrossSectionDataSet *cs) const
{
  G4String dirName(std::getenv("G4PhysListDocDir"));
  G4String physListName(std::getenv("G4PhysListName"));

	G4String pathName = dirName + "/" + physListName + "_" + HtmlFileName(cs->GetName());
	std::ofstream outCS;
	outCS.open(pathName);
	outCS << "<html>\n";
	outCS << "<head>\n";
	outCS << "<title>Description of " << cs->GetName() 
		 << "</title>\n";
	outCS << "</head>\n";
	outCS << "<body>\n";

	cs->CrossSectionDescription(outCS);

	outCS << "</body>\n";
	outCS << "</html>\n";

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4String G4CrossSectionDataStore::HtmlFileName(const G4String & in) const
{
   G4String str(in);
   // replace blanks by _  C++11 version:
   std::transform(str.begin(), str.end(), str.begin(), [](char ch) {
       return ch == ' ' ? '_' : ch;
   });
   str=str + ".html";		
   return str;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::AddDataSet(G4VCrossSectionDataSet* p)
{
  if(p->ForAllAtomsAndEnergies()) { 
    dataSetList.clear();
    nDataSetList = 0;
  }
  dataSetList.push_back(p);
  ++nDataSetList;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::AddDataSet(G4VCrossSectionDataSet* p, size_t i )  
{
  if(p->ForAllAtomsAndEnergies()) {
    dataSetList.clear();
    dataSetList.push_back(p);
    nDataSetList = 1;
  } else {
    if ( i > dataSetList.size() ) i = dataSetList.size(); 
    std::vector< G4VCrossSectionDataSet* >::iterator it = dataSetList.end() - i;
    dataSetList.insert(it , p);
    ++nDataSetList;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
