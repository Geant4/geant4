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
//      Geant4 source file 
//
//      File name: G4ParticleHPProbabilityTablesStore.cc
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//               Loic Thulliez (CEA France)      
//
//      Creation date: 4 June 2024
//
//      Description: Class to store all probability tables for
//                   different isotopes and in future also for
//                   different temperatures.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//

#include "G4ParticleHPProbabilityTablesStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleHPVector.hh"
#include "G4ParticleHPManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4DynamicParticle.hh"
#include "G4Exception.hh"
#include "G4ParticleHPIsoProbabilityTable.hh"
#include "G4ParticleHPIsoProbabilityTable_NJOY.hh"
#include "G4ParticleHPIsoProbabilityTable_CALENDF.hh"
#include <string>
#include <sstream>

G4ParticleHPProbabilityTablesStore* G4ParticleHPProbabilityTablesStore::instance = nullptr;

///--------------------------------------------------------------------------------------
G4ParticleHPProbabilityTablesStore::G4ParticleHPProbabilityTablesStore() : 
  Temperatures(0), ProbabilityTables(0), URRlimits(0), usedNjoy(0), usedCalendf(0) {
  if ( G4FindDataDir( "G4URRPTDATA" ) != nullptr ) {
    dirName = G4FindDataDir( "G4URRPTDATA" );
  } else {
    G4Exception( "G4ParticleHPProbabilityTablesStore::G4ParticleHPProbabilityTablesStore()", "hadhp01",
  		 FatalException, "Please setenv G4URRPTDATA, it is not defined." );
  }
  if ( G4ParticleHPManager::GetInstance()->GetUsedPTformat() == "njoy" ) {
    dirName += "/ProbabilityTables/njoy/";
    usedNjoy = true;
  } else if ( G4ParticleHPManager::GetInstance()->GetUsedPTformat() == "calendf" ) {
    dirName += "/ProbabilityTables/calendf/";
    usedCalendf = true;
  } else {
    G4Exception( "G4ParticleHPProbabilityTablesStore::G4ParticleHPProbabilityTablesStore()", "hadhp01",
  		 FatalException, "The format of probability tables is not set properly, "
  		 "please check variable USE_PROBABILITY_TABLE_FROM in G4ParticleHPManager." );
  }
  numIso = (G4int)G4Isotope::GetNumberOfIsotopes();
  // find all possible temperatures for all isotopes.
  G4Material* theMaterial;
  G4int Temp;
  Temperatures = new std::vector< std::vector< G4int > >;
  for ( G4int i = 0; i < numIso; ++i ) {
    std::vector< G4int > isotemperatures;
    std::map< std::thread::id, G4double > simpleMapE;
    std::map< std::thread::id, G4double > simpleMapRN;
    Temperatures->push_back( isotemperatures );
    energy_cache.push_back( simpleMapE );
    random_number_cache.push_back( simpleMapRN );
  }
  for ( std::size_t im = 0; im < G4Material::GetNumberOfMaterials(); ++im ) {
    std::vector< G4bool > isoinmat( numIso, false );
    theMaterial = G4Material::GetMaterialTable()->at(im);
    Temp = (G4int)theMaterial->GetTemperature();
    for ( G4int e = 0; e < (G4int)theMaterial->GetNumberOfElements(); ++e ) {
      for ( G4int i = 0; i < (G4int)theMaterial->GetElement(e)->GetNumberOfIsotopes(); ++i ) {
    	std::size_t indexI = theMaterial->GetElement(e)->GetIsotope(i)->GetIndex();
    	std::vector< G4int > isotemperatures = (*Temperatures)[indexI];
    	if ( std::find( isotemperatures.begin(), isotemperatures.end(), Temp ) == isotemperatures.end() ) {
    	  (*Temperatures)[indexI].push_back( Temp );
    	  isoinmat[indexI] = true;
    	} else {
    	  if ( isoinmat[indexI] ) {
    	    G4ExceptionDescription message;
    	    message << "The isotope Z=" << theMaterial->GetElement(e)->GetIsotope(i)->GetZ() 
                    << " and A=" << theMaterial->GetElement(e)->GetIsotope(i)->GetN() 
                    << " is more times in material " << theMaterial->GetName() << ".\n"
    		    << "This may cause bias in selection of target isotopes, if there are elements with different URR limits in the material.\n"
    		    << "Please make materials only with elements with different Z.";
    	    G4Exception( "G4ParticleHPProbabilityTablesStore::InitURRlimits()", "hadhp01", JustWarning, message );
    	  }
    	}
      }  // end of cycle over isotopes			
    }  // end of cycle over elements
  }  // end of cycle over materials
}

///--------------------------------------------------------------------------------------
G4ParticleHPProbabilityTablesStore::~G4ParticleHPProbabilityTablesStore() {
  for ( G4int i = 0; i < numIso; i++ ) {
    for ( std::map< G4int, G4ParticleHPIsoProbabilityTable* >::iterator it = (*ProbabilityTables)[i].begin(); 
          it != (*ProbabilityTables)[i].end(); ++it ) {
      delete it->second;
    }
  }
  delete ProbabilityTables;
  delete URRlimits;
  delete Temperatures;
}

///--------------------------------------------------------------------------------------
G4ParticleHPProbabilityTablesStore* G4ParticleHPProbabilityTablesStore::GetInstance() {
  if ( instance == nullptr ) instance = new G4ParticleHPProbabilityTablesStore;
  return instance;
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPProbabilityTablesStore::GetIsoCrossSectionPT( const G4DynamicParticle* dp, G4int MTnumber, 
  const G4Isotope* iso, const G4Element* ele, const G4Material* mat ) {
  G4double kineticEnergy = dp->GetKineticEnergy();
  std::size_t indexI = iso->GetIndex();
  G4int T = (G4int)( mat->GetTemperature() );
  std::thread::id id = std::this_thread::get_id();
  if ( energy_cache[indexI][id] == kineticEnergy ) {
    return (*ProbabilityTables)[indexI][T]->GetCorrelatedIsoCrossSectionPT( dp, MTnumber, ele, kineticEnergy, random_number_cache[indexI][id], id );
  } else {
    energy_cache[indexI][id] = kineticEnergy;
    return (*ProbabilityTables)[iso->GetIndex()][T]->GetIsoCrossSectionPT( dp, MTnumber, ele, kineticEnergy, random_number_cache[indexI], id );
  }
}

///--------------------------------------------------------------------------------------
void G4ParticleHPProbabilityTablesStore::Init() {
  ProbabilityTables = new std::vector< std::map< G4int, G4ParticleHPIsoProbabilityTable* > >;
  for ( G4int i = 0; i < numIso; i++ ) {
    G4int Z = (*(G4Isotope::GetIsotopeTable()))[i]->GetZ();
    G4int A = (*(G4Isotope::GetIsotopeTable()))[i]->GetN();
    G4int meso = (*(G4Isotope::GetIsotopeTable()))[i]->Getm();
    std::map< G4int, G4ParticleHPIsoProbabilityTable* > tempPTmap;
    for ( G4int T : (*Temperatures)[i] ) {  // create probability table for each temperature and isotope
      if ( usedNjoy ) {
  	tempPTmap.insert( { T, new G4ParticleHPIsoProbabilityTable_NJOY } );
      } else if ( usedCalendf ) {
  	tempPTmap.insert( { T, new G4ParticleHPIsoProbabilityTable_CALENDF } );
      }
      tempPTmap[T]->Init( Z, A, meso, T, dirName );
    }
    ProbabilityTables->push_back( tempPTmap );
  }
}

///--------------------------------------------------------------------------------------
void G4ParticleHPProbabilityTablesStore::InitURRlimits() {
  URRlimits = new std::vector< std::pair< G4double, G4double > >;
  // Find energy limits of the URR for each element.
  std::vector< G4bool > wasnotprintedyet( numIso, true );
  for ( std::size_t i = 0; i < G4Element::GetNumberOfElements(); ++i ) {
    G4double minE = DBL_MAX;
    G4double maxE = 0.0;
    for ( G4int ii = 0; ii < (G4int)(*(G4Element::GetElementTable()))[i]->GetNumberOfIsotopes(); ++ii ) {
      G4int Z = (*(G4Element::GetElementTable()))[i]->GetIsotope(ii)->GetZ();
      G4int A = (*(G4Element::GetElementTable()))[i]->GetIsotope(ii)->GetN();
      G4int meso = (*(G4Element::GetElementTable()))[i]->GetIsotope(ii)->Getm();
      std::size_t indexI = (*(G4Element::GetElementTable()))[i]->GetIsotope(ii)->GetIndex();
      filename = std::to_string(Z) + "_" + std::to_string(A);
      if ( meso != 0 ) filename += "_m" + std::to_string(meso);
      G4double emin = DBL_MAX;
      G4double emax = 0.0;
      G4double emin2 = DBL_MAX;
      G4double emax2 = 0.0;
      G4bool hasalreadyPT = false;
      G4bool doesnothavePTyet = false;
      for ( G4int T : (*Temperatures)[indexI] ) {
        G4String fullPathFileName = dirName + filename + "." + std::to_string(T) + ".pt.z"; 
        std::istringstream theData( std::ios::in );
        G4ParticleHPManager::GetInstance()->GetDataStream( fullPathFileName, theData );
        if ( theData.good() &&  hasalreadyPT ) {
  	  theData >> emin2 >> emax2;
  	  if ( emin2 != emin  ||  emax2 != emax ) {
  	    G4ExceptionDescription message;
  	    message << "There is mismatch between limits of unresolved resonance region for isotope " << "\n"
  		    << "with Z=" << Z << " and A=" << A << " for different temperatures, which should not happen, CHECK YOUR DATA!" << "\n"
                    << "The broadest limits will be used.";
  	    G4Exception( "G4ParticleHPProbabilityTablesStore::InitURRlimits()", "hadhp01", JustWarning, message );
  	    G4cout << "Temperature T= " << T << " emin=" << emin << " emax=" << emax << " minE=" << minE << " maxE=" << maxE << G4endl;
  	    if ( emin < minE ) minE = emin;
  	    if ( emax > maxE ) maxE = emax;
  	  } 
        } else if ( theData.good() ) {
  	  theData >> emin >> emax;
  	  if ( emin < minE ) minE = emin;
  	  if ( emax > maxE ) maxE = emax;
  	  hasalreadyPT = true;
  	  if ( doesnothavePTyet  &&  wasnotprintedyet[indexI] ) {
  	    G4ExceptionDescription message;
  	    message << "There are probability tables only for some of the used temperatures for isotope with Z=" << Z << " and A=" << A << ".\n" 
                    << "Smooth cross section will be used for other temperature(s) of this isotope, which may cause severe differences in simulation.";
  	    G4Exception( "G4ParticleHPProbabilityTablesStore::InitURRlimits()", "hadhp01", JustWarning, message );
  	    wasnotprintedyet[indexI] = false;
  	  }
  	} else if ( hasalreadyPT  &&  wasnotprintedyet[indexI] ) {
  	  G4ExceptionDescription message;
  	  message << "There are probability tables only for some of the used temperatures for isotope with Z=" << Z << " and A=" << A << ".\n" 
                  << "Smooth cross section will be used for other temperature(s) of this isotope, which may cause severe differences in simulation.";
  	  G4Exception( "G4ParticleHPProbabilityTablesStore::InitURRlimits()", "hadhp01", JustWarning, message );
  	  wasnotprintedyet[indexI] = false;
  	} else {
  	  doesnothavePTyet = true;
  	}
      }	
    }
    std::pair< G4double, G4double > simplepair( minE * eV, maxE * eV );
    URRlimits->push_back( simplepair );
  }
  // Find absolute energy limits of the URR and set them at the end of URRlimits.
  G4double absminE = (*URRlimits)[0].first;
  G4double absmaxE = (*URRlimits)[0].second;
  for ( auto iterator : (*URRlimits) ) {
    if ( iterator.first < absminE )  absminE = iterator.first;
    if ( iterator.second > absmaxE ) absmaxE = iterator.second;
  }
  std::pair< G4double, G4double > minmaxpair( absminE, absmaxE );
  URRlimits->push_back( minmaxpair );
}
