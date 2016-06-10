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
#ifndef G4FastPathHadronicCrossSection_hh
#define G4FastPathHadronicCrossSection_hh

#include "G4PhysicsFreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include <functional>
#include <utility>
#include <unordered_map>
#include <iostream>
#include <set>
#include <stdint.h>

class G4DynamicParticle;
class G4Material;
class G4CrossSectionDataStore;

//To measure performances and debug info on fast cross-section enable this
//#define FPDEBUG

//TODO: Move all logging and debug functionality to separate header
namespace G4FastPathHadronicCrossSection {
	//This data type contains the simplified representation of the
	//cross-section, by default it is a G4PhysicsVector type
	using XSParam=G4PhysicsFreeVector;
	//The key used to search in the cache.
	using G4CrossSectionDataStore_Key=std::pair<const G4ParticleDefinition*,const G4Material*>;
	//This represents the fast XS implementation.
	struct fastPathEntry{
		//fastPathEntry();
		fastPathEntry(const G4ParticleDefinition *par,const G4Material* mat,G4double min_cutoff);
		~fastPathEntry();
		inline G4double GetCrossSection(G4double ene) const { return physicsVector->Value(ene); }
		void Initialize(G4CrossSectionDataStore* );
		const G4ParticleDefinition * const particle;
		const G4Material * const material;
		const G4double min_cutoff;

		XSParam *physicsVector;
#       ifdef FPDEBUG
		//stats for debug
		G4int count;
		G4double slowpath_sum; //sum of all slowpath xs
		G4double max_delta;
		G4double min_delta;
		G4double sum_delta;
		G4double sum_delta_square;
#		endif
	};

	//A cache entry.
	struct cycleCountEntry{
		cycleCountEntry(const G4String& pname , const G4Material* mat);
		~cycleCountEntry();
		const G4String& particle;
		const G4Material * const material;

		//optional fastPathEntry
		fastPathEntry* fastPath;

		//cache per element of material test
		G4double energy;
		G4double crossSection;
#	  ifdef FPDEBUG
		uint64_t cacheHitCount;//
		uint64_t initCyclesFastPath;
		uint64_t invocationCountSlowPath;
		uint64_t totalCyclesSlowPath;
		uint64_t invocationCountFastPath;
		uint64_t totalCyclesFastPath;
		uint64_t invocationCountTriedOneLineCache;//
		uint64_t invocationCountOneLineCache;//
#	  endif
	};

	struct timing {
		unsigned long long rdtsc_start;
		unsigned long long rdtsc_stop;
	};

	struct getCrossSectionCount {
		getCrossSectionCount();
		inline void MethodCalled();
		inline void HitOneLine();
		inline void FastPath();
		inline void SlowPath();
		inline void SampleZandA();
#ifdef FPDEBUG
		uint64_t methodCalled;
		uint64_t hitOneLineCache;
		uint64_t fastPath;
		uint64_t slowPath;
		uint64_t sampleZandA;
#endif
	};

	//Hashing the key
	struct G4CrossSectionDataStore_Key_Hash {
		std::hash<uint64_t> hash_uint64_t;
		inline size_t operator()(const G4CrossSectionDataStore_Key& x) const throw() {
			return hash_uint64_t(hash_uint64_t( ((uint64_t)(x.first)) ) +  hash_uint64_t(((uint64_t)(x.second))));
		}
	};
	//Equality for two key elements
	struct G4CrossSectionDataStore_Key_EqualTo {
		inline bool operator()(const G4CrossSectionDataStore_Key& lhs, const G4CrossSectionDataStore_Key& rhs ) const {
			//TODO: Verify this: particles are singletons, materials use operator==
		        //TODO: in ref-10, G4Material::operator== becomes deleted, investigating why
			return (lhs.first==rhs.first)&&(lhs.second == rhs.second);
		}
	};
//	The cache itself
	using G4CrossSectionDataStore_Cache=std::unordered_map<G4CrossSectionDataStore_Key,cycleCountEntry*,
			G4CrossSectionDataStore_Key_Hash,G4CrossSectionDataStore_Key_EqualTo>;

	struct fastPathRequestConfig_t {
		G4CrossSectionDataStore_Key part_mat;
		G4double min_cutoff;
	};
	//Two of the elements are identical if the part_mat part is
	struct fastPathRequestConfig_Less {
		std::less<G4CrossSectionDataStore_Key> less;
		inline bool operator()(const fastPathRequestConfig_t& lhs,const fastPathRequestConfig_t& rhs ) const {
			return less(lhs.part_mat,rhs.part_mat);
		}
	};
	using G4CrossSectionDataStore_Requests=std::set<fastPathRequestConfig_t,fastPathRequestConfig_Less>;

	//Configure the caching mechanism
	struct controlFlag {
		G4bool prevCalcUsedFastPath;
		G4bool useFastPathIfAvailable;
		G4bool initializationPhase;
		controlFlag() : prevCalcUsedFastPath(false),useFastPathIfAvailable(false),initializationPhase(false) {}
	};
	//Parameters to control sampling
	struct fastPathParameters {
		fastPathParameters() {
			//default
			//TODO: are these ok?
			  queryMax = 10000;
			  sampleMin = 0.0001;
			  sampleMax = 10000;
			  sampleCount = 200000;
			  dpTol = 0.01;
		}
	  //PRUTH vars for sampling and surragate model
	  G4double queryMax;
	  G4double sampleMin;
	  G4double sampleMax;
	  G4int sampleCount;
	  G4double dpTol;
	};

	//Logging functionalities, disabled if not in FPDEBUG mode
	static inline void logInvocationTriedOneLine( cycleCountEntry* );
	static inline void logInvocationOneLine( cycleCountEntry* );
	static inline void logHit(cycleCountEntry*);
	static inline void logInvocationCountFastPath( cycleCountEntry* );
	static inline void logInvocationCountSlowPAth( cycleCountEntry* );
#ifdef FPDEBUG
	void logStartCountCycles( timing& );
	void logStopCountCycles( timing& );
#else
	inline void logStartCountCycles(timing&) {}
	inline void logStopCountCycles(timing&) {}
#endif
	static inline void logInitCyclesFastPath( cycleCountEntry* , timing& );
	static inline void logTotalCyclesFastPath( cycleCountEntry* , timing& );
	static inline void logTotalCyclesSlowPath( cycleCountEntry* , timing& );
	static inline void logTiming( cycleCountEntry* , fastPathEntry* , timing& );
}

inline std::ostream& operator<<(std::ostream& os, const G4FastPathHadronicCrossSection::fastPathEntry& fp);

//Implementation of inline functions. Note the ifdef

namespace G4FastPathHadronicCrossSection {

#ifdef FPDEBUG
inline void logInvocationTriedOneLine(cycleCountEntry* cl ) {
	if ( cl != nullptr ) ++(cl->invocationCountTriedOneLineCache);
}
inline void logInvocationOneLine( cycleCountEntry* cl ) {
	if ( cl != nullptr ) ++(cl->invocationCountOneLineCache);
}
inline void logHit(cycleCountEntry* cl) {
	if ( cl != nullptr ) ++(cl->cacheHitCount);
}
inline void logInvocationCountFastPath( cycleCountEntry* cl )
{
	if ( cl != nullptr ) ++(cl->invocationCountFastPath);
}
inline void logInvocationCountSlowPAth( cycleCountEntry* cl)
{
	if ( cl != nullptr ) ++(cl->invocationCountSlowPath);
}

inline void logInitCyclesFastPath(cycleCountEntry* cl,timing& tm)
{
	if ( cl != nullptr ) cl->initCyclesFastPath = tm.rdtsc_stop - tm.rdtsc_start;
}
inline void logTotalCyclesFastPath( cycleCountEntry* cl,timing& tm)
{
	if ( cl!=nullptr ) cl->totalCyclesFastPath = tm.rdtsc_stop - tm.rdtsc_start;
}
inline void logTotalCyclesSlowPath( cycleCountEntry* cl,timing& tm)
{
	if ( cl!=nullptr ) cl->totalCyclesSlowPath = tm.rdtsc_stop - tm.rdtsc_start;
}
inline void logTiming( cycleCountEntry* entry , fastPathEntry* fast_entry, timing& timing)
{
	if (fast_entry != nullptr ) {
		if ( entry->invocationCountFastPath == 0 ) {
			//PRUTH style initialization
			G4FastPathHadronicCrossSection::logInitCyclesFastPath(entry,timing);
			G4FastPathHadronicCrossSection::logInvocationCountFastPath(entry);
		  } else {
			//PRUTH comment to understand:
			//the first one includes the initialization... don't count it for now
			G4FastPathHadronicCrossSection::logTotalCyclesFastPath(entry,timing);
			G4FastPathHadronicCrossSection::logInvocationCountFastPath(entry);
		  }
	  	 } else {
	  	   G4FastPathHadronicCrossSection::logInvocationCountSlowPAth(entry);
	  	   G4FastPathHadronicCrossSection::logTotalCyclesSlowPath(entry,timing);
	  	 }
}
#else
inline void logInvocationTriedOneLine(cycleCountEntry*){}
inline void logInvocationOneLine( cycleCountEntry*){}
inline void logHit(cycleCountEntry*){}
inline void logInvocationCountFastPath( cycleCountEntry*){}
inline void logInvocationCountSlowPAth( cycleCountEntry*){}
inline void logInitCyclesFastPath( cycleCountEntry* , timing& ){}
inline void logTotalCyclesFastPath( cycleCountEntry* , timing& ){}
inline void logTotalCyclesSlowPath( cycleCountEntry* , timing& ){}
inline void logTiming( cycleCountEntry* , fastPathEntry* , timing& ) {}
#endif

inline void getCrossSectionCount::MethodCalled() {
#ifdef FPDEBUG
	++methodCalled;
#endif
}

inline void getCrossSectionCount::HitOneLine() {
#ifdef FPDEBUG
	++hitOneLineCache;
#endif
}

inline void getCrossSectionCount::FastPath() {
#ifdef FPDEBUG
	++fastPath;
#endif
}

inline void getCrossSectionCount::SlowPath() {
#ifdef FPDEBUG
	++slowPath;
#endif
}

inline void getCrossSectionCount::SampleZandA() {
#ifdef FPDEBUG
	++sampleZandA;
#endif
}
}//namespace

inline std::ostream& operator<<(std::ostream& os, const G4FastPathHadronicCrossSection::fastPathEntry& fp) {
	using CLHEP::MeV;
	os<<"#Particle: "<<(fp.particle!=nullptr?fp.particle->GetParticleName():"UNDEFINED")<<"\n";
	os<<"#Material: "<<(fp.material!=nullptr?fp.material->GetName():"UNDEFINED")<<"\n";
	os<<"#min_cutoff(MeV): "<<fp.min_cutoff/MeV<<"\n";
#ifdef FPDEBUG
	os<<"#DEBUG COUNTERS: count="<<fp.count<<" slowpath_sum="<<fp.slowpath_sum<<" max_delta="<<fp.max_delta;
	os<<" min_delta="<<fp.min_delta<<" sum_delta="<<fp.sum_delta<<" sum_delta_square="<<fp.sum_delta_square<<"\n";
#endif
	os<<*(fp.physicsVector)<<"\n";
	return os;
}

#endif //G4FastPathHadronicCrossSection_hh
