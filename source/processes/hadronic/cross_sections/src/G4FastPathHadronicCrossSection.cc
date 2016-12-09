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
#include "G4FastPathHadronicCrossSection.hh"
#include "G4ios.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4CrossSectionDataStore.hh"
#include <vector>
#if defined(WIN32)
  //Needed for M_LN10
  #define _USE_MATH_DEFINES // for C++
  #include <math.h>
#endif
#include <cmath>
#include <array>

#ifdef FPDEBUG
#define DBG( msg ) G4cout<< msg <<G4endl;
#define DUMP() G4cout<< <<G4endl;
#else
#define DBG(msg)
#define DUMP()
#endif

using namespace G4FastPathHadronicCrossSection;

//Utility functions used to perform fast-path calculations.
//See later for details
namespace {
	struct Point_t {
		double e;
		double xs;
	};
	int  simplify_function(G4double tolerance,
                       std::vector<Point_t> & raw_data,
                       std::vector<Point_t>  & simplified_data);
	void RemoveBias( std::vector <Point_t> &,
			std::vector <Point_t> &,
			std::vector <Point_t> &);
}

fastPathEntry::fastPathEntry(const G4ParticleDefinition* part, const G4Material* mat, G4double min) :
		particle(part),material(mat),min_cutoff(min),physicsVector(nullptr)
	{
	DBG("Initializing a fastPathEntry");
#ifdef FPDEBUG
	count = 0;
	slowpath_sum=0.;
	max_delta=0.;
	min_delta=0.;
	sum_delta=0.;
	sum_delta_square=0.;
#endif
}

fastPathEntry::~fastPathEntry()
{
	DBG("Deleting fastPathEntry");
	DBG("Dumping status for: "<<(particle?particle->GetParticleName():"PART_NONE")<<" "\
		  <<(material?material->GetName():"MAT_NONE")<<" min_cutoff:"<<min_cutoff<<" "\
		  <<" count:"<<count<<" slowpath_sum:"<<slowpath_sum<<" max_delta:"<<max_delta\
		  <<" min_delta"<<min_delta<<" sum_delta"<<sum_delta<<" sum_delta_squared:"<<sum_delta_square);
	delete physicsVector;
}

//namespace {
//	static inline G4double exp10(G4double x) {
//		return std::exp( M_LN10*x);
//	}
//}

void fastPathEntry::Initialize(G4CrossSectionDataStore* xsds)
{
	//Check this method is called when G4CrossSectionDataStore is in the correct state:
	//  FastPath is enabled and we are indeed initializing
	assert( xsds->GetFastPathControlFlags().useFastPathIfAvailable &&
			xsds->GetFastPathControlFlags().initializationPhase );
	using std::log10;
	std::vector<Point_t> data_in;
	const fastPathParameters& params = xsds->GetFastPathParameters();
	G4double xs;

	  //G4double max_query = params.queryMax;
	  //G4int count = sampleCount;
	  //G4double tol = dpTol;

	//Shift so max and min are >= 1.
	//Don't forget to shift back before computing XS
	G4double min = params.sampleMin;
	G4double max = params.sampleMax;
	G4double shift = 0.0;
	if(min < 1.0){
		shift = 1.0 - min;
	  }
	min += shift;
	max += shift;

	G4double log_max = std::log10(params.sampleMax);
	G4double log_min = std::log10(params.sampleMin);
	G4double log_step = (log_max-log_min)/(1.0*params.sampleCount);

	G4double max_xs = 0.0;

	//Utility particle to calculate XS, with 0 kin energy by default
	static const G4ThreeVector constDirection(0.,0.,1.);
	G4DynamicParticle* probingParticle = new G4DynamicParticle( particle , constDirection , 0 );

	//add the cutoff energy
	probingParticle->SetKineticEnergy(min_cutoff);
	//Sample cross-section
	xs = xsds->GetCrossSection(probingParticle,material);
	data_in.push_back({min_cutoff,xs});

	G4double currEnergy = 0.0;
	//log results
	auto exp10 = [](G4double x){ return std::exp( M_LN10*x); };
	for(G4double log_currEnergy = log_min; log_currEnergy < log_max; log_currEnergy += log_step){
		currEnergy = exp10(log_currEnergy) - shift;
		if (currEnergy <  min_cutoff) continue;
		probingParticle->SetKineticEnergy(currEnergy);
		xs=xsds->GetCrossSection(probingParticle,material);
	    //G4cout << "PRUTH: energy value " << currEnergy << ", XS value " << xs << G4endl;
	    if (xs > max_xs) max_xs = xs;
	    data_in.push_back({currEnergy,xs});
	} // --- end of loop i
	probingParticle->SetKineticEnergy(max-shift);
	xs = xsds->GetCrossSection(probingParticle,material);
	data_in.push_back({max-shift,xs});

	G4double tol = max_xs * 0.01;
	std::vector<Point_t> decimated_data;
	simplify_function(tol,  data_in,  decimated_data);
	std::vector<Point_t> debiased_data;
	RemoveBias( data_in,  decimated_data,  debiased_data);
	if ( physicsVector != nullptr ) delete physicsVector;
	physicsVector = new XSParam(decimated_data.size());
	G4int physicsVectorIndex = 0;
	for(size_t i = 0; i < decimated_data.size(); i++){
		physicsVector->PutValue(physicsVectorIndex++, decimated_data[i].e, decimated_data[i].xs);
	}
	//xsds->DumpFastPath(particle,material,G4cout);
}

cycleCountEntry::cycleCountEntry(const G4String& pname , const G4Material* mat) :
		particle(pname),material(mat),fastPath(nullptr),
		energy(-1.),crossSection(-1.)
{
	DBG("Initializing cache entry");
#ifdef FPDEBUG
	cacheHitCount = 0;
	initCyclesFastPath=0;
	invocationCountSlowPath=0;
	totalCyclesSlowPath=0;
	invocationCountFastPath=0;
	totalCyclesFastPath=0;
	invocationCountTriedOneLineCache=0;
	invocationCountOneLineCache=0;
#endif
}

cycleCountEntry::~cycleCountEntry()
{
	DBG("Deleting cache entry");
	DBG(particle<<" "<<material<<" ("<<(material?material->GetName():"MAT_NONE")<<") "<<" "\
			<<"fast path pointer:"<<fastPath<<" stored:"<<energy<<" "<<crossSection<<" "\
			<<cacheHitCount<<" "<<initCyclesFastPath<<" "<<invocationCountSlowPath<<" "\
			<<totalCyclesSlowPath<<" "<<invocationCountFastPath<<" "<<totalCyclesFastPath<<" "\
			<<invocationCountTriedOneLineCache<<" "<<invocationCountOneLineCache);
}

#ifdef FPDEBUG
namespace {
	static inline unsigned long long rdtsc() {
		unsigned hi=0,lo=0;
#if defined(__GNUC__) &&( defined(__i386__)|| defined(__x86_64__) )
		__asm__ __volatile__ ("rdtsc":"=a"(lo),"=d"(hi));
#endif
		return ((unsigned long long)lo) | ((unsigned long long)hi<<32 );
	}
}
void G4FastPathHadronicCrossSection::logStartCountCycles(timing& tm)
{
	tm.rdtsc_start=rdtsc();
}
void G4FastPathHadronicCrossSection::logStopCountCycles(timing& tm)
{
	tm.rdtsc_stop=rdtsc();
}
#endif
getCrossSectionCount::getCrossSectionCount() {
#ifdef FPDEBUG
	methodCalled = 0;
	hitOneLineCache=0;
	fastPath=0;
	slowPath=0;
	sampleZandA = 0;
#endif
}

namespace {
// Rob Fowler's simplify code

//  This is a curve simplification routine based on the Douglas-Peucker
//  algorithm.
//  Simplifying assumptions are that the input polyline is a piecewise
//  function with the x values monotonically increasing,  that the function
//  reaches an asymptote at the right (high energy) end.
//  Also, the correct error measure is the difference in y between the original
//  curve and the result.
//  In GEANT4 use, the assumption is that the calling program has identified
//  low- and high-energy cutoffs and that the vector passed in is restricted
//  to the region between the cutoffs.


// The raw_data vector comes in ordered left to right (small energy to large).
//  The simplified_data vector is initially empty.

//A.Dotti ( 16-July-2015): transform variable size C-array and use of size_t
// 						   to remove compilation warnings

int  simplify_function(G4double tolerance,
                       std::vector<Point_t> & raw_data,
                       std::vector<Point_t>  & simplified_data)
{

  int  gap_left, gap_right;  // indices of the current region

  G4double tolsq = tolerance*tolerance;  // Alternative to working with absolute values.

  std::vector<int> working_stack;
  //A stack of the points to the right of the current interval that
  // are known to be selected.

  gap_right = raw_data.size() - 1;  // index of the last element.

  gap_left = 0;

  DBG("First and last elements  " << gap_left <<"  " <<gap_right);

  simplified_data.push_back(raw_data[0]); //copy first element over.

  DBG("first point ( 0  "
		      <<simplified_data[0].e  <<",  "<<simplified_data[0].xs <<" )");

  working_stack.push_back(gap_right); // 0th element on the stack.

  while ( !working_stack.empty() )
    {  G4double   a, slope, delta;
      G4double deltasq_max= tolsq;
      int i_max;

      gap_right = working_stack.back();  //get current TOS
      i_max = gap_right;


      if ( (gap_left +1) < gap_right )  // At least three points in the range.
        {
          // co-efficients for the left to right affine line segment
          slope =  (raw_data[gap_right].xs - raw_data[gap_left].xs) /
            (raw_data[gap_right].e - raw_data[gap_left].e);
          a = raw_data[gap_left].xs - slope * raw_data[gap_left].e;

          for ( int i = gap_left +1; i <gap_right;  i++) {
            delta = raw_data[i].xs - a - slope * raw_data[i].e;
            if ( delta * delta > deltasq_max){
              deltasq_max = delta * delta;
              i_max = i;
            }
          }
        }  else {
        	DBG("      Less than 3 point interval at [ "<< gap_left <<", " <<gap_right<< " ]");
	}

      if(i_max < gap_right) {  //  Found a new point, push it on the stack
	working_stack.push_back(i_max);
	DBG("         pushing point " << i_max);
	gap_right = i_max;
      }
      else { // didn't find a new point betweek gap_left and gap_right.
	simplified_data.push_back(raw_data[gap_right]);
	DBG("inserting point ("
			   <<gap_right  <<",  "<<raw_data[gap_right].e <<", "
			   << raw_data[gap_right].xs <<" )");
	gap_left = gap_right;
	working_stack.pop_back();
	gap_right = working_stack.back();
	DBG("      new gap_right " << gap_right);
      }

    }
  DBG("Simplified curve size "<< simplified_data.size());
  return (simplified_data.size());
}

// Rob Fowler's debias code

//  This is a de-biasing routine applied after using a curve simplification
//  routine based on the Douglas-Peucker
//  algorithm.
//  Simplifying assumptions are that the input polyline is a piecewise
//  function with the x values monotonically increasing, and
//  The right error measure is the difference in y between the original
//  curve and the result.



void RemoveBias(std::vector<Point_t> &  original, std::vector<Point_t> & simplified,
		std::vector<Point_t> & result){


  const size_t originalSize = original.size();
  const size_t simplifiedSize = simplified.size();

  //Create index mapping array
  std::vector<G4int>  xindex(simplifiedSize,0);
  //G4int xindex[simplifiedSize];
  G4int lastmatch = 0;
  G4int  j = 0;

  DBG("   original and simplified vector sizes  " << originalSize <<"  "<<simplifiedSize);

  for (size_t k = 0; k <simplifiedSize; k++) {
    for (size_t i = lastmatch; i < originalSize; i++) {
      if (original[i].e == simplified[k].e) {
        xindex[j++] = i;
        lastmatch = i;
      }
    }
  }

  DBG("Matched  " << j << " values of the simplified vector");
  // Use short names here.
  G4int m = simplifiedSize;

  std::vector<G4double> GArea(m-1,0);
  //G4double GArea [m-1];
  G4double GAreatotal = 0;

  //Area of original simplified curve
  for(int i = 0; i < m-1; i++){
    G4double GAreatemp = 0;

    for(j = xindex[i]; j< xindex[i+1]; j++){
      G4double trap = (original[j+1].xs + original[j].xs) * (original[j+1].e - original[j].e)/2.0;
      GAreatemp = GAreatemp + trap;
    }
    GArea[i] = GAreatemp;
    GAreatotal = GAreatotal + GAreatemp;
  }

  DBG("   Area under the original curve " << GAreatotal);

  //aleph        Why is this not alpha?


  std::vector<G4double> aleph(m-1,0);
  //G4double aleph [m-1];
  for(int i = 0; i< m-1; i++){
    aleph[i] = (simplified[i+1].e - simplified[i].e)/2.0;
  }

  //solve for f
  std::vector<G4double> adjustedy(m-1,0);
  //G4double adjustedy [m];
  adjustedy[m-1] = simplified[m-1].xs;
  for(int i = 2; i < m+1; i++) {
    adjustedy[m-i] = (GArea[m-i]/aleph[m-i]) - adjustedy[m-i+1];
    if (adjustedy[m-i] <0.0) {
      adjustedy[m-i] = 0.0;
      DBG("   Fixing negative cross section at index " << (m-i));
    }
  }

  //error and difference tracking
  std::vector<G4double> difference(m,0.);
  //G4double difference [m];
  G4double maxdiff = 0;
  G4double adjustedarea = 0;
  G4double simplifiedarea = 0;
  for(int i = 0; i < m-1; i++){
    G4double trap;
    trap = (adjustedy[i+1]+adjustedy[i])*(simplified[i+1].e-simplified[i].e)/2.0;
    adjustedarea = adjustedarea+trap;
    trap = (simplified[i+1].xs+simplified[i].xs)*(simplified[i+1].e-simplified[i].e)/2.0;
    simplifiedarea = simplifiedarea + trap;
  }

  DBG("   Area:  Simplified curve = " <<simplifiedarea);
  DBG("    Area:  Debiased curve  = " << adjustedarea);

  for(int i = 0; i <m; i++) {
    difference[i] = simplified[i].xs-adjustedy[i];
  }
  for(int i = 0; i <m; i++){
    if(std::fabs(difference[i]) > maxdiff) {
      maxdiff = std::fabs(difference[i]);
    }
  }
  //  what is the significance of the loops above ?

  for(size_t i = 0; i < simplifiedSize; i++){
    result.push_back( {simplified[i].e , adjustedy[i] } );
  }

}


}
