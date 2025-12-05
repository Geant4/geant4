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
// Author: M.A. Cortes-Giraldo, Universidad de Sevilla
//
// History changelog prior creation of this example:
// - 13/04/2009: Messenger class added.
// - 17/10/2009: version 1.0
// - 20/11/2009: version 1.1 before publishing:
//   - Changed some names by more suitable ones
// - 02/08/2010: version 1.2-dev:
//   - Added possbility of applying axial symmetries
// - 14/09/2023: version 2.0
//   - Following Geant4 coding guidelines
// - 18/10/2025: version 3.0
//   - Creation of IAEASourceIdRegistry for thread-safe source_id assignation
//


#include "G4IAEAphspReader.hh"

#include "iaea_phsp.h"
#include "iaea_record.h"
#include "IAEASourceIdRegistry.hh"

#include <vector>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "Randomize.hh"
#include "G4Threading.hh"

#include "G4IAEAphspReaderMessenger.hh"


// =============================================================================

G4IAEAphspReader::G4IAEAphspReader(const char* filename, const G4int threads)
  :fVerbose(0)
{
  fTotalThreads = threads;
  fFileName = filename;

  InitializeMembers();
  InitializeSource(fFileName);
  //MAC  ReadAndStoreFirstParticle();
}


// =============================================================================

G4IAEAphspReader::G4IAEAphspReader(const G4String filename, const G4int threads)
  :fVerbose(0)
{
  fTotalThreads = threads;
  fFileName = filename;

  InitializeMembers();
  InitializeSource(fFileName);
  //MAC  ReadAndStoreFirstParticle();
}


// =============================================================================

G4IAEAphspReader::~G4IAEAphspReader()
{
  if (fVerbose > 0) G4cout << "Destroying G4IAEAphspReader" << G4endl;
  delete fMessenger;

  // clear and delete the vectors
  fParticleTypeVec->clear();
  fKinEVec->clear();
  fPosVec->clear();
  fMomDirVec->clear();
  fWeightVec->clear();

  if (fNumberOfExtraFloats > 0) {
    fExtraFloatVec->clear();
    fExtraFloatTypes->clear();
  }

  if (fNumberOfExtraInts > 0) {
    fExtraIntVec->clear();
    fExtraIntTypes->clear();
  }

  delete fParticleTypeVec;
  delete fKinEVec;
  delete fPosVec;
  delete fMomDirVec;
  delete fWeightVec;

  delete fExtraFloatVec;
  delete fExtraIntVec;

  delete fExtraFloatTypes;
  delete fExtraIntTypes;

  // IAEA file has to be closed
  const IAEA_I32 sourceRead = static_cast<IAEA_I32>( fSourceReadId );
  IAEA_I32 result = 0;

  if (fVerbose > 0)
    G4cout << "Destroying IAEA source #" << sourceRead << G4endl;

  iaea_destroy_source(&sourceRead, &result);
  if (result > 0) {
    IAEASourceIdRegistry::Instance().Release(fSourceReadId);
    G4cout << "File " << fFileName << ".IAEAphsp closed successfully!"
	   << G4endl;
  }
  else {
    G4Exception("G4IAEAphspReader::~G4IAEAphspReader()",
		"IAEAphspReader001", JustWarning,
		"IAEA file not closed properly");
  }

  if (fVerbose > 0) G4cout << "G4IAEAphspReader destroyed" << G4endl;
}


// =============================================================================

void G4IAEAphspReader::InitializeMembers()
{
  G4ThreeVector zeroVec;
  G4ThreeVector yAxis(0.0, 1.0, 0.0);
  G4ThreeVector zAxis(0.0, 0.0, 1.0);

  particle_time = 0.0;
  particle_position = zeroVec;

  if ( G4Threading::IsMultithreadedApplication() )
    fSourceReadId = G4Threading::G4GetThreadId();
  else
    fSourceReadId = 0;

  fOrigHistories = -1;
  fTotalParticles = -1;
  fExtraFloatTypes = new std::vector<G4int>;
  fExtraIntTypes = new std::vector<G4int>;

  fParticleTypeVec = new std::vector<G4int>;
  fKinEVec = new std::vector<G4double>;
  fPosVec = new std::vector<G4ThreeVector>;
  fMomDirVec = new std::vector<G4ThreeVector>;
  fWeightVec = new std::vector<G4double>;
  fExtraFloatVec = new std::vector< std::vector<G4double> >;
  fExtraIntVec = new std::vector< std::vector<G4long> >;

  fTotalParallelRuns = 1;
  fParallelRun = 1;
  fTimesRecycled = 0;
  fUsedOrigHistories = 0;
  fCurrentParticle = 0;
  fEndOfFile = false;
  fLastGenerated = true; //MAC to make the file 'restart' at the first event

  fGlobalPhspTranslation = zeroVec;
  fRotationOrder = 123;
  fAlpha = fBeta = fGamma = 0.;
  fIsocenterPosition = zeroVec;
  fCollimatorAngle = fGantryAngle = 0.;
  fCollimatorRotAxis = zAxis;
  fGantryRotAxis = yAxis;

  fAxialSymmetryX = false;
  fAxialSymmetryY = false;
  fAxialSymmetryZ = false;

  // Messenger class
  fMessenger = new G4IAEAphspReaderMessenger(this);
}


// =============================================================================

void G4IAEAphspReader::InitializeSource(const G4String filename)
{
  // Reserve a global ID and request it explicitly from IAEA source registry
  G4int reserved = IAEASourceIdRegistry::Instance().ReserveNextLowest();
  if (reserved < 0) {
    G4ExceptionDescription ed;
    ed << "No free IAEA source_ID's available" << G4endl;
    G4Exception("G4IAEAphspReader::InitializeSource",
		"IAEAphspReader002",FatalException, ed);
  }
  IAEA_I32 sourceRead = static_cast<IAEA_I32>(reserved);

  // Now try to open the IAEAphsp file, and check if all it's OK
  const IAEA_I32 accessRead = 1; //static_cast<IAEA_I32>(fAccessRead);
  IAEA_I32 result = 0;

  iaea_new_source(&sourceRead, const_cast<char*>(filename.data()),
		  &accessRead, &result, filename.size()+1);
  if ( sourceRead < 0 || result < 0 ) {
    IAEASourceIdRegistry::Instance().Release(reserved);
    G4ExceptionDescription msg;
    msg << "Could not open IAEA source file to read" << G4endl;
    G4Exception("G4IAEAphspReader::InitializeSource()",
		"IAEAphspReader003", FatalException, msg);
  }

  // Update fSourceReadId, just in case.
  fSourceReadId = static_cast<G4int>(sourceRead);   // should equal 'reserved'

  // After creating a new IAEA source, the value of sourceRead
  // may have change. We check this.
  G4cout << "G4IAEAphspReader ==> This object has IAEA source ID = "
	 << fSourceReadId << "." << G4endl;


  iaea_check_file_size_byte_order(&sourceRead, &result);
  if (result < 0 )
    G4Exception("G4IAEAphspReader::InitializeSource()",
		"IAEAphspReader004", FatalException,
		"Failure at iaea_check_size_byte_order()");


  // Now get the total number of particles stored in file

  IAEA_I32 particleType = -1; // To count all kind of particles
  IAEA_I64 nParticles;
  iaea_get_max_particles(&sourceRead, &particleType, &nParticles);

  if (nParticles < 0) {
    G4Exception("G4IAEAphspReader::InitializeSource()",
		"IAEAphspReader005", FatalException,
		"Failure at iaea_get_max_particles()");
  }
  else {
    fTotalParticles = static_cast<G4long>(nParticles);
    G4cout << "Total number of particles in file \"" << filename
	   << ".IAEAphsp\"" << " = " << GetTotalParticles() << G4endl;
  }


  // Now get the number of original simulated histories

  IAEA_I64 nHistories;
  iaea_get_total_original_particles(&sourceRead, &nHistories);
  fOrigHistories = static_cast<G4long>(nHistories);

  G4cout << "Number of original histories in header file \""
	 << filename << ".IAEAphsp\"" << " = " << nHistories << G4endl;


  // And finally, take all the information related to the extra variables

  IAEA_I32 nExtraFloat, nExtraInt;
  iaea_get_extra_numbers(&sourceRead, &nExtraFloat, &nExtraInt );
  fNumberOfExtraFloats = static_cast<G4int>( nExtraFloat );
  fNumberOfExtraInts = static_cast<G4int>( nExtraInt );

  G4cout << "The number of Extra Floats is " << fNumberOfExtraFloats
	 << " and the number of Extra Ints is " << fNumberOfExtraInts
	 << G4endl;

  // Initializing vectors

  IAEA_I32 extraFloatTypes[NUM_EXTRA_FLOAT], extraIntTypes[NUM_EXTRA_LONG];
  iaea_get_type_extra_variables(&sourceRead, &result,
				extraIntTypes, extraFloatTypes);

  if (fNumberOfExtraFloats > 0) {
    fExtraFloatTypes->reserve( fNumberOfExtraFloats );
    for (G4int jj = 0; jj < fNumberOfExtraFloats; jj++)
      fExtraFloatTypes->push_back( static_cast<G4int>(extraFloatTypes[jj]) );
  }

  if (fNumberOfExtraInts > 0) {
    fExtraIntTypes->reserve( fNumberOfExtraInts );
    for (G4int ii = 0; ii < fNumberOfExtraInts; ii++)
      fExtraIntTypes->push_back( static_cast<G4int>(extraIntTypes[ii]) );
  }
}


// =============================================================================
// This is the pure virtual method inherited from G4VUserPrimaryGeneratorAction
// and it is meant to be a template method.

void G4IAEAphspReader::GeneratePrimaryVertex(G4Event* evt)
{
  if (fLastGenerated) {
    RestartSourceFile();
    ReadAndStoreFirstParticle();
  }

  PrepareThisEvent();

  if (fNStat == 0) {
    ReadThisEvent();
    GeneratePrimaryParticles(evt);
  }
}


// =============================================================================

void G4IAEAphspReader::RestartSourceFile()
{
  // Restart counters and flags
  fCurrentParticle = 0;
  fEndOfFile = false;
  fLastGenerated = false;

  // Clear all the vectors
  fParticleTypeVec->clear();
  fKinEVec->clear();
  fPosVec->clear();
  fMomDirVec->clear();
  fWeightVec->clear();
  if (fNumberOfExtraFloats) fExtraFloatVec->clear();
  if (fNumberOfExtraInts) fExtraIntVec->clear();

  // Close and reopen the IAEA source with the same ID
  IAEA_I32 sourceRead = static_cast<IAEA_I32>(fSourceReadId);
  const IAEA_I32 accessRead = 1; // static_cast<IAEA_I32>( fAccessRead );
  IAEA_I32 result = 0;

  G4cout << "G4IAEAphspReader ==> Closing IAEA source ID = " << fSourceReadId
	 << " to re-open it and reset fOriginalHistoriesUsed"
	 << G4endl;
  iaea_destroy_source(&sourceRead, &result);

  iaea_new_source(&sourceRead, const_cast<char*>( fFileName.data() ),
		  &accessRead, &result, fFileName.size()+1 );
  if ( sourceRead < 0 || result < 0 )
    G4Exception("G4IAEAphspReader::RestartSourceFile()",
		"IAEAphspReader006", FatalException,
		"Could not open IAEA source file to read");

  // Ensure that source Id data member is correct.
  fSourceReadId = static_cast<G4int>(sourceRead);
  G4cout << "G4IAEAphspReader ==> IAEA source re-started with ID = "
	 << fSourceReadId << "." << G4endl;

  iaea_check_file_size_byte_order(&sourceRead, &result);

  if (result < 0 )
    G4Exception("G4IAEAphspReader::RestartSourceFile()",
		"IAEAphspReader007", FatalException,
		"Failure at iaea_check_size_byte_order()");
}



// =============================================================================

void G4IAEAphspReader::ReadAndStoreFirstParticle()
{
  // Particle properties
  IAEA_I32 type, nStat;
  IAEA_Float E, wt, x, y, z, u, v, w;
  IAEA_Float extraFloats[NUM_EXTRA_FLOAT];
  IAEA_I32 extraInts[NUM_EXTRA_LONG]; 

  IAEA_I32 sourceRead = static_cast<IAEA_I32>(fSourceReadId);

  // -------------------------------------------------------------------
  // Compute first and last particle in case that parallel run commands
  // may have not been issued
  // -------------------------------------------------------------------

  if (fParallelRun == 1)
    ComputeFirstLastParticle();

  // ------------------------------
  //  Go to the suitable particle
  // ------------------------------
  IAEA_I32 chunk = static_cast<IAEA_I32>((fParallelRun-1)*fTotalThreads);

  if ( G4Threading::IsMultithreadedApplication() )
    chunk += static_cast<IAEA_I32>(G4Threading::G4GetThreadId()+1);
  else
    chunk += 1;

  // G4cout << "DEBUG!! fParallelRun= " << fParallelRun << G4endl;
  // G4cout << "DEBUG!! fTotalThreads= " << fTotalThreads << G4endl;
  // G4cout << "DEBUG!! chunk= " << chunk << G4endl;

  IAEA_I32 totalChunks =
    static_cast<IAEA_I32>(fTotalParallelRuns*fTotalThreads);
  // G4cout << "DEBUG!! totalChunks= " << totalChunks << G4endl;

  IAEA_I32 result;
  iaea_set_parallel(&sourceRead, 0, &chunk, &totalChunks, &result);
  // G4cout << "DEBUG!!!   " << sourceRead << "  " << chunk
  //  	 << "   " << totalChunks << "    " << result << G4endl;

  if (result < 0) {
    G4ExceptionDescription ed;
    ed << "ERROR placing the cursor within the phsp file [iaea_set_parallel()]"
       << G4endl;
    G4Exception("G4IAEAphspReader::ReadAndStoreFirstParticle()",
		"IAEAphspReader008", FatalException, ed);
  }

  fCurrentParticle = fFirstParticle;  // To keep track of ordering in phsp file

  // -------------------------------------------------
  //  Obtain all the information needed from the file
  // -------------------------------------------------

  // read IAEA particle
  fCurrentParticle++; // to keep the position of this particle in the PSF
  iaea_get_particle(&sourceRead, &nStat, &type, 
		    &E, &wt, &x, &y, &z, &u, &v, &w, extraFloats, extraInts);
  // This function increases the number returne
  // by iaea_get_used_original_particles

  if (fVerbose > 1) {
    G4cout << std::setprecision(6) << G4endl
	   << "G4IAEAphspReader: Reading particle # "
	   << fCurrentParticle << "   type= " << type
	   << "   n_stat= " << nStat
	   << G4endl
	   << "\t\t E= " << E << "   wt= " << wt
	   << G4endl
	   << "\t\t x= " << x << "   y= " << y << "   z= " << z
	   << "   u= " << u << "   v= " << v << "   w= " << w
	   << G4endl;
  }

  if (nStat == -1)
    G4Exception("G4IAEAphspReader::ReadAndStoreFirstParticle()",
		"IAEAphspReader009", FatalException,
		"Cannot find source file");
  else if (nStat == 0)
    fNStat = 1; // needed to set this first particle as new event
  else
    fNStat = nStat;
    // important to calculate correlations properly

  // -------------------------------------------------
  //  Store the information into the data members
  // -------------------------------------------------

  fParticleTypeVec->push_back(static_cast<G4int>(type) );
  fKinEVec->push_back(static_cast<G4double>(E));

  G4ThreeVector pos, momDir;
  pos.set(static_cast<G4double>(x),
	  static_cast<G4double>(y),
	  static_cast<G4double>(z));

  momDir.set(static_cast<G4double>(u),
	     static_cast<G4double>(v),
	     static_cast<G4double>(w));

  fPosVec->push_back(pos);
  fMomDirVec->push_back(momDir);

  fWeightVec->push_back(static_cast<G4double>(wt));

  if (fNumberOfExtraFloats > 0) {
    std::vector<G4double> vExtraFloats;
    vExtraFloats.reserve(fNumberOfExtraFloats);
    for (G4int jj = 0; jj < fNumberOfExtraFloats; jj++)
      vExtraFloats.push_back( static_cast<G4double>(extraFloats[jj]) );

    fExtraFloatVec->push_back(vExtraFloats);
  }

  if (fNumberOfExtraInts > 0) {
    std::vector<G4long> vExtraInts;
    vExtraInts.reserve(fNumberOfExtraInts);
    for (G4int ii = 0; ii < fNumberOfExtraInts; ii++)
      vExtraInts.push_back( static_cast<G4long>(extraInts[ii]) );

    fExtraIntVec->push_back(vExtraInts);
  }

  //  Update fUsedOrigHistories
  // ---------------------------------

  IAEA_I64 nUsedOriginal;
  iaea_get_used_original_particles(&sourceRead, &nUsedOriginal);
  fUsedOrigHistories = static_cast<G4long>(nUsedOriginal);
  if (fVerbose > 0)
    G4cout << "G4IAEAphspReader::fUsedOrigHistories = "
	   << fUsedOrigHistories << " from IAEAphsp source #" << sourceRead
	   << G4endl;

  //  Safety check in case that only one particle is stored
  //  or only one particle can be taken into consideration.
  // -------------------------------------------------------

  if (fCurrentParticle == fLastParticle) 
    fEndOfFile = true;
}


// =============================================================================

void G4IAEAphspReader::PrepareThisEvent()
{
  fNStat--;  // A new event begins

  // Erase all the elements but the last in the vectors
  // In case that the size of the lists is just 1, there's no need
  // to clean up.

  G4int listSizes = fParticleTypeVec->size();

  if (listSizes > 1) {
    fParticleTypeVec->erase(fParticleTypeVec->begin(),
			    fParticleTypeVec->end()-1);
    fKinEVec->erase(fKinEVec->begin(), fKinEVec->end()-1);
    fPosVec->erase(fPosVec->begin(), fPosVec->end()-1);
    fMomDirVec->erase(fMomDirVec->begin(), fMomDirVec->end()-1);
    fWeightVec->erase(fWeightVec->begin(), fWeightVec->end()-1);

    if (fNumberOfExtraFloats > 0)
      fExtraFloatVec->erase(fExtraFloatVec->begin(), fExtraFloatVec->end()-1);

    if (fNumberOfExtraInts > 0)
      fExtraIntVec->erase(fExtraIntVec->begin(), fExtraIntVec->end()-1);
  }
}


// =============================================================================

void G4IAEAphspReader::ReadThisEvent()
{
  // Particle properties
  IAEA_I32 type, nStat;
  IAEA_Float E, wt, x, y, z, u, v, w;
  IAEA_I32 extraInts[NUM_EXTRA_LONG]; 
  IAEA_Float extraFloats[NUM_EXTRA_FLOAT];

  IAEA_I32 sourceRead = static_cast<IAEA_I32>(fSourceReadId);

  // -------------------------------------------------
  //  Obtain all the information needed from the file
  // -------------------------------------------------

  while (fNStat == 0 && !fEndOfFile) {
    //  Read next IAEA particle
    // -------------------------

    fCurrentParticle++;
    iaea_get_particle(&sourceRead, &nStat, &type, &E, &wt,
		      &x, &y, &z, &u, &v, &w, extraFloats, extraInts);
    if (fVerbose > 1)
      G4cout << std::setprecision(6) << G4endl
	     << "G4IAEAphspReader: Reading particle # "
	     << fCurrentParticle << "   type= " << type
	     << "   n_stat= " << nStat
	     << G4endl
	     << "\t\t E= " << E << "   wt= " << wt
	     << G4endl
	     << "\t\t x= " << x << "   y= " << y << "   z= " << z
	     << "   u= " << u << "   v= " << v << "   w= " << w
	     << G4endl;

    if (nStat == -1)
      G4Exception("G4IAEAphspReader::ReadThisEvent()",
		  "IAEAphspReader010", FatalException,
		  "Cannot find source file");
    else
      fNStat += nStat;  // statistical book-keeping


    //  Store the information into the data members
    // ---------------------------------------------

    fParticleTypeVec->push_back(static_cast<G4int>(type) );
    fKinEVec->push_back(static_cast<G4double>(E));

    G4ThreeVector pos, momDir;
    pos.set(static_cast<G4double>(x), 
	    static_cast<G4double>(y),
	    static_cast<G4double>(z));

    momDir.set(static_cast<G4double>(u),
	    static_cast<G4double>(v),
	    static_cast<G4double>(w));

    fPosVec->push_back(pos);
    fMomDirVec->push_back(momDir);

    fWeightVec->push_back(static_cast<G4double>(wt));

    if (fNumberOfExtraFloats > 0) {
      std::vector<G4double> vExtraFloats;
      vExtraFloats.reserve(fNumberOfExtraFloats);
      for (G4int jj = 0; jj < fNumberOfExtraFloats; jj++)
	vExtraFloats.push_back( static_cast<G4double>(extraFloats[jj]) );

      fExtraFloatVec->push_back(vExtraFloats);
    }

    if (fNumberOfExtraInts > 0) {
      std::vector<G4long> vExtraInts;
      vExtraInts.reserve(fNumberOfExtraInts);
      for (G4int ii = 0; ii < fNumberOfExtraInts; ii++)
	vExtraInts.push_back( static_cast<G4long>(extraInts[ii]) );

      fExtraIntVec->push_back(vExtraInts);
    }

    //  Update fUsedOrigHistories
    // ---------------------------------

    IAEA_I64 nUsedOriginal;
    iaea_get_used_original_particles(&sourceRead, &nUsedOriginal);
    fUsedOrigHistories = static_cast<G4long>(nUsedOriginal);
    if (fVerbose > 0)
      G4cout << "G4IAEAphspReader::fUsedOrigHistories = "
	     << fUsedOrigHistories << " from IAEAphsp source #" << sourceRead
	     << G4endl;

    //  Check whether the end of chunk has been reached
    // -------------------------------------------------
    if (fCurrentParticle == fLastParticle)
      fEndOfFile = true;
  }
}


// =============================================================================

void G4IAEAphspReader::GeneratePrimaryParticles(G4Event* evt)
{
  // -----------------------------------------------------------
  // Only when the EOF has been reached and fNStat is still 0 
  // all the particles must be read.
  // Otherwise, don't read the last particle.
  // -----------------------------------------------------------

  G4int listSize = fParticleTypeVec->size();

  if (fEndOfFile && fNStat == 0) {
    // Read all the particles, so this flag switches on
    fLastGenerated = true;

    // Throw a warning due to the following restarting
    G4cout << "WARNING: Last particle to read from IAEAphsp reached at event "
	   << evt->GetEventID() << "."
	   << G4endl;
  }
  else {
    // The last particle belongs to a later event
    listSize--;
  }

  // ---------------------------------------------------------------------
  // Generate the random rotations if any rotational symmetry is applied.
  // This is applied when recycling particles.
  // In order to keep correlations, the same rotation must be applied to
  // all the particles of the same original history during same re-use.
  // ---------------------------------------------------------------------

  std::vector<G4double> randomRotations;

  if (fAxialSymmetryZ || fAxialSymmetryX || fAxialSymmetryY) {
    randomRotations.reserve(fTimesRecycled+1);
    for (G4int ii = 0; ii <= fTimesRecycled; ii++) {
      G4double randomAngle = G4UniformRand();
      randomAngle *= (360.*deg);
      randomRotations.push_back(randomAngle);
    }
  }

  // ---------------------------------------------
  // loop over all the particles obtained from PSF
  // ---------------------------------------------
  for (G4int ii = 0; ii < listSize; ii++) {

    // First: Particle Definition
    // --------------------------
    
    G4ParticleDefinition * partDef = 0;
    switch((*fParticleTypeVec)[ii]) {
    case 1:
      partDef = G4Gamma::Definition();
      break;
    case 2:
      partDef = G4Electron::Definition();
      break;
    case 3:
      partDef = G4Positron::Definition();
      break;
    case 4:
      partDef = G4Neutron::Definition();
      break;
    case 5:
      partDef = G4Proton::Definition();
      break;
    default:
      G4ExceptionDescription ED;
      ED << "Particle code read at event #" << evt->GetEventID()
	 << "is" << (*fParticleTypeVec)[ii]
	 << " - this is not supported by the IAEAphsp format." << G4endl;
      G4Exception("G4IAEAphspReader::GeneratePrimaryParticles()",
		  "IAEAphspReader011", EventMustBeAborted, ED);
      return;
    }

    // Second: Particle position, time and momentum
    // --------------------------------------------

    particle_position = (*fPosVec)[ii];
    particle_position *= cm;        // IAEA file stores in cm

    G4double partMass = partDef->GetPDGMass();
    G4double partTotE = static_cast<G4double>((*fKinEVec)[ii])*MeV + partMass;
    G4double partMom = std::sqrt( partTotE*partTotE - partMass*partMass );

    G4ThreeVector partMomVec;
    partMomVec = (*fMomDirVec)[ii];
    partMomVec *= partMom;

    // Third: Translation and rotations
    // ---------------------------------

    // Translation is performed before rotations
    particle_position += fGlobalPhspTranslation;
    PerformRotations(partMomVec);

    // -------------------------------------------------
    //  Creation of the new primary particle and vertex
    // -------------------------------------------------

    // loop to take care of recycling
    for (G4int jj = 0; jj <= fTimesRecycled; jj++)  {
      // Apply the rotational symmetries if they are applicable
      if (fAxialSymmetryZ) {
	particle_position.rotateZ(randomRotations[jj]);
	partMomVec.rotateZ(randomRotations[jj]);
      }
      else if (fAxialSymmetryX) {
	particle_position.rotateX(randomRotations[jj]);
	partMomVec.rotateX(randomRotations[jj]);
      }
      else if (fAxialSymmetryY) {
	particle_position.rotateY(randomRotations[jj]);
	partMomVec.rotateY(randomRotations[jj]);
      }

      // Create the new primary particle
      G4PrimaryParticle * particle =
	new G4PrimaryParticle(partDef,
			      partMomVec.x(),partMomVec.y(),partMomVec.z());

      particle->SetWeight( ((*fWeightVec)[ii])/(fTimesRecycled+1) );

      // Create the new primary vertex and set the primary to it
      G4PrimaryVertex * vertex =
	new G4PrimaryVertex(particle_position, particle_time);
      vertex->SetPrimary(particle);

      if (fVerbose > 1)
	G4cout << std::setprecision(6) << G4endl
	       << "G4IAEAphspReader ==> Vertex produced: "
	       << " Event # " << evt->GetEventID()
	       << "   ParticleName: " << partDef->GetParticleName()
	       << G4endl
	       << "\t\t kinEnergy[MeV]= " << particle->GetKineticEnergy()/MeV
	       << "   weight= " << particle->GetWeight()
	       << G4endl
	       << "\t\t x[cm]= " << vertex->GetX0()/cm
	       << "   y[cm]= " << vertex->GetY0()/cm
	       << "   z[cm]= " << vertex->GetZ0()/cm
	       << "   u= " << particle->GetMomentumDirection().x()
	       << "   v= " << particle->GetMomentumDirection().y()
	       << "   w= " << particle->GetMomentumDirection().z()
	       << G4endl;

      // And finally set the vertex to this event
      evt->AddPrimaryVertex(vertex);
    }
  }
}


// =============================================================================

void G4IAEAphspReader::PerformRotations(G4ThreeVector & mom)
{
  PerformGlobalRotations(mom);
  PerformHeadRotations(mom);
}


// =============================================================================

void G4IAEAphspReader::PerformGlobalRotations(G4ThreeVector & momentum)
{
  switch (fRotationOrder) {
  case 123:
    particle_position.rotateX(fAlpha);
    particle_position.rotateY(fBeta);
    particle_position.rotateZ(fGamma);
    momentum.rotateX(fAlpha);
    momentum.rotateY(fBeta);
    momentum.rotateZ(fGamma);
    break;

  case 132:
    particle_position.rotateX(fAlpha);
    particle_position.rotateZ(fGamma);
    particle_position.rotateY(fBeta);
    momentum.rotateX(fAlpha);
    momentum.rotateZ(fGamma);
    momentum.rotateY(fBeta);
    break;

  case 213:
    particle_position.rotateY(fBeta);
    particle_position.rotateX(fAlpha);
    particle_position.rotateZ(fGamma);
    momentum.rotateY(fBeta);
    momentum.rotateX(fAlpha);
    momentum.rotateZ(fGamma);
    break;

  case 231:
    particle_position.rotateY(fBeta);
    particle_position.rotateZ(fGamma);
    particle_position.rotateX(fAlpha);
    momentum.rotateY(fBeta);
    momentum.rotateZ(fGamma);
    momentum.rotateX(fAlpha);
    break;

  case 312:
    particle_position.rotateZ(fGamma);
    particle_position.rotateX(fAlpha);
    particle_position.rotateY(fBeta);
    momentum.rotateZ(fGamma);
    momentum.rotateX(fAlpha);
    momentum.rotateY(fBeta);
    break;

  case 321:
    particle_position.rotateZ(fGamma);
    particle_position.rotateY(fBeta);
    particle_position.rotateX(fAlpha);
    momentum.rotateZ(fGamma);
    momentum.rotateY(fBeta);
    momentum.rotateX(fAlpha);
    break;

  default:
    G4ExceptionDescription ED;
    ED << "Invalid combination in G4IAEAphspReader::fRotationOrder -> "
       << "NO global rotations are done." << G4endl;
    G4Exception("G4IAEAphspReader::PerformGlobalRotations()",
		"IAEAphspReader012", FatalErrorInArgument, ED);
  }
}


// =============================================================================

void G4IAEAphspReader::PerformHeadRotations(G4ThreeVector & momentum)
{
  // First we need the position related to the isocenter
  particle_position -= fIsocenterPosition;

  // And now comes both rotations whose axis contain the isocenter

  // First is the rotation of the collimator around its defined axis
  particle_position.rotate(fCollimatorRotAxis, fCollimatorAngle);
  momentum.rotate(fCollimatorRotAxis, fCollimatorAngle);

  // The rotation of the complete machine is applied afterwards
  particle_position.rotate(fGantryRotAxis, fGantryAngle);
  momentum.rotate(fGantryRotAxis, fGantryAngle);

  // And finally, restore the position to the global coordinate system
  particle_position += fIsocenterPosition;
}



// =============================================================================

void G4IAEAphspReader::SetParallelRun(const G4int parallelRun)
{
  if (parallelRun > fTotalParallelRuns || parallelRun < 1) {
    G4ExceptionDescription ED;
    ED << "SetParallelRun argument must be an integer between 1 "
       << "and fTotalParallelRuns" << G4endl;
    G4Exception("G4IAEAphspReader::SetParallelRun()",
		"IAEAphspReader013", FatalErrorInArgument, ED);
  }
  fParallelRun = parallelRun;

  if (fVerbose > 0)
    G4cout << "G4IAEAphspReader::fParallelRun = " << fParallelRun
	   << G4endl;

  ComputeFirstLastParticle();
}


// =============================================================================

void G4IAEAphspReader::ComputeFirstLastParticle()
{
  // ---------------------------------------------------------
  // Get the limits within this instance reads from the file,
  // defined by first and last particle.
  // ---------------------------------------------------------

  G4cout << "G4IAEAphspReader: Defined " << fTotalParallelRuns
	 << " chunks in the file"
	 << " (i.e. number of independent simulations in parallel)" << G4endl;

  G4cout << "G4IAEAphspReader ==> This run is using the chunk #"
         << fParallelRun << " of the " << fTotalParallelRuns
         << " defined in the phsp file." << G4endl;

  // Below, operations are done exactly the same way as seen in
  // iaea_set_parallel() method, for consistency.
  // This means that all integers are obtained by truncation of floats.
  // Thus, it is expected that
  // particlesChunk*(fTotalParallelRuns*fTotalThreads) <= fTotalParticles

  IAEA_I64 particlesChunk = fTotalParticles/(fTotalParallelRuns*fTotalThreads);
  IAEA_I64 firstParticle = particlesChunk*(fParallelRun-1)*fTotalThreads;

  if ( G4Threading::IsMultithreadedApplication() ) {
    IAEA_I32 threadID = static_cast<IAEA_I32>(G4Threading::G4GetThreadId());
    firstParticle += particlesChunk*threadID;
  }

  IAEA_I64 lastParticle = firstParticle + particlesChunk;

  fFirstParticle = static_cast<G4long>(firstParticle);
  fLastParticle = static_cast<G4long>(lastParticle);

  // In case of roundings instead of truncation, the following may happen.
  if (fLastParticle > fTotalParticles)
    fLastParticle = fTotalParticles;

  // Note: The first particle is at position #1, not #0.
  G4cout << "G4IAEAphspReader: "
	 << "This thread is reading the particles from place #"
	 << fFirstParticle+1 << " to #" << fLastParticle
	 << " within the phsp file" << G4endl;
}


// =============================================================================

void G4IAEAphspReader::SetCollimatorRotationAxis(const G4ThreeVector & axis)
{
  G4double mod = axis.mag();
  if (mod > 0) {
    // Normalization
    G4double scale = 1.0/mod;
    G4double ux = scale*axis.getX();
    G4double uy = scale*axis.getY();
    G4double uz = scale*axis.getZ();
    fCollimatorRotAxis.setX(ux);
    fCollimatorRotAxis.setY(uy);
    fCollimatorRotAxis.setZ(uz);
  }
  else {
    G4ExceptionDescription ED;
    ED << "Trying to set as CollimatorRotationAxis a null vector, "
       << "thus the previous value will remain." << G4endl;
    G4Exception("G4IAEAphspReader::SetCollimatorRotationAxis()",
		"IAEAphspReader014", JustWarning, ED);
  } 
}


// =============================================================================

void G4IAEAphspReader::SetGantryRotationAxis(const G4ThreeVector & axis)
{
  G4double mod = axis.mag();
  if ( mod > 0 ) {
    // Normalization
    G4double scale = 1.0/mod;
    G4double ux = scale*axis.getX();
    G4double uy = scale*axis.getY();
    G4double uz = scale*axis.getZ();
    fGantryRotAxis.setX(ux);
    fGantryRotAxis.setY(uy);
    fGantryRotAxis.setZ(uz);
  }
  else {
    G4ExceptionDescription ED;
    ED << "Trying to set as GantryRotationAxis a null vector, "
       << "thus the previous value will remain." << G4endl;
    G4Exception("G4IAEAphspReader::SetGantryRotationAxis()",
		"IAEAphspReader015", JustWarning, ED);
  }
}


// =============================================================================

G4long G4IAEAphspReader::GetTotalParticlesOfType(const G4String type) const
{
  IAEA_I32 particleType;

  if (type == "ALL") particleType = -1;
  else if (type == "PHOTON") particleType = 1;
  else if (type == "ELECTRON") particleType = 2;
  else if (type == "POSITRON") particleType = 3;
  else if (type == "NEUTRON") particleType = 4;
  else if (type == "PROTON") particleType = 5;
  else particleType = 0;

  if (particleType == 0) {
    G4ExceptionDescription ED;
    ED << "Particle type \"" << type << "\" is not valid."
       << " (Please do ensure that all characters are uppercase)"
       << G4endl;
    G4Exception("G4IAEAphspReader::GetTotalParticlesOfType()",
		"IAEAphspReader016", JustWarning, ED);
    return -1;
  }

  IAEA_I64 nParticles;
  const IAEA_I32 sourceRead = static_cast<IAEA_I32>( fSourceReadId );
  iaea_get_max_particles(&sourceRead, &particleType, &nParticles);

  if (nParticles < 0)
    G4Exception("G4IAEAphspReader::GetTotalParticlesOfType()",
		"IAEAphspReader017", JustWarning,
		"IAEA header file not found");
  else
    G4cout << "G4IAEAphspReader: "
	   << "Total number of " << type << " in filename "
	   << fFileName << ".IAEAphsp" << " = " << nParticles << G4endl;

  return ( static_cast<G4long>(nParticles) );
}


// =============================================================================

G4double G4IAEAphspReader::GetConstantVariable(const G4int index) const
{
  const IAEA_I32 sourceRead = static_cast<IAEA_I32>( fSourceReadId );
  const IAEA_I32 idx = static_cast<IAEA_I32>( index );
  IAEA_Float constVariable;
  IAEA_I32 result;

  iaea_get_constant_variable(&sourceRead, &idx, &constVariable, &result);

  if (result == -1)
    G4Exception("G4IAEAphspReader::GetConstantVariable()",
		"IAEAphspReader018", JustWarning,
		"IAEA header file source not found");
  else if (result == -2)
    G4Exception("G4IAEAphspReader::GetConstantVariable()",
		"IAEAphspReader019", JustWarning,
		"index out of range");
  else if (result == -3)
    G4Exception("G4IAEAphspReader::GetConstantVariable()",
		"IAEAphspReader020", JustWarning,
		"parameter indicated by 'index' is not a constant");

  return ( static_cast<G4double>(constVariable) );
}
