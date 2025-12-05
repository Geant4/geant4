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

#ifndef G4IAEAphspReader_h
#define G4IAEAphspReader_h 1

#include "G4VPrimaryGenerator.hh"

#include <vector>

#include "globals.hh"
#include "G4ThreeVector.hh"


class G4Event;
class G4IAEAphspReaderMessenger;


class G4IAEAphspReader :public G4VPrimaryGenerator
{

public:

  G4IAEAphspReader(const char* filename, const G4int threads = 1);
  G4IAEAphspReader(const G4String filename, const G4int threads = 1);
  // 'filename' must include the path if needed, but NOT the extension
  ~G4IAEAphspReader() override;
  
  void GeneratePrimaryVertex(G4Event* evt) override;   // Mandatory

  inline void SetVerbose(const G4int verb)
  {
    fVerbose = verb;
    if (fVerbose > 0)
      G4cout << "G4IAEAphspReader::fVerbose = " << fVerbose << G4endl;
  }

  inline void SetTotalParallelRuns(const G4int nParallelRuns)
  {
    fTotalParallelRuns = nParallelRuns;
    if (fVerbose > 0)
      G4cout << "G4IAEAphspReader::fTotalParallelRuns = " << fTotalParallelRuns
	     << G4endl;
  }

  void SetParallelRun(const G4int parallelRun);
  inline void SetTotalThreads(const G4int threads) {fTotalThreads = threads;}
  inline void SetTimesRecycled(const G4int ntimes) {fTimesRecycled = ntimes;}

  inline void SetGlobalPhspTranslation(const G4ThreeVector & pos)
  {fGlobalPhspTranslation = pos;}
  inline void SetRotationOrder(const G4int ord)  { fRotationOrder = ord; }
  inline void SetRotationX(const G4double alpha) { fAlpha = alpha; }
  inline void SetRotationY(const G4double beta)  { fBeta = beta; }
  inline void SetRotationZ(const G4double gamma) { fGamma = gamma; }
  inline void SetIsocenterPosition(const G4ThreeVector & pos)
  {fIsocenterPosition = pos;}
  void SetCollimatorRotationAxis(const G4ThreeVector & axis);
  void SetGantryRotationAxis(const G4ThreeVector & axis);
  inline void SetCollimatorAngle(const G4double ang) {fCollimatorAngle = ang;}
  inline void SetGantryAngle(const G4double ang) {fGantryAngle = ang;}

  inline void SetAxialSymmetryX(const G4bool value) 
  {
    fAxialSymmetryX = value;
    if (value) {
      fAxialSymmetryY = false;
      fAxialSymmetryZ = false;
    }
  }
  inline void SetAxialSymmetryY(const G4bool value)
  {
    fAxialSymmetryY = value;
    if (value) {
      fAxialSymmetryZ = false;
      fAxialSymmetryX = false;
    }
  }
  inline void SetAxialSymmetryZ(const G4bool value)
  {
    fAxialSymmetryZ = value;
    if (value) {
      fAxialSymmetryX = false;
      fAxialSymmetryY = false;
    }
  }

  inline G4String GetFileName() const         {return fFileName;}
  inline G4int GetSourceReadId() const        {return fSourceReadId;}
  inline G4long GetOrigHistories() const      {return fOrigHistories;}
  inline G4long GetUsedOrigHistories() const  {return fUsedOrigHistories;}
  inline G4long GetTotalParticles() const     {return fTotalParticles;}
  inline G4int GetNumberOfExtraFloats() const {return fNumberOfExtraFloats;}
  inline G4int GetNumberOfExtraInts() const   {return fNumberOfExtraInts;}
  inline std::vector<G4int>* GetExtraFloatTypes() const
  {return fExtraFloatTypes;}
  inline std::vector<G4int>* GetExtraIntTypes() const
  {return fExtraIntTypes;}
  G4long GetTotalParticlesOfType(const G4String type) const;
  G4double GetConstantVariable(const G4int index) const;

  inline std::vector<G4int>* GetParticleTypeVec() const
  {return fParticleTypeVec;}
  inline std::vector<G4double>* GetKinEVec() const
  {return fKinEVec;}
  inline std::vector<G4ThreeVector>* GetPosVec() const
  {return fPosVec;}
  inline std::vector<G4ThreeVector>* GetMomDirVec() const
  {return fMomDirVec;}
  inline std::vector<G4double>* GetWeightVec() const
  {return fWeightVec;}
  inline std::vector< std::vector<G4double> >* GetExtraFloatVec() const
  {return fExtraFloatVec;}
  inline std::vector< std::vector<G4long> >* GetExtraIntVec() const
  {return fExtraIntVec;}

  inline G4int GetTotalParallelRuns() const {return fTotalParallelRuns;}
  inline G4int GetParallelRun() const       {return fParallelRun;}
  inline G4int GetTotalThreads() const      {return fTotalThreads;}
  inline G4long GetFirstParticle() const    {return fFirstParticle;}
  inline G4long GetLastParticle() const     {return fLastParticle;}
  inline G4int GetTimesRecycled() const     {return fTimesRecycled;}

  inline G4ThreeVector GetGlobalPhspTranslation() const
  {return fGlobalPhspTranslation;}
  inline G4int GetRotationOrder() const {return fRotationOrder;}
  inline G4double GetRotationX() const {return fAlpha;}
  inline G4double GetRotationY() const {return fBeta;}
  inline G4double GetRotationZ() const {return fGamma;}
  inline G4ThreeVector GetIsocenterPosition() const
  {return fIsocenterPosition;}
  inline G4double GetCollimatorAngle() const {return fCollimatorAngle;}
  inline G4double GetGantryAngle() const {return fGantryAngle;}
  inline G4ThreeVector GetCollimatorRotationAxis() const
  {return fCollimatorRotAxis;}
  inline G4ThreeVector GetGantryRotationAxis() const {return fGantryRotAxis;}

  inline G4bool GetAxialSymmetryX() const {return fAxialSymmetryX;}
  inline G4bool GetAxialSymmetryY() const {return fAxialSymmetryY;}
  inline G4bool GetAxialSymmetryZ() const {return fAxialSymmetryZ;}


private:

  G4IAEAphspReader() = default;

  void InitializeMembers();
  void InitializeSource(const G4String filename);
  void ComputeFirstLastParticle();
  void ReadAndStoreFirstParticle();
  void PrepareThisEvent();
  void ReadThisEvent();
  void GeneratePrimaryParticles(G4Event* evt);
  void PerformRotations(G4ThreeVector& mom);
  void PerformGlobalRotations(G4ThreeVector& mom);
  void PerformHeadRotations(G4ThreeVector& mom);
  void RestartSourceFile();


  // ========== Data members ==========

private:

  // ----------------------
  // FILE GLOBAL PROPERTIES
  // ----------------------

  G4String fFileName;
  // Must include the path, but NOT the IAEA extension

  G4int fSourceReadId;
  // The Id the file source has for the IAEA routines.
  // This value is set by IAEA routines, but should correspond to thread Id.

  // static const G4int fAccessRead = 1;
  // A value needed to open the file in the IAEA codes

  G4long fOrigHistories;
  // Number of original histories which generated the phase space file

  G4long fTotalParticles;
  // Number of particles stored in the phase space file

  G4int fNumberOfExtraFloats, fNumberOfExtraInts;
  // Number of extra variables stored for each particle

  std::vector<G4int>* fExtraFloatTypes; 
  std::vector<G4int>* fExtraIntTypes;
  // Identification to classify the different extra variables

  // ---------------------
  // PARTICLE PROPERTIES
  // ---------------------

  std::vector<G4int>* fParticleTypeVec;
  std::vector<G4double>* fKinEVec;
  std::vector<G4ThreeVector>* fPosVec;
  std::vector<G4ThreeVector>* fMomDirVec;
  std::vector<G4double>* fWeightVec;
  std::vector< std::vector<G4double> >* fExtraFloatVec;
  std::vector< std::vector<G4long> >* fExtraIntVec;

  // -------------------
  // COUNTERS AND FLAGS
  // -------------------

  G4int fTotalParallelRuns;
  // For independent parallel runs, number of fragments in which the
  // PSF is divided.

  G4int fParallelRun;
  // Sets the fragment of PSF from which the particles must be read.

  G4int fTotalThreads;
  // Stores the total number of threads being used. Set via G4RunManager.

  G4long fFirstParticle;
  // First particle to read.
  // Value given by the number of independent parallel runs and threads.

  G4long fLastParticle;
  // Last particle to read.
  // Value given by the number of independent parallel runs and threads.

  G4int fTimesRecycled;
  // Set the number of times that each particle is recycled (not repeated)

  G4int fNStat;
  // Decides how many events should pass before throwing a new particle

  G4long fUsedOrigHistories;
  // Variable that stores the number of original histories read so far

  G4long fCurrentParticle;
  // Number to store the current particle position in PSF

  G4bool fEndOfFile;
  // Flag active when the file has reached the end

  G4bool fLastGenerated;
  // Flag active only when the last particle has been simulated

  // ------------------------
  // SPATIAL TRANSFORMATIONS
  // ------------------------

  G4ThreeVector fGlobalPhspTranslation;
  // Global translation performed to particles

  G4int fRotationOrder;
  // Variable to decide first, second and third rotations
  // For example, 132 means rotations using X, Z and Y global axis

  G4double fAlpha, fBeta, fGamma;
  // Angles of rotations around global axis

  G4ThreeVector fIsocenterPosition;
  // Position of the isocenter if needed

  G4double fCollimatorAngle, fGantryAngle;
  G4ThreeVector fCollimatorRotAxis, fGantryRotAxis;
  // Angles and axis of isocentric rotations in the machine
  // The collimator ALWAYS rotates first.

  // --------------------
  // ROTATIONAL SYMMETRY
  // --------------------

  // Boolean data members to apply rotational symmetry around XYZ axis
  // Only one can be set to true.
  G4bool fAxialSymmetryX;
  G4bool fAxialSymmetryY;
  G4bool fAxialSymmetryZ;

  // ----------------
  // MESSENGER CLASS
  // ----------------

  G4IAEAphspReaderMessenger* fMessenger;

  // ----------
  // VERBOSITY
  // ----------
  G4int fVerbose;

};

#endif
