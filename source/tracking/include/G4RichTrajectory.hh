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
// $Id: G4RichTrajectory.hh 90219 2015-05-21 08:07:27Z gcosmo $
//
//---------------------------------------------------------------
//
// G4RichTrajectory.hh
//
// class description:
//   This class extends G4Trajectory, which includes the following:
//     1) List of trajectory points which compose the trajectory,
//     2) static information of particle which generated the 
//        trajectory, 
//     3) trackID and parent particle ID of the trajectory.
//   The extended information, only publicly accessible through AttValues,
//   includes:
//     4) physical volume and next physical volume;
//     5) creator process;
//   ...and much more.
//
// Contact:
//   Questions and comments on G4Trajectory should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@kekvax.kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//   and on the extended code to:
//     John Allison   (e-mail: John.Allison@manchester.ac.uk)
//     Joseph Perl    (e-mail: perl@slac.stanford.edu)
//
// ---------------------------------------------------------------

#ifndef G4RICHTRAJECTORY_HH
#define G4RICHTRAJECTORY_HH

#include "trkgdefs.hh"
#include "G4Trajectory.hh"
#include "G4TouchableHandle.hh"

typedef std::vector<G4VTrajectoryPoint*>  RichTrajectoryPointsContainer;

class G4RichTrajectory : public G4Trajectory {
  
public: // with description
  
  // Constructors/destructors
  G4RichTrajectory();
  G4RichTrajectory(const G4Track* aTrack);
  G4RichTrajectory(G4RichTrajectory &);
  virtual ~G4RichTrajectory();

private:
  G4RichTrajectory& operator= (const G4RichTrajectory &);

  // Operators
public:
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  inline int operator == (const G4RichTrajectory& right) const
  {return (this==&right);} 
  
  // Other (virtual) member functions
  void ShowTrajectory(std::ostream& os=G4cout) const;
  void DrawTrajectory() const;
  void AppendStep(const G4Step* aStep);
  void MergeTrajectory(G4VTrajectory* secondTrajectory);
  int GetPointEntries() const { return fpRichPointsContainer->size(); }
  G4VTrajectoryPoint* GetPoint(G4int i) const 
   { return (*fpRichPointsContainer)[i]; }

  // Get methods for HepRep style attributes
  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;

private:

  // Extended information (only publicly accessible through AttValues)...
  RichTrajectoryPointsContainer* fpRichPointsContainer;
  G4TouchableHandle fpInitialVolume;
  G4TouchableHandle fpInitialNextVolume;
  const G4VProcess* fpCreatorProcess;
  G4int             fCreatorModelID;
  G4TouchableHandle fpFinalVolume;
  G4TouchableHandle fpFinalNextVolume;
  const G4VProcess* fpEndingProcess;
  G4double          fFinalKineticEnergy;

};

extern G4TRACKING_DLL G4ThreadLocal
G4Allocator<G4RichTrajectory> *aRichTrajectoryAllocator;

inline void* G4RichTrajectory::operator new(size_t)
{
  if (!aRichTrajectoryAllocator)
  { aRichTrajectoryAllocator = new G4Allocator<G4RichTrajectory>; }
  return (void*)aRichTrajectoryAllocator->MallocSingle();
}

inline void G4RichTrajectory::operator delete(void* aRichTrajectory)
{
  aRichTrajectoryAllocator->FreeSingle((G4RichTrajectory*)aRichTrajectory);
}

#endif
