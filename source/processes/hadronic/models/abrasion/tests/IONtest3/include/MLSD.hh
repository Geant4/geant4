#ifndef MLSD_h
#define MLSD_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

//#include "MLGeometryConstruction.hh"
#include "G4Step.hh"
#include "MLHit.hh"

class MLGeometryConstruction;
class G4HCofThisEvent;
////////////////////////////////////////////////////////////////////////////////
//
class MLSD : public G4VSensitiveDetector
{
  public:
  
      MLSD (G4String, MLGeometryConstruction* );
     ~MLSD ();

      void Initialize (G4HCofThisEvent*);
      G4bool ProcessHits (G4Step*,G4TouchableHistory*);
      void EndOfEvent (G4HCofThisEvent*) {};
      void clear () {};
      void DrawAll () {};
      void PrintAll () {};

  private:
      MLHitsCollection         *MLCollection;
      MLGeometryConstruction   *geometry;

  friend class MLSteppingAction;
};
////////////////////////////////////////////////////////////////////////////////
#endif



