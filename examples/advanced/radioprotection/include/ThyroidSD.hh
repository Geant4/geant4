//
//    ********************************
//    *                              *  
//    *     ThyroidSD.hh             *
//    *                              *
//    ********************************

#ifndef ThyroidSD_h
#define ThyroidSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ThyroidHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ThyroidSD : public G4VSensitiveDetector
{
 public:
      ThyroidSD(G4String name);
      ~ThyroidSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();
 private:
     ThyroidHitsCollection *m_pRightThyroidHitsCollection;
	
    
};
#endif









