#ifndef LXeScintSD_h
#define LXeScintSD_h 1

#include "LXeScintHit.hh"

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

class LXeScintSD : public G4VSensitiveDetector
{
public:
  LXeScintSD(G4String name);
  ~LXeScintSD();
  
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  void EndOfEvent(G4HCofThisEvent* HCE);
  void clear();
  void DrawAll();
  void PrintAll();
  
private:
  LXeScintHitsCollection* scintCollection;
  
};

#endif

