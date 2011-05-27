#ifndef Hits_h
#define Hits_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"



class Hits : public G4VHit
{
  public:

      Hits();
     ~Hits();
      Hits(const Hits&);
      const Hits& operator=(const Hits&);
      int operator==(const Hits&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  
  //funzioni che gestiscono le informazioni salvate nell'hit
  void SetTrackID  (G4int track)      { trackID = track; };  
  void SetEdep     (G4double de)      { edep = de; };
  void SetPos      (G4ThreeVector xyz){ pos = xyz; };
  void SetParentID (G4int id) {parentID = id;};
    
  G4int GetTrackID()    { return trackID; };
  G4double GetEdep()    { return edep; };      
  G4ThreeVector GetPos(){ return pos; };
  G4int GetParentID() {return parentID;}
      
  private:
  //parametri che costituiscono l'hit
  G4int         trackID;
  G4double      edep;
  G4ThreeVector pos;
  G4int parentID;
};


typedef G4THitsCollection<Hits> HitsCollection;
//HitsCollection rappresenta un nome convenzionale per indicare il template
//G4THitsCollection applicato alla classe <Hits> !!

extern G4Allocator<Hits> HitAllocator;


inline void* Hits::operator new(size_t)
{
  void *aHit;
  aHit = (void *) HitAllocator.MallocSingle();
  return aHit;
}


inline void Hits::operator delete(void *aHit)
{
  HitAllocator.FreeSingle((Hits*) aHit);
}


#endif
