///////////////////
//G4RayTrajectory.hh
///////////////////

#ifndef G4RayTrajectory_h
#define G4RayTrajectory_h 1

class G4Step;

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>
#include <g4rw/tpordvec.h>
#include "globals.hh"
#include "G4Track.hh"
#include "G4RayTrajectoryPoint.hh"


class G4RayTrajectory : public G4VTrajectory
{
   public:

   G4RayTrajectory(); 
   G4RayTrajectory(G4RayTrajectory & right);
   virtual ~G4RayTrajectory();

   inline void* operator new(size_t);
   inline void  operator delete(void*);
   inline int operator == (const G4RayTrajectory& right){return (this==&right);}

   virtual void AppendStep(const G4Step*);
   virtual void ShowTrajectory() const;
   virtual void DrawTrajectory(G4int i_mode=0) const {;}
   virtual int GetPointEntries() const {return positionRecord->entries();}
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
   { return (*positionRecord)[i]; };
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);  

   private:

   G4RWTPtrOrderedVector<G4RayTrajectoryPoint>* positionRecord;
};


extern G4Allocator<G4RayTrajectory> G4RayTrajectoryAllocator;

inline void* G4RayTrajectory::operator new(size_t)
{
   void* aTrajectory;
   aTrajectory = (void*)G4RayTrajectoryAllocator.MallocSingle();
   return aTrajectory;
}

inline void G4RayTrajectory::operator delete(void* aTrajectory)
{
   G4RayTrajectoryAllocator.FreeSingle((G4RayTrajectory*)aTrajectory);
}


#endif

