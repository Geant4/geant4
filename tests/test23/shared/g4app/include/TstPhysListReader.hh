#ifndef TstPhysListReader_H
#define TstPhysListReader_H 1

#include "TstReader.hh"

class TstPhysListReader : public TstReader
{

   public:  
      // ctor & dtor

      TstPhysListReader() : TstReader() {}
      virtual ~TstPhysListReader() {}

      // G4String     GetPhysListName()   const { return fPhysListName; }

      virtual G4bool IsDiscreteProcess() const { return false; }
      virtual G4bool IsAtRestProcess()   const { return false; }
      virtual G4bool IsPhysicsList()     const { return true; }
   
   protected:
   
      virtual void ProcessLine( G4String line );

//   private:
//   
//      G4String fPhysListName;

};

#endif
