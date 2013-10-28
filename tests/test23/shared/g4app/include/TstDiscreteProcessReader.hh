#ifndef TstProcessReader_H
#define TstProcessReader_H 1

#include "TstReader.hh"

class TstDiscreteProcessReader : public TstReader
{

   public:
   
      // ctor & dtor
      TstDiscreteProcessReader() : TstReader() {}
      virtual ~TstDiscreteProcessReader() {}
      
      // G4String        GetProcessName() const { return fProcessName; }

      virtual G4bool IsDiscreteProcess() const { return true; }
      virtual G4bool IsAtRestProcess()   const { return false; }
      virtual G4bool IsPhysicsList()     const { return false; }
      
   protected:
   
      virtual void ProcessLine( G4String line );
   
   // private:
  
      // G4String fProcessName;

};

#endif
