#ifndef FTFPWrapper_HH
#define FTFPWrapper_HH 1

#include "ProcessWrapper.hh"

class FTFPWrapper : public ProcessWrapper
{

   public:
      
      FTFPWrapper(const G4String& name="FTFPWrapper",
                        G4ProcessType processType = fHadronic );
      ~FTFPWrapper() {}
      
      virtual void Compose();

};



#endif
