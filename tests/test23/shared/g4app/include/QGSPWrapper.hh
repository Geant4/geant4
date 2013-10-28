#ifndef QGSPWrapper_HH
#define QGSPWrapper_HH 1

#include "ProcessWrapper.hh"

class QGSPWrapper : public ProcessWrapper
{

   public:
   
      QGSPWrapper( const G4String& name = "QGSPWrapper", 
                         G4ProcessType processType = fHadronic );
      ~QGSPWrapper() {}
      
      virtual void Compose();

};

#endif
