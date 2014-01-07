#ifndef Test19ExecProcessLevel_H
#define Test19ExecProcessLevel_H 1

#include "ExecProcessLevel.hh"

class Test19ExecProcessLevel : public ExecProcessLevel
{

   public: 
   
      Test19ExecProcessLevel( const TstReader* pset ) : ExecProcessLevel(pset) { InitProcess(pset); InitSetup(pset); InitBeam(pset);}
   
   protected:
      
      virtual void InitProcess( const TstReader* );


};

#endif
