#ifndef Tst75ExecProcessLevel_H
#define Tst75ExecProcessLevel_H 1

#include "ExecProcessLevel.hh"

class Tst75ExecProcessLevel : public ExecProcessLevel
{

   public: 
   
      Tst75ExecProcessLevel( const TstReader* ); 
   
   protected:
      
      virtual void InitProcess( const TstReader* );


};

#endif
