#ifndef TBBHELPER_HH
#define TBBHELPER_HH

void tbbSlaveBuildGeometryAndPhysicsVector();
void tbbSlaveDestroyGeometryAndPhysicsVector();
void tbbdebugmsg( const char* file, int where, const char* msg);

#ifdef TBBDEBUG
#include <sstream>
#define TBBMSG( msg ) { \
  std::ostringstream os; \
  os << msg; \
  tbbdebugmsg(__FILE__,__LINE__,os.str().c_str());
  }
#else
#define TBBMSG( msg )
#endif
#endif //TBBHELPER_HH
