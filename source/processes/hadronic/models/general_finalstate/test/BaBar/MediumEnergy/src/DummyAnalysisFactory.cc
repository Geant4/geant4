#include <AIDA/IAnalysisFactory.h>

#ifndef G4ANALYSIS_USE
#include "G4ios.hh"
extern "C" {
IAnalysisFactory* AIDA_createAnalysisFactory();
}

//////////////////////////////////////////////////////////////////////////////
IAnalysisFactory* AIDA_createAnalysisFactory(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  G4cout <<"Warning: Dummy analysis factory loaded" << G4endl;
  return 0;
}
#endif

