////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
#include "MLVersion.hh"
#include "globals.hh"

#include "RPTofstream.hh"
////////////////////////////////////////////////////////////////////////////////
//
RPTofstream & operator << (RPTofstream &RPTFile, MLVersion &q)
{
  RPTFile <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<"MULTI-LAYER SHIELDING SIMULATION SOFTWARE (MULASSIS)" <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<"Code Development Information" <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<G4endl;

  RPTFile <<"MULASSIS has been developed by    : " <<q.authors <<G4endl;
  RPTFile <<"  with the assistance of          : " <<q.coauthors <<G4endl;
  RPTFile <<"MULASSIS development sponsored by : " <<q.sponsor <<G4endl;
  RPTFile <<"MULASSIS version                  : " <<q.version <<G4endl;
#ifdef SPENVIS
  RPTFile <<"SPENVIS version                   : " <<q.SPENVISver <<G4endl;
  RPTFile <<"MULASSIS implementation           : " <<q.implementation <<G4endl;
#endif
  RPTFile <<"File creation date & time         : " <<G4endl;
  RPTFile <<G4endl;

  return RPTFile;
}
#endif

