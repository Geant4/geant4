// Rich advanced example for Geant4
// FilterTypeSpec.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef FilterTypeSpec_h
#define FilterTypeSpec_h 1
//  different types of filters.
static const G4int MaxNumberOfFilterTypes=1;
enum FilterType{GlassD263};
static const G4String FilterTypeString[]={"GlassD263"};
//
// in the G4Example only 1 type of filter used. In the
// LHCb implementation several types of filters used. 
#endif 
