///////////////////////////////////////////////////////////////////////////////
// File: CCaloOrganization.hh
// Description: Packing, unpacking and other related utilities for 
//              calorimetric numbering schema
///////////////////////////////////////////////////////////////////////////////

#ifndef CCaloOrganization_h
#define CCaloOrganization_h

#include "g4std/map"
#include "globals.hh"

class CCaloOrganization{

public:
  CCaloOrganization(){};
  virtual ~CCaloOrganization(){};
	 
public:
  unsigned int packindex(int det, int z, int eta, int phi) const;
  unsigned int packindex(int det, int depth, int z, int eta, int phi) const;
  void unpackindex(const unsigned int& idx, int& det, int& z, int& eta, int& phi) const;
  void unpackindex(const unsigned int& idx, int& det, int& depth, int& z, int& eta, int& phi) const;
    
public:
      
  // additional tool:
  int getUnitWithMaxEnergy(G4std::map<int,float,G4std::less<int> >& themap);

  float  energyInMatrix(int nCellInEta, int nCellInPhi, 
			int crystalWithMaxEnergy, 
			G4std::map<int,float,G4std::less<int> >& themap); 


};

#endif
