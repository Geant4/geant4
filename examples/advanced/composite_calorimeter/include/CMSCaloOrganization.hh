///////////////////////////////////////////////////////////////////////////////
// File: CMSCaloOrganization.hh
// Date: 29.10.99 
// Description: Definition of sensitive unit numbering schema
// Modifications
///////////////////////////////////////////////////////////////////////////////

#ifndef CMSCaloOrganization_h
#define CMSCaloOrganization_h

#include <map>

class CMSCaloOrganization{

public:
  CMSCaloOrganization(){};
  virtual ~CMSCaloOrganization(){};
	 
public:
  unsigned int packindex(int det, int z, int eta, int phi) const;
  unsigned int packindex(int det, int depth, int z, int eta, int phi) const;
  void unpackindex(const unsigned int& idx, int& det, int& z, int& eta, int& phi) const;
  void unpackindex(const unsigned int& idx, int& det, int& depth, int& z, int& eta, int& phi) const;
    
public:
      
  // additional tool:
  int getUnitWithMaxEnergy(map<int,float,less<int> >& themap);

  float  energyInMatrix(int nCellInEta, int nCellInPhi, 
			int crystalWithMaxEnergy, 
			map<int,float,less<int> >& themap); 


};

#endif
