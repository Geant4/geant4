//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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
