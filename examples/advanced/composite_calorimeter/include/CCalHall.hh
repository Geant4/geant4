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
//
///////////////////////////////////////////////////////////////////////////////
// File: CCalHall.hh
// Description: Equipped to construct the geometry of the 96 Test Beam
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalHall_h
#define CCalHall_h 1

#include "CCalDetector.hh"

class CCalHall: public CCalDetector {
public:
  //Constructor and Destructor
  CCalHall(const G4String &name);
  virtual ~CCalHall();

  //Get Methods
  G4String getMaterial()                  const {return genMaterial;}
  double   getDy_2Hall()                  const {return dy_2Hall;}
  double   getDx_2Hall()                  const {return dx_2Hall;}

protected:
  virtual int readFile();
  virtual void constructDaughters();

private:
  G4String genMaterial;            //General material
  double   dy_2Hall;               //Half width     of the Experimental Hall
  double   dx_2Hall;               //Half thickness of the Experimental Hall
};

#endif
