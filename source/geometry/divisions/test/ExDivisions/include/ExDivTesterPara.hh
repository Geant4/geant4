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
// Class description:
//
// Virtual mother class of the division testers for each solid

// History:
// 13.08.03 - P.Arce Initial version
// ********************************************************************

#ifndef ExDivTesterPara_H
#define ExDivTesterPara_H 1

class G4VSolid;
class G4VPhysicalVolume;
class HepTransform3D;

#include "ExVDivTester.hh"

class ExDivTesterPara : public ExVDivTester
{ 
  public:  

  ExDivTesterPara( PVType& pvtype, std::vector<G4String>& extraPars );
  virtual ~ExDivTesterPara(){};

  virtual void GenerateScanPoints();
  virtual void BuildParentSolids();
  virtual void BuildChildrenSolids();

};

#endif

