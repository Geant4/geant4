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
//
// R&D: Vladimir.Grichine@cern.ch
//      Simone.Gilardoni@cern.ch
//  
//  Class for writing data from MSU data base:
//  http://srd.sinp.msu.ru/INTAS/cr4.html
//  the data are double differential cross-sections
//
//

#ifndef G4MSUHadData_HH
#define G4MSUHadData_HH 1

#include "G4VHadDataWriting.hh"
#include "globals.hh"
#include <vector>
#include <map>
#include "G4HadFileSpec.hh"

class G4PhysicsTable;
class G4PhysicsVector;
class G4DataVector;
class G4Isotope;
class G4Element;
class G4Material;



class G4MSUHadData: public G4VHadDataWriting
{
public:

  G4MSUHadData(G4String&,G4HadFileSpec&);
 
  ~G4MSUHadData();


protected:

  void WriteDataFile(G4HadFileSpec&);
  virtual  void FillDoubleDiffXSC(G4String&);


private:
  static G4String fMSUenergy;
  static G4String fMSUangle;
  static G4String fMSUddXsc;
};

#endif









