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


#ifndef G4EXFORHadData_HH
#define G4EXFORHadData_HH 1

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



class G4EXFORHadData: public G4VHadDataWriting
{
public:

  G4EXFORHadData(G4String&,G4HadFileSpec&);
 
  ~G4EXFORHadData();


protected:

  void WriteDataFile(G4HadFileSpec&);
  virtual  void FillDoubleDiffXSC(G4String&);


private:
  static G4String EXFORenergy;
  static G4String EXFORangle;
  static G4String EXFORddXsc;
};

#endif









