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


#ifndef G4ReadHadDiffMomentumCXSC_HH
#define G4ReadHadDiffMomentumCXSC_HH 1

#include "G4HadDataReading.hh"
#include "globals.hh"
#include <vector>
#include <map>

class G4PhysicsTable;
class G4PhysicsVector;
class G4DataVector;
class G4Isotope;
class G4Element;
class G4Material;



class G4ReadHadDiffMomentumCXSC: public G4HadDataReading
{
public:

  G4ReadHadDiffMomentumCXSC(G4String&, G4Isotope*,G4String&);
  G4ReadHadDiffMomentumCXSC(G4String&, G4Element*,G4String&);
  G4ReadHadDiffMomentumCXSC(G4String&, G4Material*,G4String&);

  ~G4ReadHadDiffMomentumCXSC();


protected:

  void ReadDiffMomentumCXSC(G4String&, G4Isotope*,G4String&);
  void ReadDiffMomentumCXSC(G4String&, G4Element*,G4String&);
  void ReadDiffMomentumCXSC(G4String&, G4Material*,G4String&);


private:

};

#endif





