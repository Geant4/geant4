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


#ifndef G4ReadHadInelasticXSC_HH
#define G4ReadHadInelasticXSC_HH 1

#include "G4HadFileSpec.hh"
#include "G4HadDataReading.hh"
#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"


class G4DataVector;
class G4Isotope;
class G4Element;
class G4Material;



class G4ReadHadInelasticXSC: public G4HadDataReading
{
public:

  G4ReadHadInelasticXSC(G4HadFileSpec&);

  ~G4ReadHadInelasticXSC();

protected:

  void ReadInelasticXSC(G4HadFileSpec&);
 
private:

};

#endif





