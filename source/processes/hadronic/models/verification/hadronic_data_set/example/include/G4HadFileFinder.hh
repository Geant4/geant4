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


#ifndef G4HadFileFinder_HH
#define G4HadFileFinder_HH 1

#include "globals.hh"


class G4HadFileFinder 
{
public:

  G4HadFileFinder();
  G4HadFileFinder(G4String&, G4String&);

  ~G4HadFileFinder();


private:

  void Crawldir(char*, G4String&);
  void G4StripFile(G4String&, char*, G4String&);

};
#endif
