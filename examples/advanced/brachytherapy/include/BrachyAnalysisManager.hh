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

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"

class ITree;
class IHistogramFactory;
class IAnalysisFactory;
class ITupleFactory;
class ITuple;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class BrachyAnalysisManager
{
public:
  
  ~BrachyAnalysisManager();

  void book();
  
  void finish();

  static BrachyAnalysisManager* getInstance();
  
  void analyse(G4double,G4double,G4float);
  void  hist(G4double,G4double,G4float);
  void Spectrum(G4double);

private:

  G4double xx,zz;
  G4float  en; 
  G4double  x,z;
  static BrachyAnalysisManager* instance;

private:
  BrachyAnalysisManager();

private:

  IAnalysisFactory  *aFact;
  ITree             *theTree;
  IHistogramFactory *histFact;
  ITupleFactory     *tupFact;

};

#endif




