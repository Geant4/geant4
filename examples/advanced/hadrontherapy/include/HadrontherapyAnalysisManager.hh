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
#ifdef G4ANALYSIS_USE
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
# include <AIDA/AIDA.h>

namespace AIDA{
  class ITree; 
  class IDataPoint;
  class IAnalysisFactory;
  class ITreeFactory;
};

class HadrontherapyAnalysisManager
{
private:
  HadrontherapyAnalysisManager();

public:
  ~HadrontherapyAnalysisManager();
  static HadrontherapyAnalysisManager* getInstance();
  void book();

  void energyDeposit3D(G4int PointNumber, G4int voxelX, G4int voxelY, 
                       G4int voxelZ, G4double edep);  
//Store the energy deposit of the sensitive detector in a DataSet
  
  void finish();

private:
  static HadrontherapyAnalysisManager* instance;

private:
  AIDA::IAnalysisFactory*  aFact;
  AIDA::ITree*             theTree;
  AIDA::ITreeFactory      *treeFact;
  AIDA::IDataPointSetFactory *  dataPointFactory;  
  AIDA::IDataPointSet *  energyDepositDataPoint; 
};

#endif
#endif



