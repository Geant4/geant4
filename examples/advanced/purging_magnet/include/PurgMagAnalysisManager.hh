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
// Code developed by:
//  S.Larsson
//
//    ************************************
//    *                                  *
//    *    PurgMagAnalysisManager.hh     *
//    *                                  *
//    ************************************
//
// $Id: PurgMagAnalysisManager.hh,v 1.2 2004/06/18 09:17:43 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#ifdef G4ANALYSIS_USE
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IAnalysisFactory.h"

namespace AIDA
{
  class ITree;
  class IHistogramFactory;
  class IAnalysisFactory;
  class ITupleFactory;
  class ITuple;
  class ITreeFactory;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class PurgMagAnalysisManager

{
public:
  
  ~PurgMagAnalysisManager();
  
  void book();
  
  void finish();
  
  static PurgMagAnalysisManager* getInstance();
  
  void fill_Tuple_Electrons(G4double,G4double,G4double,G4double,G4double,G4double,G4double);
  void fill_Tuple_Gamma(G4double,G4double,G4double,G4double,G4double,G4double,G4double);
  void fill_Tuple_Positrons(G4double,G4double,G4double,G4double,G4double,G4double,G4double);
  
private:
  
  G4double ex,ey,ez,ee,epx,epy,epz,gx,gy,gz,ge,gpx,gpy,gpz,px,py,pz,pe,ppx,ppy,ppz;
  
  static PurgMagAnalysisManager* instance;
  
private:
  PurgMagAnalysisManager();
  
private:
  
  AIDA::IAnalysisFactory*  aFact;
  AIDA::ITree*             theTree;
  AIDA::IHistogramFactory *histFact;
  AIDA::ITupleFactory     *tupFact;
  AIDA::ITreeFactory      *treeFact;
  AIDA::ITuple *ntuple1;
  AIDA::ITuple *ntuple2;
  AIDA::ITuple *ntuple3;
};

#endif
#endif



