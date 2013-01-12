//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#ifndef G4ecpssrBaseKxsModel_hh
#define G4ecpssrBaseKxsModel_hh 1

#include "G4VecpssrKModel.hh"
#include "globals.hh"
#include <map>
#include <vector>

#include "G4CrossSectionDataSet.hh"


class G4ecpssrBaseKxsModel : public G4VecpssrKModel
{
public:

  G4ecpssrBaseKxsModel();

  ~G4ecpssrBaseKxsModel();
			     
  
  G4double CalculateCrossSection(G4int, G4double, G4double);//according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)
  
  G4double  ExpIntFunction(G4int n,G4double x);//Exponential Integral Function
  
private:
  
  G4ecpssrBaseKxsModel(const G4ecpssrBaseKxsModel&);
  G4ecpssrBaseKxsModel & operator = (const G4ecpssrBaseKxsModel &right);

  G4double FunctionFK(G4double k, G4double theta);

  G4double LogLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
   
  G4double LinLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
   
  G4double QuadInterpolator(G4double e11, 
 		            G4double e12, 
			    G4double e21, 
			    G4double e22, 
			    G4double x11,
			    G4double x12, 
			    G4double x21, 
			    G4double x22, 
			    G4double t1, 
			    G4double t2, 
			    G4double t, 
			    G4double e);

  typedef std::map<double, std::map<double, double> > TriDimensionMap;

  TriDimensionMap FKData;
  std::vector<double> dummyVec;

  typedef std::map<double, std::vector<double> > VecMap;
  VecMap aVecMap;

  G4int verboseLevel;

  G4CrossSectionDataSet* tableC1;
  G4CrossSectionDataSet* tableC2;
  G4CrossSectionDataSet* tableC3;
};
  
#endif
