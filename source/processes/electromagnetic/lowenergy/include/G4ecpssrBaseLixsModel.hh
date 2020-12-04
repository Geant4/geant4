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
//
// Author: Haifa Ben Abdelouahed
//         
//
// History:
// -----------
//  23 Apr 2008   H. Ben Abdelouahed   1st implementation
//  28 Apr 2008   MGP        Major revision according to a design iteration
//  29 Apr 2009   ALF Updated Desing for Integration
//  01 Sep 2009   ALF Updated to G4AnalyticalEcpssrLiCrossSection
//  28 Oct 2011   ALF Changed name G4AnalyticalEcpssrLiCrossSection to G4ecpssrBaseLixsModel
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, Cross section, p and alpha ionisation, L shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------


#ifndef G4ecpssrBaseLixsModel_hh
#define G4ecpssrBaseLixsModel_hh 1

#include "G4VecpssrLiModel.hh"
#include "globals.hh"
#include <map>
#include <vector>

class G4ecpssrBaseLixsModel : public G4VecpssrLiModel

{
public:

  G4ecpssrBaseLixsModel();

  ~G4ecpssrBaseLixsModel();
			     
  G4double CalculateL1CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident);//according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)

  G4double CalculateL2CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident);//according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)

  G4double CalculateL3CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident);//according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)
				    
  G4double CalculateVelocity(G4int subShell, G4int zTarget,G4double massIncident, G4double energyIncident); 
			      
  G4double  ExpIntFunction(G4int n,G4double x);//Exponential Integral Function

   

private:


  G4ecpssrBaseLixsModel(const G4ecpssrBaseLixsModel&);
  G4ecpssrBaseLixsModel & operator = (const G4ecpssrBaseLixsModel &right);

  G4double FunctionFL1(G4double k, G4double theta);
  
  G4double FunctionFL2(G4double k, G4double theta);


  G4double LogLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
   
  G4double LinLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
 
  G4double LinLinInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
  
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

  TriDimensionMap FL1Data;
  
  TriDimensionMap FL2Data;
  std::vector<double> dummyVec1;
  std::vector<double> dummyVec2;



  typedef std::map<double, std::vector<double> > VecMap;
  VecMap aVecMap1;
  VecMap aVecMap2;

  G4int verboseLevel;

};

#endif
