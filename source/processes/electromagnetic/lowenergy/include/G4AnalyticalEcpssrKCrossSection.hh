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
// $Id: G4AnalyticalEcpssrKCrossSection.hh,v 1.1 2010-06-06 23:40:35 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4ANALYTICALECPSSRKCROSSSECTION_HH
#define G4ANALYTICALECPSSRKCROSSSECTION_HH 1

#include "G4VecpssrKModel.hh"
#include "globals.hh"
#include <map>
#include <vector>

#include "G4DNACrossSectionDataSet.hh"


class G4AnalyticalEcpssrKCrossSection : public G4VecpssrKModel
{
public:

  G4AnalyticalEcpssrKCrossSection();

  ~G4AnalyticalEcpssrKCrossSection();
			     
  
  G4double CalculateCrossSection(G4int, G4double, G4double);//according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)
  
  G4double  ExpIntFunction(G4int n,G4double x);//Exponential Integral Function
  
private:
  
  G4AnalyticalEcpssrKCrossSection(const G4AnalyticalEcpssrKCrossSection&);
  G4AnalyticalEcpssrKCrossSection & operator = (const G4AnalyticalEcpssrKCrossSection &right);

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

  G4DNACrossSectionDataSet* tableC1;
  G4DNACrossSectionDataSet* tableC2;
  G4DNACrossSectionDataSet* tableC3;
};
  
#endif
