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
// $Id: G4FinalStateIonisationBorn.hh,v 1.4 2008-07-14 20:47:34 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4FINALSTATEIONISATIONBORN_HH
#define G4FINALSTATEIONISATIONBORN_HH 1
 
#include "G4FinalStateProduct.hh"
#include "G4WaterIonisationStructure.hh"
#include "G4CrossSectionIonisationBornPartial.hh"

class G4Track;
class G4Step;

class G4FinalStateIonisationBorn
{
 public:
   
   G4FinalStateIonisationBorn();
   
   ~G4FinalStateIonisationBorn();
   
   const G4FinalStateProduct& GenerateFinalState(const G4Track& track, const G4Step& step);
   
private:
   
   G4double lowEnergyLimitDefault;
   G4double highEnergyLimitDefault;
   std::map<G4String,G4double,std::less<G4String> > lowEnergyLimit;
   std::map<G4String,G4double,std::less<G4String> > highEnergyLimit;

   G4FinalStateProduct product;

   G4WaterIonisationStructure waterStructure;

   G4CrossSectionIonisationBornPartial cross;
   
   G4double RandomizeEjectedElectronEnergy(G4ParticleDefinition * aParticleDefinition, G4double incomingParticleEnergy, G4int shell) ;

   void RandomizeEjectedElectronDirection(G4ParticleDefinition * aParticleDefinition, G4double incomingParticleEnergy, G4double
                                           outgoingParticleEnergy, G4double & cosTheta, G4double & phi );
   
   double  DifferentialCrossSection(G4ParticleDefinition * aParticleDefinition, G4double k, G4double energyTransfer, G4int shell);

   G4double LogLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
   
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
   TriDimensionMap eDiffCrossSectionData[6];
   TriDimensionMap pDiffCrossSectionData[6];
   std::vector<double> eTdummyVec;
   std::vector<double> pTdummyVec;

   typedef std::map<double, std::vector<double> > VecMap;
   VecMap eVecm;
   VecMap pVecm;
 
};

#endif
