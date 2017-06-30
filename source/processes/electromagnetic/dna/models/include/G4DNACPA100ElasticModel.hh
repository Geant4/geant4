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
// CPA100 elastic model class for electrons
//
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries, 
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class: 
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//

#ifndef G4DNACPA100ElasticModel_h
#define G4DNACPA100ElasticModel_h 1

#include <map>
#include "G4DNACrossSectionDataSet.hh"
#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LogLogInterpolation.hh"
//#include "G4DNACPA100LogLogInterpolation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"

class G4DNACPA100ElasticModel : public G4VEmModel
{

public:

  G4DNACPA100ElasticModel(const G4ParticleDefinition* p = 0, 
              const G4String& nam = "DNACPA100ElasticModel");

  virtual ~G4DNACPA100ElasticModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material* material,
     const G4ParticleDefinition* p,
     G4double ekin,
     G4double emin,
     G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
     const G4MaterialCutsCouple*,
     const G4DynamicParticle*,
     G4double tmin,
     G4double maxEnergy);
     
  //---
  // kept for backward compatibility
  inline void SetKillBelowThreshold (G4double threshold);   
  G4double GetKillBelowThreshold () { return killBelowEnergy; }   
  //---
  
  inline void SelectStationary(G4bool input); 

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:
  
  G4bool statCode;

  // Water density table
  const std::vector<G4double>* fpMolWaterDensity;

  G4double killBelowEnergy;  
  G4double lowEnergyLimit;  
  G4double highEnergyLimit; 
  G4bool isInitialised;
  G4int verboseLevel;
  
  // Cross section
  
  typedef std::map<G4String,G4String,std::less<G4String> > MapFile;
  MapFile tableFile;

  typedef std::map<G4String,G4DNACrossSectionDataSet*,std::less<G4String> > MapData;
  MapData tableData;
  
  // Final state

  //G4double DifferentialCrossSection(G4ParticleDefinition * aParticleDefinition, G4double k, G4double theta);

  G4double Theta(G4ParticleDefinition * aParticleDefinition, G4double k, G4double integrDiff);
  
  G4double LinLinInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);

  G4double LinLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
   
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

  TriDimensionMap eDiffCrossSectionData;
  std::vector<double> eTdummyVec;

  typedef std::map<double, std::vector<double> > VecMap;
  VecMap eVecm;
   
  G4double RandomizeCosTheta(G4double k);
   
  //
   
  G4DNACPA100ElasticModel & operator=(const  G4DNACPA100ElasticModel &right);
  G4DNACPA100ElasticModel(const  G4DNACPA100ElasticModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNACPA100ElasticModel::SelectStationary (G4bool input)
{ 
    statCode = input; 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//---
// kept for backward compatibility

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNACPA100ElasticModel::SetKillBelowThreshold (G4double /*threshold*/) 
{ 

  G4ExceptionDescription errMsg;
  errMsg << "*** WARNING : "
      << "G4DNACPA100ElasticModel::SetKillBelowThreshold"
      << "is deprecated, the kill threshold won't be taken into account";

  G4Exception ("*** WARNING : G4DNACPA100ElasticModel::SetKillBelowThreshold DEPRECATED","",JustWarning,"") ;
            
}   

#endif
