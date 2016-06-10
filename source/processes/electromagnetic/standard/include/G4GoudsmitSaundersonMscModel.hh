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
// $Id: G4GoudsmitSaundersonMscModel.hh 67990 2013-03-13 10:56:28Z gcosmo $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4GoudsmitSaundersonMscModel
//
// Author:        Omrane Kadri
//
// Creation date: 20.02.2009
//
// Modifications:
// 04.03.2009 V.Ivanchenko cleanup and format according to Geant4 EM style
// 12.05.2010 O.Kadri: adding Qn1 and Qn12 as private doubles
//
// Class description:
//
// Multiple scattering model using classical Goudsmit-Saunderson model
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//REFERENCES:
//Ref.1:E. Benedito et al.,"Mixed simulation ... cross-sections", NIMB 174 (2001) pp 91-110;
//Ref.2:I. Kawrakow et al.,"On the condensed ... transport",NIMB 142 (1998) pp 253-280;
//Ref.3:I. Kawrakow et al.,"On the representation ... calculations",NIMB 134 (1998) pp 325-336;
//Ref.4:I. Kawrakow et al.,"The EGSnrc code ... Transport",NRCC Report PIRS-701, Sept. 21, 2006;
//Ref.5:F. Salvat et al.,"ELSEPA--Dirac partial ...molecules", Comp. Phys. Comm. 165 (2005) pp 157-190;
//Ref.6:G4UrbanMscModel G4_v9.1Ref09; 
//Ref.7:G4WentzelVIModel G4_v9.3.
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#ifndef G4GoudsmitSaundersonMscModel_h
#define G4GoudsmitSaundersonMscModel_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VMscModel.hh"
#include "G4PhysicsTable.hh"
#include "globals.hh"

class G4DataVector;
class G4ParticleChangeForMSC;
class G4LossTableManager;
class G4GoudsmitSaundersonTable;

class G4GoudsmitSaundersonMscModel : public G4VMscModel
{
public:

  G4GoudsmitSaundersonMscModel(const G4String& nam = "GoudsmitSaunderson");

  virtual ~G4GoudsmitSaundersonMscModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  void StartTracking(G4Track*);

  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition* particle,
					      G4double KineticEnergy,
					      G4double AtomicNumber, G4double, 
					      G4double, G4double);

  virtual G4ThreeVector& SampleScattering(const G4ThreeVector&, 
					  G4double safety);

  virtual G4double ComputeTruePathLengthLimit(const G4Track& track,
					      G4double& currentMinimalStep);

  virtual G4double ComputeGeomPathLength(G4double truePathLength);

  virtual G4double ComputeTrueStepLength(G4double geomStepLength);

private:  
  void SampleCosineTheta(G4double,G4double,G4double &,G4double &);

  void CalculateIntegrals(const G4ParticleDefinition* ,G4double,G4double,
                          G4double&,G4double&);

  void LoadELSEPAXSections();

  inline void SetParticle(const G4ParticleDefinition* p);

  inline G4double GetLambda(G4double);

  //  hide assignment operator
  G4GoudsmitSaundersonMscModel & operator=(const  G4GoudsmitSaundersonMscModel &right);
  G4GoudsmitSaundersonMscModel(const  G4GoudsmitSaundersonMscModel&);

  G4double lowKEnergy;
  G4double highKEnergy;
  G4double currentKinEnergy;
  G4double currentRange; 

  G4double smallstep,tlimitminfix,skindepth;
  G4double fr,rangeinit,masslimite,tgeom;
  G4double par1,par2,par3,zPathLength,truePathLength;
  G4double tausmall,taulim,tlimit,tlimitmin,geommin,geombig;
  G4double charge,lambdalimit;
  G4double tPathLength,stepmin ;
  G4double lambda0,lambda1,lambda11;
  G4double mass;
  G4int    currentMaterialIndex;

  G4bool   firstStep;
  G4bool   inside;
  G4bool   insideskin;

  G4GoudsmitSaundersonTable* GSTable;
  G4LossTableManager*        theManager;
  const G4ParticleDefinition* particle;
  G4ParticleChangeForMSC*     fParticleChange;
  const G4MaterialCutsCouple* currentCouple;

  static G4double ener[106];
  static G4double TCSE[103][106];
  static G4double FTCSE[103][106];
  static G4double TCSP[103][106];
  static G4double FTCSP[103][106];

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4GoudsmitSaundersonMscModel::SetParticle(const G4ParticleDefinition* p)
{
  if (p != particle) {
    particle = p;
    mass = p->GetPDGMass();
    charge = p->GetPDGCharge()/CLHEP::eplus;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

