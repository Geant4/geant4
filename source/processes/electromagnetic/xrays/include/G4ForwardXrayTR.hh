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
// $Id: G4ForwardXrayTR.hh 108508 2018-02-15 15:54:35Z gcosmo $
//
// G4ForwardXrayTR
//
// Class for description
//
// Class for forward X-ray transition radiation generated
// by relativistic charged particle crossed interface between material 1
// and material 2 (1 -> 2)

// History:
// 22.09.97, V. Grichine (Vladimir.Grichine@cern.ch)
// 26.01.00, V.Grichine, new constructor and protected DM for fast sim. models
// 10.03.03, V.Ivanchenko migrade to "cut per region"
// 03.06.03, V.Ivanchenko fix compilation warnings

#ifndef G4FORWARDXRAYTR_H
#define G4FORWARDXRAYTR_H


#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4OpticalPhoton.hh"
#include "G4TransportationManager.hh"

#include "G4TransitionRadiation.hh"
#include "G4PhysicsTable.hh"
#include "G4Gamma.hh"
#include "G4PhysicsLogVector.hh"

class G4ForwardXrayTR : public G4TransitionRadiation
{
public:

  // Constructors

  explicit G4ForwardXrayTR(  const G4String& matName1,    //  G4Material* pMat1,
		    const G4String& matName2,    //  G4Material* pMat2,     
                    const G4String& processName="XrayTR"     );
  
  explicit G4ForwardXrayTR(  const G4String& processName="XrayTR"     );

  // Destructor //  virtual

  virtual ~G4ForwardXrayTR();

///////////////////////    Methods    /////////////////////////////////

  void BuildXrayTRtables();

  G4double GetMeanFreePath(const G4Track&, G4double,
			   G4ForceCondition* condition) override;

  G4VParticleChange* PostStepDoIt( const G4Track& aTrack,
           const G4Step&  aStep    ) override;

  G4double GetEnergyTR(G4int iMat, G4int jMat, G4int iTkin) const;
  
  G4double GetThetaTR(G4int iMat, G4int jMat, G4int iTkin) const;     


///////////////////// Angle distribution  /////////////////////////////
//

G4double SpectralAngleTRdensity( G4double energy,
                                 G4double varAngle ) const override;

G4double AngleDensity( G4double energy,
                       G4double varAngle ) const;

G4double EnergyInterval( G4double energy1,
                         G4double energy2,
                         G4double varAngle ) const;

G4double AngleSum( G4double varAngle1,
                   G4double varAngle2     )  const;

/////////////////////////  Energy distribution ///////////////////////////////

G4double SpectralDensity( G4double energy,
                          G4double x       ) const;

G4double AngleInterval( G4double energy,
                        G4double varAngle1,
                        G4double varAngle2   ) const;

G4double EnergySum( G4double energy1,
                    G4double energy2     )   const;


///////////////////////////   Access functions  ////////////////////////////

  G4PhysicsTable* GetAngleDistrTable();
  G4PhysicsTable* GetEnergyDistrTable();

  static G4int    GetSympsonNumber();
  static G4int            GetBinTR();

  static G4double GetMinProtonTkin();
  static G4double GetMaxProtonTkin();
  static G4int           GetTotBin();


protected:  // for access from X-ray TR fast simulation models

  // private :  ///////////////  Data members   ///////////////////////////

G4ParticleDefinition* fPtrGamma;  // pointer to TR photon

const std::vector<G4double>* fGammaCutInKineticEnergy;
                                     // TR photon cut in energy array
G4double  fGammaTkinCut;            // Tkin cut of TR photon in current mat.

G4PhysicsTable* fAngleDistrTable;
G4PhysicsTable* fEnergyDistrTable;

G4PhysicsLogVector* fProtonEnergyVector;

static G4int fSympsonNumber;                // Accuracy of Sympson integration 

static G4double fTheMinEnergyTR;            //  static min TR energy
static G4double fTheMaxEnergyTR;            //  static max TR energy
       G4double fMinEnergyTR;               //  min TR energy in material
       G4double fMaxEnergyTR;               //  max TR energy in material
static G4double  fTheMaxAngle;               //  max theta of TR quanta
static G4double  fTheMinAngle;               //  max theta of TR quanta
       G4double  fMaxThetaTR;               //  max theta of TR quanta
static G4int          fBinTR;               //  number of bins in TR vectors

static G4double fMinProtonTkin;             // min Tkin of proton in tables
static G4double fMaxProtonTkin;             // max Tkin of proton in tables
static G4int    fTotBin;             // number of bins in log scale
       G4double fGamma;             // current Lorentz factor

static G4double fPlasmaCof;               // physical consts for plasma energy
static G4double fCofTR;

G4double fSigma1;                       // plasma energy Sq of matter1
G4double fSigma2;                       // plasma energy Sq of matter2

private:
  // Operators

  G4ForwardXrayTR(const G4ForwardXrayTR& right) = delete;

  G4ForwardXrayTR& operator=(const G4ForwardXrayTR& right) = delete;

  // G4int operator==(const G4ForwardXrayTR& right)const;
  // G4int operator!=(const G4ForwardXrayTR& right)const;

};    // end of G4ForwardXrayTR class ---------------------------

#endif   // G4FORWARDXRAYTR_H



