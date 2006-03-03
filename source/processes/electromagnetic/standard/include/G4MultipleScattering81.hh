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
// $Id: G4MultipleScattering81.hh,v 1.2 2006-03-03 19:40:36 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//------------- G4MultipleScattering81 physics process -------------------------
//               by Laszlo Urban, March 2001
//
// 07-08-01 new methods Store/Retrieve PhysicsTable
// 23-08-01 new angle and z distribution,energy dependence reduced,
//          Store,Retrieve methods commented out temporarily, L.Urban
// 11-09-01 G4MultipleScattering81x put as default: G4MultipleScattering81
//          Store,Retrieve methods reactived (mma)
// 13-09-01 Unused TrueToGeomTransformation method deleted,
//          class description (L.Urban)
// 19-09-01 come back to previous process name msc
// 17-04-02 NEW angle distribution + boundary algorithm modified, L.Urban
// 22-04-02 boundary algorithm modified -> important improvement in timing !!!!
//          (L.Urban)
// 24-05-02 changes in data members, L.Urban
// 30-10-02 changes in data members, L.Urban
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 05-02-03 changes in data members, L.Urban
// 28-03-03 Move to model design (V.Ivanchenko)
// 18-04-03 Change name (V.Ivanchenko)
// 16-06-03: ShortLived are not applicable any more (V.Ivanchenko)
// 17-08-04 name of data member facxsi changed to factail together
//          with the corresponding set function (L.Urban)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 15-04-05 optimize internal interfaces (V.Ivanchenko)
// 02-10-05 new algorithm for step limitation, new data members (L.Urban)
// 05-10-05 value of data member tlimitmin has been changed (L.Urban)
// 23-10-05 new Boolean data member prec (false ~ 7.1 like, true new step
//          limit in TruePathLengthLimit, L.Urban)
// 25-10-05 prec renamed to steppingAlgorithm, set function triggers
//          'default' facrange too, true - 0.02, false - 0.2 (L.Urban)
// 26-10-05 the above is put in the function MscStepLimitation() (mma)
// 05-11-05 new data member rangecut (L.Urban)
// 13-11-05 some code cleaning (L.Urban)
// 07-12-05 GeomLimit is protected instead of public
// 11-12-05 data menber rangecut removed (L.Urban)
// 19-01-07 tlimitmin = facrange*50*micrometer, i.e. it depends on the
//          value of facrange (L.Urban)
// 16-02-06 set function for data member factail (L.Urban)
//
//------------------------------------------------------------------------------
//
// $Id: G4MultipleScattering81.hh,v 1.2 2006-03-03 19:40:36 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// class description
//
//  The class simulates the multiple scattering for any kind
//  of charged particle.
//
// class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MultipleScattering81_h
#define G4MultipleScattering81_h 1

#include "G4VMultipleScattering81.hh"

class G4Navigator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MultipleScattering81 : public G4VMultipleScattering81

{
public:    // with description

  G4MultipleScattering81(const G4String& processName="msc");

 ~G4MultipleScattering81();

  // returns true for charged particles, false otherwise
  G4bool IsApplicable (const G4ParticleDefinition& p);

  G4double TruePathLengthLimit(const G4Track&,
                               G4double&,
                               G4double) {return DBL_MAX;};

  // Print few lines of informations about the process: validity range,
  void PrintInfo();

  // set boolean flag steppingAlgorithm
  // ( true/false : standard or 7.1 style process)
  void MscStepLimitation(G4bool algorithm, G4double factor = -1.);

  // geom. step length distribution should be sampled or not
  void Setsamplez(G4bool value) { samplez = value;};

  // to reduce the energy/step dependence
  void Setdtrl(G4double value) { dtrl = value;};

  // 'soften' step limitation above Tkinlimit
  void SetTkinlimit(G4double value) { tkinlimit = value;};

  // Steplimit = facrange*max(range,lambda)
  void SetFacrange(G4double val) { facrange=val;};

  // connected with step size reduction due to geometry
  void SetFacgeom(G4double val) { facgeom=val;};

  // parameter governs the tail of the angular distribution
  void SetFactail(G4double val) { factail=val;};

protected:

  // This function initialise models
  void InitialiseProcess(const G4ParticleDefinition*);

private:        // data members

  G4double lowKineticEnergy;
  G4double highKineticEnergy;
  G4int    totBins;

  G4double tkinlimit;
  G4double facrange;
  G4double facgeom;
  G4double dtrl;
  G4double factail;

  G4bool   steppingAlgorithm;
  G4bool   samplez;
  G4bool   isInitialized;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
