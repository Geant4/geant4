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
// $Id: G4MultipleScattering.hh,v 1.10 2004/12/01 19:37:13 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
//
//------------- G4MultipleScattering physics process --------------------------
//               by Laszlo Urban, March 2001
//
// 07-08-01 new methods Store/Retrieve PhysicsTable
// 23-08-01 new angle and z distribution,energy dependence reduced,
//          Store,Retrieve methods commented out temporarily, L.Urban
// 11-09-01 G4MultipleScatteringx put as default: G4MultipleScattering
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
//
//------------------------------------------------------------------------------
//
// $Id: G4MultipleScattering.hh,v 1.10 2004/12/01 19:37:13 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $

// class description
//
//  The class simulates the multiple scattering for any kind
//  of charged particle.
//
// class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MultipleScattering_h
#define G4MultipleScattering_h 1

#include "G4VMultipleScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MultipleScattering : public G4VMultipleScattering

{
public:    // with description

  G4MultipleScattering(const G4String& processName="msc");

  virtual ~G4MultipleScattering();

  G4bool IsApplicable (const G4ParticleDefinition& p)
             {return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());};
     // returns true for charged particles, false otherwise

  G4double TruePathLengthLimit(const G4Track&  track,
                                      G4double& lambda,
                                      G4double  currentMinimalStep);

  void PrintInfoDefinition();
     // Print few lines of informations about the process: validity range,

  void Setsamplez(G4bool value)               {samplez = value;};
     // geom. step length distribution should be sampled or not

  void Setdtrl(G4double value)                {dtrl = value;};
     // to reduce the energy/step dependence

  void Setfactail(G4double value)              {factail = value;}
  void SetFacrange(G4double val)              {facrange=val;
        nsmallstep = G4int(std::log((cf+facrange-1.)/facrange)/std::log(cf))+1;}
      // Steplimit after boundary crossing = facrange*range
      // estimated nb of steps at boundary nsmallstep = 1/facrange


  void SetNuclCorrPar(G4double val)            {NuclCorrPar = val;};
  void SetFactPar(G4double val)                {FactPar = val;};
     // corrs to transport cross section for high energy

protected:

  void InitialiseProcess(const G4ParticleDefinition*);
     // This function initialise models

private:

 //  hide assignment operator as  private
  G4MultipleScattering & operator = (const G4MultipleScattering &right);
  G4MultipleScattering ( const G4MultipleScattering &);

private:        // data members

  G4double lowKineticEnergy;
  G4double highKineticEnergy;
  G4int    totBins;

  G4double facrange;
  G4double tlimit;
  G4double tlimitmin;
  G4double dtrl;
  G4double NuclCorrPar;
  G4double FactPar;
  G4double factail;
  G4double cf;

  G4int    stepnolastmsc;
  G4int    nsmallstep;

  G4bool   samplez;
  G4bool   boundary;

};

#endif



