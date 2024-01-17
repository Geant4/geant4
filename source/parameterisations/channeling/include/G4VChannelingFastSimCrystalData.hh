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

#ifndef G4VChannelingFastSimCrystalData_h
#define G4VChannelingFastSimCrystalData_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VSolid.hh"
#include <unordered_map>

#include "G4ChannelingFastSimInterpolation.hh"

/** \file G4VChannelingFastSimCrystalData.hh
* \brief Definition of the G4VChannelingFastSimCrystalData class
* The class contains the data and properties related to the crystal lattice as well as
* functions to simulate of important physical processes, i.e. coulomb scattering on
* screened atomic potential, on single electrons and ionization energy losses;
* functions of electric fields, nuclear and electron densities and minimum energy
* of ionization (the corresponding interpolation coefficients are in
* G4ChannelingFastSimInterpolation).
* The functions related to the crystal geometry (transformation of coordinates and angles
* from the reference system of the bounding box of the local volume to
* the crystal lattice co-rotating reference system and vice versa) and
* initialization function SetMaterialProperties are created as virtual to make
* material data input and geometry functions flexible for modification.
*/

class G4VChannelingFastSimCrystalData{
public:

    G4VChannelingFastSimCrystalData();
    virtual ~G4VChannelingFastSimCrystalData();

    ///electric fields produced by crystal lattice
    G4double Ex(G4double x,G4double y) {return (fElectricFieldX->GetIF(x,y))*(-fZ2/fPV);}
    G4double Ey(G4double x,G4double y) {return (fElectricFieldY->GetIF(x,y))*(-fZ2/fPV);}

    ///electron density function
    G4double ElectronDensity(G4double x,G4double y)
    {
        G4double nel0=fElectronDensity->GetIF(x,y);
        if(nel0<0.) {nel0=0.;}//exception, errors of interpolation functions
        return nel0;
    }
    ///minimum energy of ionization function
    G4double MinIonizationEnergy(G4double x,G4double y)
        {return fMinIonizationEnergy->GetIF(x,y);}
    ///nuclear density function (normalized to average nuclear density)
    G4double NuclearDensity(G4double x,G4double y, G4int ielement)
        {return std::abs(fNucleiDensity[ielement]->GetIF(x,y));}
        //abs to describe exception, errors of interpolation functions,
        //don't put it =0, otherwise division on 0 in CoulombAtomicScattering

    ///Calculate the value of the Lindhard angle (!!! the value for a straight crystal)
    G4double GetLindhardAngle(G4double etotal, G4double mass);
    ///Calculate the value of the Lindhard angle (!!! the value for a straight crystal)
    G4double GetLindhardAngle();//return the Lindhard angle value calculated in
                                //SetParticleProperties

    ///Calculate simulation step (standard value for channeling particles and
    ///reduced value for overbarrier particles)
    G4double GetSimulationStep(G4double tx,G4double ty);
    ///Calculate maximal simulation step (standard value for channeling particles)
    G4double GetMaxSimulationStep(G4double etotal, G4double mass);

    ///get particle velocity/c
    G4double GetBeta(){return fBeta;}

    G4int GetNelements() {return fNelements;}
    G4int GetModel() {return iModel;}//=1 for planes, =2 for axes

    ///get bending angle of the crystal planes/axes
    ///(default BendingAngle=0 => straight crystal);
    G4double GetBendingAngle(){return fBendingAngle;}
    ///fBendingAngle MAY BE NOT THE SAME AS THE BENDING ANGLE OF THE CRYSTAL VOLUME:
    ///THE VOLUME OF A BENT CRYSTAL MAY BE G4Box, while the planes/axes inside may be bent

    G4double GetMiscutAngle(){return fMiscutAngle;}

    ///get crystal curvature
    ///for crystalline undulator the curvature is a function, otherwise it's a constant
    G4double GetCurv(G4double z){return fCU ? -fCUK2*GetCUx(z) : fCurv;}

    ///get crystalline undulator wave function
    G4double GetCUx(G4double z){return fCUAmplitude*std::cos(fCUK*z+fCUPhase);}
    ///get crystalline undulator wave 1st derivative function
    G4double GetCUtetax(G4double z){
        return fCU ? -fCUAmplitudeK*std::sin(fCUK*z+fCUPhase) : 0;}

    ///find and upload crystal lattice input files, calculate all the basic values
    ///(to do only once)
    virtual void SetMaterialProperties(const G4Material* crystal,
                                       const G4String &lattice) = 0;

    ///set geometry parameters from current logical volume
    void SetGeometryParameters(const G4LogicalVolume *crystallogic);

    ///set bending angle of the crystal planes/axes
    ///(default fBendingAngle=0 => straight crystal);
    ///only non-negative values! crystal is bent in the positive direction of x
    void SetBendingAngle(G4double tetab, const G4LogicalVolume *crystallogic);
    ///fBendingAngle MAY BE NOT THE SAME AS THE BENDING ANGLE OF THE CRYSTAL VOLUME
    ///THE VOLUME OF A BENT CRYSTAL MAY BE G4Box, while the planes/axes inside may be bent

    ///set miscut angle (default fMiscutAngle=0), acceptable range +-1 mrad,
    ///otherwise geometry routines may be unstable
    void SetMiscutAngle(G4double tetam, const G4LogicalVolume *crystallogic);

    ///set crystalline undulator parameters: amplitude, period and phase
    /// (default: all 3 value = 0)
    /// function to use in Detector Construction
    void SetCrystallineUndulatorParameters(G4double amplitude,
                                           G4double period,
                                           G4double phase,
                                           const G4LogicalVolume *crystallogic);

    ///set crystalline undulator parameters (internal function of the model)
    ///for convenience we put amplitude, period and phase in a G4ThreeVector
    void SetCUParameters(const G4ThreeVector &amplitudePeriodPhase,
                         const G4LogicalVolume *crystallogic);

    ///recalculate all the important values
    ///(to do both at the trajectory start and after energy loss)
    void SetParticleProperties(G4double etotal,
                               G4double mp,
                               G4double charge,
                               G4bool ifhadron);

    ///calculate the coordinates in the co-rotating reference system
    ///within a channel (periodic cell)
    ///(connected with crystal planes/axes either bent or straight)
    virtual G4ThreeVector CoordinatesFromBoxToLattice(const G4ThreeVector &pos0) = 0;

    ///calculate the coordinates in the Box reference system
    ///(connected with the bounding box of the volume)
    virtual G4ThreeVector CoordinatesFromLatticeToBox(const G4ThreeVector &pos) = 0;

    ///change the channel if necessary, recalculate x o y
    virtual G4ThreeVector ChannelChange(G4double& x, G4double& y, G4double& z) = 0;

    ///return correction of the longitudinal coordinate
    /// (along current plane/axis vs "central plane/axis")
    G4double GetCorrectionZ(){return fCorrectionZ;}

    ///calculate the horizontal angle in the co-rotating reference system
    ///within a channel (periodic cell)
    ///(connected with crystal planes/axes either bent or straight)
    virtual G4double AngleXFromBoxToLattice(G4double tx, G4double z)=0;

    ///calculate the horizontal angle in the Box reference system
    ///(connected with the bounding box of the volume)
    virtual G4double AngleXFromLatticeToBox(G4double tx, G4double z)=0;

    ///auxialiary function to transform the horizontal angle
    virtual G4double AngleXShift(G4double z)=0;

    ///multiple and single scattering on screened potential
    G4ThreeVector CoulombAtomicScattering(
                                 G4double effectiveStep,
                                 G4double step,
                                 G4int ielement);
    ///multiple and single scattering on electrons
    G4ThreeVector CoulombElectronScattering(G4double eMinIonization,
                                            G4double electronDensity,
                                            G4double step);
    ///ionization losses
    G4double IonizationLosses(G4double dz, G4int ielement);

  void SetVerbosity(G4int ver){fVerbosity = ver;}

protected:
    ///classes containing interpolation coefficients
    //horizontal electric field data
    G4ChannelingFastSimInterpolation* fElectricFieldX{nullptr};
    //vertical electric field data
    G4ChannelingFastSimInterpolation* fElectricFieldY{nullptr};
    //electron density data
    G4ChannelingFastSimInterpolation* fElectronDensity{nullptr};
    //minimal energy of ionization data
    G4ChannelingFastSimInterpolation* fMinIonizationEnergy{nullptr};
    //nuclear density distributions data
    std::vector <G4ChannelingFastSimInterpolation*> fNucleiDensity;

    ///values related to the crystal geometry
    G4ThreeVector fHalfDimBoundingBox;//bounding box half dimensions
    G4int fBent=0;//flag of bent crystal,
                  //=0 for straight and =1 for bent, by default straight crystal

    G4double fBendingAngle=0.;// angle of bending of the crystal planes/axes
                              //inside the crystal volume
                          //MAY BE NOT THE SAME AS THE BENDING ANGLE OF THE CRYSTAL VOLUME
                          //THE VOLUME OF A BENT CRYSTAL MAY BE G4Box,
                          //while the planes/axes inside may be bent
    G4double fBendingR = 0.; // bending radius of the crystal planes/axes
    G4double fBending2R=0.; // =2*fBendingR
    G4double fBendingRsquare=0.; // =fBendingR**2
    G4double fCurv=0.; //=1/fBendingR bending curvature of the crystal planes/axes

    G4double fMiscutAngle = 0.;// miscut angle, can be of either sign or 0;
                               //safe values |ThetaMiscut|<0.001
    G4double fCosMiscutAngle=1.;// = std::cos(fMiscutAngle), to economy operations
    G4double fSinMiscutAngle=0.;// = std::sin(fMiscutAngle), to economy operations

    G4double fCorrectionZ = 1.;//correction of the longitudinal coordinate
    //(along current plane/axis vs "central plane/axis"), 1 is default value
    //(for "central plane/axis" or a straight crystal)

    G4bool fCU = false;//flag of crystalline undulator geometry
                       //(periodically bent crystal)
    G4double fCUAmplitude=0.; //Amplitude of a crystalline undulator
    G4double fCUK=0.;    //2*pi/period of a crystalline undulator
    G4double fCUPhase=0.;//Phase of a crystalline undulator
    G4double fCUAmplitudeK=0.;//fCUAmplitude*fCUK
    G4double fCUK2=0.;   //fCUK^2

    ///values related to the crystal lattice
    G4int fNelements=1;//number of nuclear elements in a crystal
    G4int iModel=1;// model type (iModel=1 for interplanar potential,
                 //iModel=2 for the interaxial one)

    G4double fVmax=0;  // the height of the potential well
    G4double fVmax2=0; // =2*fVmax
    G4double fVMinCrystal=0;// non-zero minimal potential inside the crystal,
                         // necessary for angle recalculation for entrance/exit
                         //through the crystal lateral surface

    G4double fChangeStep=0;// fChannelingStep = fChangeStep/fTetaL

    std::vector <G4double> fI0; //Mean excitation energy

    std::vector <G4double> fRF;//Thomas-Fermi screening radius

    ///angles necessary for multiple and single coulomb scattering

    //minimal scattering angle by coulomb scattering on nuclei
    //defined by shielding by electrons
    std::vector <G4double> fTeta10;//(in the Channeling model
                                   //teta1=fTeta10/fPz*(1.13+fK40/vz**2)
    //maximal scattering angle by coulomb scattering on nuclei defined by nucleus radius
    std::vector <G4double> fTetamax0;//(in the Channeling model tetamax=fTetamax0/fPz)
    std::vector <G4double> fTetamax2;//=tetamax*tetamax
    std::vector <G4double> fTetamax12;//=teta1*teta1+tetamax*tetamax
    std::vector <G4double> fTeta12; //= teta1*teta1

    ///coefficients necessary for multiple and single coulomb scattering
    std::vector <G4double> fK20; //a useful coefficient, fK2=fK20/fPV/fPV
    std::vector <G4double> fK2;  //a useful coefficient,
                                 //fK2=(fZ2*alpha*hdc)**2*4.*pi*fN0*(fZ1/fPV)**2
    std::vector <G4double> fK40; //a useful coefficient, fK40=3.76D0*(alpha*fZ1)**2
    G4double fK30=0;//a useful coefficient, fK3=fK30/fPV/fPV
    G4double fK3=0;//a useful coefficient,  fK3=2.*pi*alpha*hdc/electron_mass_c2/(fPV)**2

    std::vector <G4double> fKD;  //a useful coefficient for dE/dx

    ///coefficients for multiple scattering suppression
    std::vector <G4double> fPu11;//a useful coefficient for exponent containing u1
    std::vector <G4double> fPzu11;//a useful coefficient for exponent containing u1
    std::vector <G4double> fBB;//a useful coefficient
    std::vector <G4double> fE1XBbb;//a useful coefficient
    std::vector <G4double> fBBDEXP;//a useful coefficient
  
  //Variable to control printout
  G4int fVerbosity = 1;


private:

    //exponential integral
    G4double expint(G4double x);    
    
    ///private variables

    std::unordered_map<G4int, G4double> fMapBendingAngle;//the map fBendingAngle
                                                           //for different logical volumes

    std::unordered_map<G4int, G4double> fMapMiscutAngle;//the map fMiscutAngle
                                                        //for different logical volumes

    std::unordered_map<G4int, G4ThreeVector> fMapCUAmplitudePeriodPhase;//the map of
                                                        //AmplitudePeriodPhase
                                                        //for different logical volumes

    G4double fChannelingStep=0;// simulation step under the channeling conditions =
                             //channeling oscillation length/fNsteps
                             // channeling oscillation length: Biryukov book Eq. (1.24)

    ///energy depended values
    G4double fPz=0; // particle momentum absolute value
    G4double fPV=0; // pv
    G4double fTetaL=0;  //Lindhard angle
    G4double fBeta=0; //particle (velocity/c)
    G4double fV2=0; // particle (velocity/c)^2
    G4double fGamma=0; //Lorentz factor
    G4double fMe2Gamma=0; // me^2*fGamma
    G4double fTmax=0; // max ionization losses

    ///particle properties flags
    G4bool fHadron=false;//=true (for hadrons); =false (for leptons)
    G4double fZ2=0; //particle charge

};

#endif
    
