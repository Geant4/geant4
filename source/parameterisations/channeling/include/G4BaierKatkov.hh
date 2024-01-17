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
#ifndef G4BaierKatkov_h
#define G4BaierKatkov_h 1

#include "globals.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include <vector>
#include "G4ThreeVector.hh"

#include "G4VFastSimulationModel.hh"

/** \file G4BaierKatkov.hh
* \brief Definition of the G4BaierKatkov class
* This class is designed for the calculation of radiation probability, radiation point
* and the parameters of the photon produced as well as spectrum accumulation using
* the Baier-Katkov integral:
* V. N. Baier, V. M. Katkov, and V. M. Strakhovenko,
* Electromagnetic Processes at High Energies in Oriented Single Crystals
* (World Scientific, Singapore, 1998).
*/

class G4BaierKatkov
{
public:
    // default constructor
    G4BaierKatkov();

    // destructor
    ~G4BaierKatkov() = default;

    /**
       You may call DoRadiation at each step of your trajectory
       CAUTION: please ensure that your steps are physically small enough for calculation
       of the radiation type you are interested in
       CAUTION: do ResetRadIntegral() before the start of a new trajectory

       1) change some model defaults if necessary
         (SetSinglePhotonRadiationProbabilityLimit,
          SetNSmallTrajectorySteps, SetSpectrumEnergyRange)
       2) call DoRadiation at each step of your trajectory
       3) if DoRadiation returns TRUE, this means that a photon is produced (not added
          as a secondary yet) and its parameters are calculated.
       4) You may generate a new photon using GeneratePhoton either with
          the parameters calculated in DoRadiation or your own parameters.
          CAUTION: By now GeneratePhoton works only for a FastSim model
       5) Use GetPhotonEnergyInSpectrum() and GetTotalSpectrum() to return calculated
          total spectrum (all the photons altogether)
          Caution: is not normalized on the event number
       6) Get the charged particle parameters in the radiation point:
          GetParticleNewTotalEnergy(),
          GetParticleNewAngleX(), GetParticleNewAngleY(),
          GetNewGlobalTime(),
          GetParticleNewCoordinateXYZ()
     */

    ///get functions

    /// get maximal radiation probability to preserve single photon radiation
    G4double GetSinglePhotonRadiationProbabilityLimit()
    {return fSinglePhotonRadiationProbabilityLimit;}

    ///CAUTION! : use the get functions below ONLY AFTER the call of DoRadiation
    /// and ONLY IF IT RETURNS true

    ///total probability of radiation: needs calculation of DoRadiation first
    G4double GetTotalRadiationProbability(){return fTotalRadiationProbability;}
    ///get new parameters of the particle
    ///(the parameters at the point of radiation emission)
    ///needs calculation of DoRadiation first
    G4double GetParticleNewTotalEnergy(){return fNewParticleEnergy;}
    G4double GetParticleNewAngleX(){return fNewParticleAngleX;}
    G4double GetParticleNewAngleY(){return fNewParticleAngleY;}
    G4double GetNewGlobalTime(){return fNewGlobalTime;}
    const G4ThreeVector& GetParticleNewCoordinateXYZ(){return fNewParticleCoordinateXYZ;}

  ///get photon energies (x-value in spectrum)
  const std::vector<G4double>& GetPhotonEnergyInSpectrum()
  {return fPhotonEnergyInSpectrum;}
  
  ///get fTotalSpectrum after finishing the trajectory part with DoRadiation
  const std::vector<G4double>& GetTotalSpectrum(){return fTotalSpectrum;}

    ///set functions

    ///set maximal radiation probability to preserve single photon radiation
    void SetSinglePhotonRadiationProbabilityLimit(G4double wmax)
    {fSinglePhotonRadiationProbabilityLimit = wmax;}

    ///number of steps in a trajectory small piece before
    ///the next call of the radiation integral
    void SetNSmallTrajectorySteps(G4double nSmallTrajectorySteps)
    {fNSmallTrajectorySteps = nSmallTrajectorySteps;}

    ///reinitialize intermediate integrals fFa, fSs, fSc, fSsx, fSsy, fScx, fScy;
    ///reset radiation integral internal variables to defaults;
    ///reset the trajectory and radiation probability along the trajectory
    void ResetRadIntegral();

    ///setting the number of photons in sampling of Baier-Katkov Integral
    ///(MC integration by photon energy and angles <=> photon momentum)
    void SetSamplingPhotonsNumber(G4int nPhotons){fNMCPhotons = nPhotons;}

    ///setting the number of radiation angles 1/gamma, defining the width of
    ///the angular distribution of photon sampling in the Baier-Katkov Integral
    void SetRadiationAngleFactor(G4double radiationAngleFactor)
    {fRadiationAngleFactor = radiationAngleFactor;}

    ///CAUTION, the bins width is logarithmic
    ///Do not worry if the maximal energy > particle energy.
    ///This elements of spectrum with non-physical energies
    ///will not be processed (they will be 0).
    void SetSpectrumEnergyRange(G4double emin,
                                G4double emax,
                                G4int numberOfBins);
    /// SetSpectrumEnergyRange also calls ResetRadIntegral()

    void SetMinPhotonEnergy(G4double emin){SetSpectrumEnergyRange(emin,
                                                                  fMaxPhotonEnergy,
                                                                  fNBinsSpectrum);}
    void SetMaxPhotonEnergy(G4double emax){SetSpectrumEnergyRange(fMinPhotonEnergy,
                                                                  emax,
                                                                  fNBinsSpectrum);}
    void SetNBinsSpectrum(G4int nbin){SetSpectrumEnergyRange(fMinPhotonEnergy,
                                                             fMaxPhotonEnergy,
                                                             nbin);}

    /// Increase the statistic of virtual photons in a certain energy region
    /// CAUTION! : don't do it before SetSpectrumEnergyRange or SetMinPhotonEnergy
    void AddStatisticsInPhotonEnergyRegion(G4double emin, G4double emax,
                                           G4int timesPhotonStatistics);

    /// Virtual collimator masks the selection of photon angles in fTotalSpectrum
    /// Virtual collimator doesn't influence on Geant4 simulations.
    void SetVirtualCollimator(G4double virtualCollimatorAngularDiameter)
         {fVirtualCollimatorAngularDiameter=virtualCollimatorAngularDiameter;}

    /// add the new elements of the trajectory, calculate radiation in a crystal
    /// see complete description in G4BaierKatkov::DoRadiation
    /// calls RadIntegral and all the necessary functions
    /// sets the parameters of a photon produced (if any)
    /// using SetPhotonProductionParameters()
    /// returns true in the case of photon generation, false if not
    G4bool DoRadiation(G4double etotal, G4double mass,
                       G4double angleX, G4double angleY,
                       G4double angleScatteringX, G4double angleScatteringY,
                       G4double step, G4double globalTime,
                       G4ThreeVector coordinateXYZ,
                       G4bool flagEndTrajectory=false);

    /// generates secondary photon belonging to fastStep with variables
    /// photon energy, momentum direction, coordinates and global time
    /// CALCULATED IN DoRadiation => USE IT ONLY AFTER DoRadiation returns true
    void GeneratePhoton(G4FastStep &fastStep);

private:

    ///set functions

    ///function setting the photon sampling parameters in the Baier-Katkov integral;
    ///only the maximal energy is set, while fMinPhotonEnergy is used as a minimal energy;
    ///the angles set the angular distribution (the tails are infinite)
    void SetPhotonSamplingParameters(G4double ekin,
                                   G4double minPhotonAngleX, G4double maxPhotonAngleX,
                                   G4double minPhotonAngleY, G4double maxPhotonAngleY);

    ///main functions:

    ///generation of the photons in sampling of Baier-Katkov Integral
    ///(MC integration by photon energy and angles <=> by photon momentum)
    void GeneratePhotonSampling();

    ///Baier-Katkov method: calculation of integral, spectrum, full probability;
    ///returns the total radiation probability;
    ///calculates the radiation spectrum on this trajectory piece
    G4double RadIntegral(G4double etotal, G4double mass,
                         std::vector<G4double> &vectorParticleAnglesX,
                         std::vector<G4double> &vectorParticleAnglesY,
                         std::vector<G4double> &vectorScatteringAnglesX,
                         std::vector<G4double> &vectorScatteringAnglesY,
                         std::vector<G4double> &vectorSteps,
                         G4int imin);

    ///set photon production parameters (returns false if no photon produced)
    ///accumulates fTotalSpectrum
    ///CAUTION: it is an accessory function of DoRadiation, do not use it separately
    G4bool SetPhotonProductionParameters(G4double etotal, G4double mass);

    G4int FindVectorIndex(std::vector<G4double> &myvector, G4double value);

    G4double fTotalRadiationProbability = 0.;
    G4double fSinglePhotonRadiationProbabilityLimit=0.25;//Maximal radiation
                                         //probability to preserve single photon radiation

    //number of steps in a trajectory piece before the next call of the radiation integral
    G4int fNSmallTrajectorySteps=10000;
    ///trajectory element No (the first element of the array feeded in RadIntegral)
    G4int fImin0 = 0;
    ///Monte Carlo statistics of photon sampling in Baier-Katkov with 1 trajectory
    G4int fNMCPhotons =150;
    ///the number of bins in photon spectrum
    G4int fNBinsSpectrum = 110;
    G4double fMinPhotonEnergy = 0.1*CLHEP::MeV;//min energy in spectrum output
    G4double fMaxPhotonEnergy = 1*CLHEP::GeV;  //max energy in spectrum output
    G4double fLogEmaxdEmin = 1.;// = log(fMaxPhotonEnergy/fMinPhotonEnergy),
                               // 1/normalizing coefficient in
                               // 1/E distribution between
                               // fMinPhotonEnergy and fMaxPhotonEnergy
                               // is used only for spectrum output, not for simulations
                               //(we take bremsstrahlung for photon sampling)
    G4double fLogEdEmin = 1.;   // = log(E/fMinPhotonEnergy), the same as fLogEmaxdEmin
                               // but with the particle energy as the maximal limit

    G4double fVirtualCollimatorAngularDiameter=1.;//default, infinite angle
    std::vector<G4bool> fInsideVirtualCollimator;

    ///data of the phootn energy range with additional statistics
    std::vector<G4double> fLogAddRangeEmindEmin;//=G4Log(emin/fMinPhotonEnergy)
    std::vector<G4double> fLogAddRangeEmaxdEmin;//=G4Log(emax/fMinPhotonEnergy)
    std::vector<G4int> fTimesPhotonStatistics;

    ///number of trajectories
    //(at each of the Baier-Katkov Integral is calculated for the same photons)
    G4int fItrajectories = 0;

    G4double fEph0=0;   //energy of the photon produced
    G4ThreeVector PhMomentumDirection;   //momentum direction of the photon produced

    ///Radiation integral variables
    G4double fMeanPhotonAngleX =0.;        //average angle of radiated photon direction
                                           //in sampling, x-plane
    G4double fParamPhotonAngleX=1.e-3*CLHEP::rad; //a parameter radiated photon
                                           //sampling distribution, x-plane
    G4double fMeanPhotonAngleY =0.;        //average angle of radiated photon direction
                                           //in sampling, y-plane
    G4double fParamPhotonAngleY=1.e-3*CLHEP::rad; //a parameter radiated photon
                                           //sampling distribution, y-plane
    G4double fRadiationAngleFactor = 1.; // number of radiation angles 1/gamma:
                                         // more fRadiationAngleFactor =>
                                         // higher fParamPhotonAngleX and Y

    ///new particle parameters (the parameters at the point of radiation emission)
    G4double fNewParticleEnergy=0;
    G4double fNewParticleAngleX=0;
    G4double fNewParticleAngleY=0;
    G4double fNewGlobalTime=0;
    G4ThreeVector fNewParticleCoordinateXYZ;

    ///sampling of the energy and the angles of a photon emission
    ///(integration variables, Monte Carlo integration)
    std::vector<G4double> fPhotonEnergyInIntegral;
    std::vector<G4double> fPhotonAngleInIntegralX;
    std::vector<G4double> fPhotonAngleInIntegralY;
    std::vector<G4double> fPhotonAngleNormCoef;
    ///spectrum bin index for each photon
    std::vector<G4double> fIBinsSpectrum;
    ///the vector of the discrete CDF of the radiation of sampling photons
    std::vector<G4double> fPhotonProductionCDF;

    ///vectors of the trajectory
    std::vector<G4double> fParticleAnglesX;
    std::vector<G4double> fParticleAnglesY;
    std::vector<G4double> fScatteringAnglesX;
    std::vector<G4double> fScatteringAnglesY;
    std::vector<G4double> fSteps;
    std::vector<G4double> fGlobalTimes;
    std::vector<G4ThreeVector> fParticleCoordinatesXYZ;

    ///intermediate integrals (different for each photon energy value)!!!
    std::vector<G4double> fFa;//phase
    std::vector<G4double> fSs;
    std::vector<G4double> fSc;
    std::vector<G4double> fSsx;
    std::vector<G4double> fSsy;
    std::vector<G4double> fScx;
    std::vector<G4double> fScy;

    ///output
    std::vector<G4double> fPhotonEnergyInSpectrum; //energy values in spectrum

    std::vector<G4int> fNPhotonsPerBin; //number of photons per spectrum bin
                                   //(accumulating during total run)

    std::vector<G4double> fSpectrum; //spectrum normalized by the total
                        //radiation probability of one particle at one call of RadIntegral

    std::vector<std::vector<G4double>> fAccumSpectrum; //accumulate Spectrum during
                                             //the part of a trajectory

    std::vector<G4double> fAccumTotalSpectrum; //spectrum normalized by the total
                                          //radiation probability summed
                                          //for all the particles (is not divided
                                          //of one particle number fNPhotonsPerBin)

    std::vector<G4double> fTotalSpectrum; //spectrum normalized by
                                     //the total radiation probability summed
                                     //for all the particles
                                     //(is divided by the photon number fNPhotonsPerBin)
                                     //multiplied by the number of trajectories
                                     //(fItrajectories)

    std::vector<G4double> fImax0; //trajectory element numbers at the end of each
                            //small piece; G4double just for security of some operations
    ///total radiation probability along this trajectory
    std::vector<G4double> fTotalRadiationProbabilityAlongTrajectory;
};

#endif
