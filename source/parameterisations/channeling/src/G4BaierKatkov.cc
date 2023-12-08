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

#include "G4BaierKatkov.hh"

#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Log.hh"

G4BaierKatkov::G4BaierKatkov()
{
    //sets the default spectrum energy range of integration and
    //calls ResetRadIntegral()
    SetSpectrumEnergyRange(0.1*MeV,1.*GeV,110);

    //Do not worry if the maximal energy > particle energy
    //this elements of spectrum with non-physical energies
    //will not be processed (they will be 0)

    G4cout << " "<< G4endl;
    G4cout << "G4BaierKatkov model is activated."<< G4endl;
    G4cout << " "<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BaierKatkov::ResetRadIntegral()
{
    fAccumSpectrum.clear();

    //reinitialize intermediate integrals with zeros
    fFa.resize(fNMCPhotons);
    fSs.resize(fNMCPhotons);
    fSc.resize(fNMCPhotons);
    fSsx.resize(fNMCPhotons);
    fSsy.resize(fNMCPhotons);
    fScx.resize(fNMCPhotons);
    fScy.resize(fNMCPhotons);
    std::fill(fFa.begin(), fFa.end(), 0.);
    std::fill(fSs.begin(), fSs.end(), 0.);
    std::fill(fSc.begin(), fSc.end(), 0.);
    std::fill(fSsx.begin(), fSsx.end(), 0.);
    std::fill(fSsy.begin(), fSsy.end(), 0.);
    std::fill(fScx.begin(), fScx.end(), 0.);
    std::fill(fScy.begin(), fScy.end(), 0.);

    //Reset radiation integral internal variables to defaults
    fMeanPhotonAngleX =0.;        //average angle of
                                  //radiated photon direction in sampling, x-plane
    fParamPhotonAngleX=1.e-3*rad; //a parameter of
                                  //radiated photon sampling distribution, x-plane
    fMeanPhotonAngleY =0.;        //average angle of
                                  //radiated photon direction in sampling, y-plane
    fParamPhotonAngleY=1.e-3*rad; //a parameter of
                                  //radiated photon sampling distribution, y-plane

    fImin0 = 0;//set the first vector element to 0

    //reset the trajectory
    fParticleAnglesX.clear();
    fParticleAnglesY.clear();
    fScatteringAnglesX.clear();
    fScatteringAnglesY.clear();
    fSteps.clear();
    fGlobalTimes.clear();
    fParticleCoordinatesXYZ.clear();

    //resets the vector of element numbers at the trajectory start
    fImax0.clear();
    //sets 0 element of the vector of element numbers
    fImax0.push_back(0.);

    //resets the radiation probability
    fTotalRadiationProbabilityAlongTrajectory.clear();
    //sets the radiation probability at the trajectory start
    fTotalRadiationProbabilityAlongTrajectory.push_back(0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BaierKatkov::SetSpectrumEnergyRange(G4double emin,
                                           G4double emax,
                                           G4int numberOfBins)
{
    fMinPhotonEnergy = emin;
    fMaxPhotonEnergy = emax;
    fNBinsSpectrum = numberOfBins;

    fLogEmaxdEmin = G4Log(fMaxPhotonEnergy/fMinPhotonEnergy);

    //in initializing fNPhotonsPerBin
    fNPhotonsPerBin.resize(fNBinsSpectrum);
    std::fill(fNPhotonsPerBin.begin(), fNPhotonsPerBin.end(), 0);

    //initializing the Spectrum
    fSpectrum.resize(fNBinsSpectrum);
    std::fill(fSpectrum.begin(), fSpectrum.end(), 0);

    //initializing the fAccumTotalSpectrum
    fAccumTotalSpectrum.resize(fNBinsSpectrum);
    std::fill(fAccumTotalSpectrum.begin(), fAccumTotalSpectrum.end(), 0);

    //initializing the fTotalSpectrum
    fTotalSpectrum.resize(fNBinsSpectrum);
    std::fill(fTotalSpectrum.begin(), fTotalSpectrum.end(), 0);

    fPhotonEnergyInSpectrum.clear();
    for (G4int j=0;j<fNBinsSpectrum;j++)
    {
        //bin position (mean between 2 bin limits)
        fPhotonEnergyInSpectrum.push_back(fMinPhotonEnergy*
                                       (std::exp(fLogEmaxdEmin*j/fNBinsSpectrum)+
                                        std::exp(fLogEmaxdEmin*(j+1)/fNBinsSpectrum))/2.);
    }

    fItrajectories = 0;

    ResetRadIntegral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BaierKatkov::AddStatisticsInPhotonEnergyRegion(G4double emin,
                                                      G4double emax,
                                                      G4int timesPhotonStatistics)
{

    if(timesPhotonStatistics<=1)
    {
        G4cout << "G4BaierKatkov model, "
                  "function AddStatisticsInPhotonEnergyRegion("
                              << emin/CLHEP::MeV << " MeV, "
                              << emax/CLHEP::MeV << " MeV, "
                              << timesPhotonStatistics << ")" << G4endl;
        G4cout << "Warning: the statistics factor cannot be <=1." << G4endl;
        G4cout << "The statistics was not added." << G4endl;
        G4cout << " "<< G4endl;
    }
    else if(fMinPhotonEnergy>emin)
    {
        G4cout << "G4BaierKatkov model, "
                  "function AddStatisticsInPhotonEnergyRegion("
                              << emin/CLHEP::MeV << " MeV, "
                              << emax/CLHEP::MeV << " MeV, "
                              << timesPhotonStatistics << ")" << G4endl;
        G4cout << "Warning: the minimal energy inserted is less then "
                  "the minimal energy cut of the spectrum: "
               << fMinPhotonEnergy/CLHEP::MeV << " MeV." << G4endl;
        G4cout << "The statistics was not added." << G4endl;
        G4cout << " "<< G4endl;
    }
    else if(emax-emin<DBL_EPSILON)
    {
        G4cout << "G4BaierKatkov model, "
                  "function AddStatisticsInPhotonEnergyRegion("
                              << emin/CLHEP::MeV << " MeV, "
                              << emax/CLHEP::MeV << " MeV, "
                              << timesPhotonStatistics << ")" << G4endl;
        G4cout << "Warning: the maximal energy <= the minimal energy." << G4endl;
        G4cout << "The statistics was not added." << G4endl;
        G4cout << " "<< G4endl;
    }
    else
    {
        G4bool setrange = true;
        G4double logAddRangeEmindEmin = G4Log(emin/fMinPhotonEnergy);
        G4double logAddRangeEmaxdEmin = G4Log(emax/fMinPhotonEnergy);

        G4int nAddRange = (G4int)fTimesPhotonStatistics.size();
        for (G4int j=0;j<nAddRange;j++)
        {
            if((logAddRangeEmindEmin>=fLogAddRangeEmindEmin[j]&&
                logAddRangeEmindEmin< fLogAddRangeEmaxdEmin[j])||
               (logAddRangeEmaxdEmin> fLogAddRangeEmindEmin[j]&&
                logAddRangeEmaxdEmin<=fLogAddRangeEmaxdEmin[j])||
               (logAddRangeEmindEmin<=fLogAddRangeEmindEmin[j]&&
                logAddRangeEmaxdEmin>=fLogAddRangeEmaxdEmin[j]))
            {
                G4cout << "G4BaierKatkov model, "
                          "function AddStatisticsInPhotonEnergyRegion("
                                      << emin/CLHEP::MeV << " MeV, "
                                      << emax/CLHEP::MeV << " MeV, "
                                      << timesPhotonStatistics << ")" << G4endl;
               G4cout << "Warning: the energy range intersects another "
                         "added energy range." << G4endl;
               G4cout << "The statistics was not added." << G4endl;
               G4cout << " "<< G4endl;
               setrange = false;
               break;
            }
        }
        if (setrange)
        {
            fLogAddRangeEmindEmin.push_back(logAddRangeEmindEmin);
            fLogAddRangeEmaxdEmin.push_back(logAddRangeEmaxdEmin);
            fTimesPhotonStatistics.push_back(timesPhotonStatistics);

            G4cout << "G4BaierKatkov model: increasing the statistics of photon sampling "
                      "in Baier-Katkov with a factor of "
                   << timesPhotonStatistics << G4endl;
            G4cout << "in the energy spectrum range: ("
                   << emin/CLHEP::MeV << " MeV, "
                   << emax/CLHEP::MeV << " MeV)" << G4endl;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BaierKatkov::SetPhotonSamplingParameters(G4double ekin,
                                                G4double minPhotonAngleX,
                                                G4double maxPhotonAngleX,
                                                G4double minPhotonAngleY,
                                                G4double maxPhotonAngleY)
{
    fLogEdEmin = G4Log(ekin/fMinPhotonEnergy);
    fMeanPhotonAngleX  = (maxPhotonAngleX+minPhotonAngleX)/2.;
    fParamPhotonAngleX = (maxPhotonAngleX-minPhotonAngleX)/2.;
    fMeanPhotonAngleY  = (maxPhotonAngleY+minPhotonAngleY)/2.;
    fParamPhotonAngleY = (maxPhotonAngleY-minPhotonAngleY)/2.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BaierKatkov::GeneratePhotonSampling()
{
    fPhotonEnergyInIntegral.clear();
    fPhotonAngleInIntegralX.clear();
    fPhotonAngleInIntegralY.clear();
    fPhotonAngleNormCoef.clear();
    fInsideVirtualCollimator.clear();
    fIBinsSpectrum.clear();

    G4double ksi=0.;
    std::vector<G4int> moreStatistics;
    moreStatistics.resize(fTimesPhotonStatistics.size());
    std::fill(moreStatistics.begin(), moreStatistics.end(), 0);
    G4int nAddRange = (G4int)fTimesPhotonStatistics.size();

    //sampling of the energy of a photon emission
    //(integration variables, Monte Carlo integration)
    for (G4int j=0;j<fNMCPhotons;j++)
    {
        ksi = G4UniformRand()*fLogEdEmin;
        fIBinsSpectrum.push_back((G4int)std::trunc(
                                     ksi*fNBinsSpectrum/fLogEmaxdEmin));
        //we consider also the energy outside the spectrum  output range
        //(E>Emax => fLogEdEmin>fLogEmaxdEmin)
        //in this case we don't count the photon in the spectrum output
        if(fIBinsSpectrum[j]<fNBinsSpectrum) {fNPhotonsPerBin[fIBinsSpectrum[j]]+=1;}

        fPhotonEnergyInIntegral.push_back(fMinPhotonEnergy*std::exp(ksi));

        fPhotonAngleNormCoef.push_back(1.);

        for (G4int j2=0;j2<nAddRange;j2++)
        {
            if(ksi>fLogAddRangeEmindEmin[j2]&&
               ksi<fLogAddRangeEmaxdEmin[j2])
            {
               //calculating the current statistics in this region
               //to increase it proportionally
               moreStatistics[j2]+=1;
               fPhotonAngleNormCoef[j]/=fTimesPhotonStatistics[j2];
               break;
            }
        }
    }

    for (G4int j2=0;j2<nAddRange;j2++)
    {
        G4int totalAddRangeStatistics = moreStatistics[j2]*fTimesPhotonStatistics[j2];
        for (G4int j=moreStatistics[j2];j<totalAddRangeStatistics;j++)
        {
            ksi = fLogAddRangeEmindEmin[j2]+
                  G4UniformRand()*(std::min(fLogAddRangeEmaxdEmin[j2],fLogEdEmin)-
                                   fLogAddRangeEmindEmin[j2]);
            fIBinsSpectrum.push_back((G4int)std::trunc(
                                         ksi*fNBinsSpectrum/fLogEmaxdEmin));
         /*   //we consider also the energy outside the spectrum  output range
            //(E>Emax => fLogEdEmin>fLogEmaxdEmin)
            //in this case we don't count the photon in the spectrum output
            if(fIBinsSpectrum.back()<fNBinsSpectrum)
                {fNPhotonsPerBin[fIBinsSpectrum.back()]+=1;}*/

            fPhotonEnergyInIntegral.push_back(fMinPhotonEnergy*std::exp(ksi));

            fPhotonAngleNormCoef.push_back(1./fTimesPhotonStatistics[j2]);
        }
    }

    G4double rho=1.;
    const G4double rhocut=15.;//radial angular cut of the distribution
    G4double norm=std::atan(rhocut*rhocut)*
                            CLHEP::pi*fParamPhotonAngleX*fParamPhotonAngleY;

    //sampling of the angles of a photon emission
    //(integration variables, Monte Carlo integration)
    G4int nmctotal = (G4int)fPhotonEnergyInIntegral.size();
    for (G4int j=0;j<nmctotal;j++)
    {
        //photon distribution with long tails (useful to not exclude particle angles
        //after a strong single scattering)
        //at ellipsescale < 1 => half of statistics of photons
        do
        {
            rho = std::sqrt(std::tan(CLHEP::halfpi*G4UniformRand()));
        }
        while (rho>rhocut);

        ksi = G4UniformRand();
        fPhotonAngleInIntegralX.push_back(fMeanPhotonAngleX+
                                         fParamPhotonAngleX*
                                          rho*std::cos(CLHEP::twopi*ksi));
        fPhotonAngleInIntegralY.push_back(fMeanPhotonAngleY+
                                         fParamPhotonAngleY*
                                          rho*std::sin(CLHEP::twopi*ksi));
        fPhotonAngleNormCoef[j]*=(1.+rho*rho*rho*rho)*norm;

        //test if the photon with these angles enter the virtual collimator
        //(doesn't influence the Geant4 simulations,
        //but only the accumulation of fTotalSpectrum
        fInsideVirtualCollimator.push_back(fVirtualCollimatorAngularDiameter >
                                           std::sqrt(fPhotonAngleInIntegralX[j]*
                                                     fPhotonAngleInIntegralX[j]+
                                                     fPhotonAngleInIntegralY[j]*
                                                     fPhotonAngleInIntegralY[j]));
    }
    //reinitialize the vector of radiation CDF for each photon
    fPhotonProductionCDF.resize(nmctotal+1);//0 element equal to 0
    std::fill(fPhotonProductionCDF.begin(), fPhotonProductionCDF.end(), 0.);

    //if we have additional photons
    if (nmctotal>fNMCPhotons)
    {
        //reinitialize intermediate integrals with zeros again
        fFa.resize(nmctotal);
        fSs.resize(nmctotal);
        fSc.resize(nmctotal);
        fSsx.resize(nmctotal);
        fSsy.resize(nmctotal);
        fScx.resize(nmctotal);
        fScy.resize(nmctotal);
        std::fill(fFa.begin(), fFa.end(), 0.);
        std::fill(fSs.begin(), fSs.end(), 0.);
        std::fill(fSc.begin(), fSc.end(), 0.);
        std::fill(fSsx.begin(), fSsx.end(), 0.);
        std::fill(fSsy.begin(), fSsy.end(), 0.);
        std::fill(fScx.begin(), fScx.end(), 0.);
        std::fill(fScy.begin(), fScy.end(), 0.);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BaierKatkov::RadIntegral(G4double etotal, G4double mass,
                                    std::vector<G4double> &vectorParticleAnglesX,
                                    std::vector<G4double> &vectorParticleAnglesY,
                                    std::vector<G4double> &vectorScatteringAnglesX,
                                    std::vector<G4double> &vectorScatteringAnglesY,
                                    std::vector<G4double> &vectorSteps,
                                    G4int imin)
{
    //preliminary values are defined:

    G4double om=0.;// photon energy
    G4double eprime=0., eprime2=0.; //E'=E-omega eprime2=eprime*eprime
    G4double omprime=0.,omprimed2=0.;//om'=(E*om/E'), omprimed2=omprime/2
    G4double vxin =0.,vyin=0.,vxno=0.,vyno=0.;

    G4double dzmod=0.;
    G4double fa1=0.,faseBefore=0.,faseBeforedz=0.,faseBeforedzd2=0.;
    G4double faseAfter=0.,fa2dfaseBefore2=0.;

    G4double skJ=0, skIx=0., skIy=0.;
    G4double sinfa1=0.,cosfa1=0.;
    G4double i2=0.,j2=0.;// Ix^2+Iy^2 and of BK Jvector^2

    std::size_t nparts=vectorParticleAnglesX.size();
    G4int kmin = imin;
    if(imin==0) {kmin=1;}//skipping 0 trajectory element

    //total radiation probability for each photon
    G4double totalRadiationProbabilityPhj = 0.;

    //total radiation probability along this trajectory (fill with 0 only new elements)
    fTotalRadiationProbabilityAlongTrajectory.resize(nparts);

    //reset Spectrum
    std::fill(fSpectrum.begin(), fSpectrum.end(), 0.);

    //intermediate vectors to reduce calculations
    std::vector<G4double> axt;//acceleration of a charged particle in a horizontal plane
    axt.resize(nparts);
    std::vector<G4double> ayt;//acceleration of a charged particle in a vertical plane
    ayt.resize(nparts);
    std::vector<G4double> dz;//step in in MeV^-1
    dz.resize(nparts);
    //setting values interesting for us
    for(std::size_t k=kmin;k<nparts;k++)
    {
        dz[k]=vectorSteps[k]/CLHEP::hbarc;// dz in MeV^-1

        // accelerations
        axt[k]=(vectorParticleAnglesX[k]-vectorScatteringAnglesX[k]-
                vectorParticleAnglesX[k-1])/dz[k];
        ayt[k]=(vectorParticleAnglesY[k]-vectorScatteringAnglesY[k]-
                vectorParticleAnglesY[k-1])/dz[k];
    }

    //intermediate variables to reduce calculations:
    //the calculations inside for (G4int j=0;j<fNMCPhotons;j++)
    //{for(G4int k=kmin;k<nparts;k++){...
    //are the main cpu time consumption
    G4double e2 = etotal*etotal;
    G4double gammaInverse2 = mass*mass/(etotal*etotal);// 1/gamma^2 of
                                                       //the radiating charge particle
    G4double coefNormLogdNMC = fLogEdEmin/fNMCPhotons;//here fNMCPhotons is correct,
                                          //additional photons have been already
                                          //taken into account in weights
    G4double coefNorm = CLHEP::fine_structure_const/(8*(CLHEP::pi2))*coefNormLogdNMC;
    G4double e2pluseprime2 = 0.;//e2pluseprime2 =e2+eprime2
    G4double coefNormom2deprime2 = 0.; //coefNormom2deprime2 = coefNorm*om2/eprime2;
    G4double gammaInverse2om = 0.; //gammaInverse2*om

    std::size_t nmctotal = fPhotonEnergyInIntegral.size();
    for (std::size_t j=0;j<nmctotal;j++)
        //the final number of photons may be different from fNMCPhotons
    {
        om = fPhotonEnergyInIntegral[j];
        eprime=etotal-om; //E'=E-omega
        eprime2 = eprime*eprime;
        e2pluseprime2 =e2+eprime2;
        omprime=etotal*om/eprime;//om'=(E*om/E')
        omprimed2=omprime/2;
        coefNormom2deprime2 = coefNorm*om*om/eprime2;
        gammaInverse2om = gammaInverse2*om;

        for(std::size_t k=kmin;k<nparts;k++)
        {
            //the angles of a photon produced (with incoherent scattering)
            vxin = vectorParticleAnglesX[k]-fPhotonAngleInIntegralX[j];
            vyin = vectorParticleAnglesY[k]-fPhotonAngleInIntegralY[j];
            //the angles of a photon produced (without incoherent scattering)
            vxno = vxin-vectorScatteringAnglesX[k];
            vyno = vyin-vectorScatteringAnglesY[k];

            //phase difference before scattering
            faseBefore=omprimed2*(gammaInverse2+vxno*vxno+vyno*vyno);//phi' t<ti//MeV

            faseBeforedz = faseBefore*dz[k];
            faseBeforedzd2 = faseBeforedz/2.;
            fFa[j]+=faseBeforedz; //
            fa1=fFa[j]-faseBeforedzd2;//
            dzmod=2*std::sin(faseBeforedzd2)/faseBefore;//MeV^-1

            //phi''/faseBefore^2
            fa2dfaseBefore2 = omprime*(axt[k]*vxno+ayt[k]*vyno)/(faseBefore*faseBefore);

            //phase difference after scattering
            faseAfter=omprimed2*(gammaInverse2+vxin*vxin+vyin*vyin);//phi' ti+O//MeV

            skJ=1/faseAfter-1/faseBefore-fa2dfaseBefore2*dzmod;//MeV^-1
            skIx=vxin/faseAfter-vxno/faseBefore+dzmod*(axt[k]/faseBefore-
                                                       vxno*fa2dfaseBefore2);
            skIy=vyin/faseAfter-vyno/faseBefore+dzmod*(ayt[k]/faseBefore-
                                                       vyno*fa2dfaseBefore2);

            sinfa1 = std::sin(fa1);
            cosfa1 = std::cos(fa1);

            fSs[j]+=sinfa1*skJ;//sum sin integral J of BK
            fSc[j]+=cosfa1*skJ;//sum cos integral J of BK
            fSsx[j]+=sinfa1*skIx;// sum sin integral Ix of BK
            fSsy[j]+=sinfa1*skIy;// sum sin integral Iy of BK
            fScx[j]+=cosfa1*skIx;// sum cos integral Ix of BK
            fScy[j]+=cosfa1*skIy;// sum cos integral Iy of BK

            i2=fSsx[j]*fSsx[j]+fScx[j]*fScx[j]+fSsy[j]*fSsy[j]+fScy[j]*fScy[j];//MeV^-2
            j2=fSs[j]*fSs[j]+fSc[j]*fSc[j];//MeV^-2

            //updating the total radiation probability along the trajectory
            totalRadiationProbabilityPhj = coefNormom2deprime2*fPhotonAngleNormCoef[j]*
                    (i2*e2pluseprime2+j2*gammaInverse2om);
            fTotalRadiationProbabilityAlongTrajectory[k] += totalRadiationProbabilityPhj;
        }

        fPhotonProductionCDF[j+1] = fTotalRadiationProbabilityAlongTrajectory.back();

        //filling spectrum (adding photon probabilities to a histogram)
        //we consider also the energy outside the spectrum  output range
        //(E>Emax => fLogEdEmin>fLogEmaxdEmin)
        //in this case we don't count the photon in the spectrum output
        //we fill the spectrum only in case of the angles inside the virtual collimator
        if((fIBinsSpectrum[j]<fNBinsSpectrum)&&fInsideVirtualCollimator[j])
            {fSpectrum[fIBinsSpectrum[j]] += totalRadiationProbabilityPhj/
                    (om*coefNormLogdNMC);}

    } // end cycle

    fAccumSpectrum.push_back(fSpectrum);

    return fTotalRadiationProbabilityAlongTrajectory.back();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4BaierKatkov::FindVectorIndex(std::vector<G4double> &myvector, G4double value)
{
    auto iteratorbegin = myvector.begin();
    auto iteratorend   = myvector.end();

    //vector index (for non precise values lower_bound gives upper value)
    auto loweriterator = std::lower_bound(iteratorbegin, iteratorend, value);
    //return the index of the vector element
    return (G4int)std::distance(iteratorbegin, loweriterator);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4BaierKatkov::SetPhotonProductionParameters(G4double etotal, G4double mass)
{
    //flag if a photon was produced
    G4bool flagPhotonProduced = false;

    //how many small pieces of a trajectory are
    //before radiation (all if not radiation)
    std::size_t nodeNumber = fAccumSpectrum.size()-1;

    //ksi is a random number
    //Generally ksi = G4UniformRand() is ok, but
    //we use as a correction for the case
    //when the radiation probability becomes too high (> 0.1)
    G4double ksi = -G4Log(G4UniformRand());

    if (ksi< fTotalRadiationProbabilityAlongTrajectory.back())  // photon produced
    {
      G4double ksi1 = G4UniformRand()*fTotalRadiationProbabilityAlongTrajectory.back();

      //randomly choosing the photon to be produced from the sampling list
      //according to the probabilities calculated in the Baier-Katkov integral
      G4int iphoton = FindVectorIndex(fPhotonProductionCDF,ksi1)-1;//index of
                                                                  //a photon produced

      //energy of the photon produced
      fEph0 = fPhotonEnergyInIntegral[iphoton];
      //anagles of the photon produced
      G4double photonAngleX = fPhotonAngleInIntegralX[iphoton];
      G4double photonAngleY = fPhotonAngleInIntegralY[iphoton];

      G4double momentumDirectionZ = 1./
          std::sqrt(1.+std::pow(std::tan(photonAngleX),2)+
                       std::pow(std::tan(photonAngleY),2));

      //momentum direction vector of the photon produced
      PhMomentumDirection.set(momentumDirectionZ*std::tan(photonAngleX),
                              momentumDirectionZ*std::tan(photonAngleY),
                              momentumDirectionZ);

      //random calculation of the radiation point index (iNode)
      //ksi = G4UniformRand()*fTotalRadiationProbabilityAlongTrajectory.back();

      //sort fTotalRadiationProbabilityAlongTrajectory
      //(increasing but oscillating function => non-monotonic)
      std::vector<G4double> temporaryVector;
      temporaryVector.assign(fTotalRadiationProbabilityAlongTrajectory.begin(),
                             fTotalRadiationProbabilityAlongTrajectory.end());
      std::sort(temporaryVector.begin(), temporaryVector.end());

      //index of the point of radiation ("poststep") in sorted vector
      G4int iNode = FindVectorIndex(temporaryVector,ksi);

      //index of the point of radiation ("poststep") in unsorted vector
      auto it = std::find_if(fTotalRadiationProbabilityAlongTrajectory.begin(),
                             fTotalRadiationProbabilityAlongTrajectory.end(),
                             [&](G4double value)
                       {return std::abs(value - temporaryVector[iNode]) < DBL_EPSILON;});
      iNode = (G4int)std::distance(fTotalRadiationProbabilityAlongTrajectory.begin(), it);

      //the piece of trajectory number (necessary only for the spectrum output)
      nodeNumber = FindVectorIndex(fImax0,iNode*1.)-1;

      //set new parameters of the charged particle
      //(returning the particle back to the radiation point, i.e.
      //remembering the new parameters for corresponding get functions)
      fNewParticleEnergy = etotal-fEph0;
      fNewParticleAngleX = fParticleAnglesX[iNode];
      fNewParticleAngleY = fParticleAnglesY[iNode];
      fNewGlobalTime = fGlobalTimes[iNode];
      fNewParticleCoordinateXYZ = fParticleCoordinatesXYZ[iNode];

      //particle angle correction from momentum conservation
      //(important for fEph0 comparable to E,
      // may kick off a particle from channeling)
      G4double pratio = fEph0/std::sqrt(etotal*etotal-mass*mass);
      fNewParticleAngleX -= std::asin(pratio*std::sin(photonAngleX));
      fNewParticleAngleY -= std::asin(pratio*std::sin(photonAngleY));

      flagPhotonProduced = true;
    }

    //accumulation during entire code run
    for (G4int j=0;j<fNBinsSpectrum;j++)
    {
        fAccumTotalSpectrum[j] += fAccumSpectrum[nodeNumber][j];
        //in the case of missing photon in spectrum, probability is 0
        //(may happen only at the beginning of simulations)
        if(fNPhotonsPerBin[j]==0)
        {
           fTotalSpectrum[j] = 0;
        }
        else
        {
           fTotalSpectrum[j] = fAccumTotalSpectrum[j]/fNPhotonsPerBin[j]*fItrajectories;
        }
    }

    return flagPhotonProduced;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4BaierKatkov::DoRadiation(G4double etotal, G4double mass,
                                  G4double angleX, G4double angleY,
                                  G4double angleScatteringX, G4double angleScatteringY,
                                  G4double step, G4double globalTime,
                                  G4ThreeVector coordinateXYZ,
                                  G4bool flagEndTrajectory)
{
    /**
    DoRadiation description

    The real trajectory is divided onto parts.
    Every part is accumulated until the radiation cannot be considered as single photon,
    i.e. the radiation probability exceeds some threshold
    fSinglePhotonRadiationProbabilityLimit.
    Typically fSinglePhotonRadiationProbabilityLimit should be between 0.01 and 0.1.
    Then the photon generation as well as its parameters are simulated
    in SetPhotonProductionParameters()

    Finally if radiation took place, this function sets the new particle parameters
    (the parameters at the point of radiation emission)
    fNewParticleEnergy; fNewParticleAngleX;
    fNewParticleAngleY; NewStep; fNewGlobalTime; fNewParticleCoordinateXYZ;

    In order to control if the trajectory part reaches this limit,
    one needs to divide it into smaller pieces (consisting of
    fNSmallTrajectorySteps elements, typically several hundred) and to calculate the total
    radiation probability along the trajectory part after each small piece accomplished.
    If it exceeds the fSinglePhotonRadiationProbabilityLimit,
    the trajectory part is over and the radiation can be generated.

    In order to calculate the radiation probability the Baier-Katkov integral is called
    after each small piece of the trajectory. It is summarized with the integral along
    the previous small pieces for this part of the trajectory.

    To speed-up the simulations, the photon angles are generated only once
    at the start of the trajectory part.
    */

    //flag if a photon was produced
    G4bool flagPhotonProduced = false;

    //adding the next trajectory element to the vectors
    fParticleAnglesX.push_back(angleX);
    fParticleAnglesY.push_back(angleY);
    fScatteringAnglesX.push_back(angleScatteringX);
    fScatteringAnglesY.push_back(angleScatteringY);
    fSteps.push_back(step);
    fGlobalTimes.push_back(globalTime);
    fParticleCoordinatesXYZ.push_back(coordinateXYZ);

    G4double imax = fSteps.size();
    if((imax==fImin0+fNSmallTrajectorySteps)||flagEndTrajectory)
    {
        //set the angular limits at the start of the trajectory part
        if(fImin0==0)
        {
            //radiation within the angle = +-fRadiationAngleFactor/gamma
            G4double radiationAngleLimit=fRadiationAngleFactor*mass/etotal;

            SetPhotonSamplingParameters(etotal-mass,
            *std::min_element(fParticleAnglesX.begin(),
                              fParticleAnglesX.end())-radiationAngleLimit,
            *std::max_element(fParticleAnglesX.begin(),
                              fParticleAnglesX.end())+radiationAngleLimit,
            *std::min_element(fParticleAnglesY.begin(),
                              fParticleAnglesY.end())-radiationAngleLimit,
            *std::max_element(fParticleAnglesY.begin(),
                              fParticleAnglesY.end())+radiationAngleLimit);

            //calculation of angles of photon emission
            //(these angles are integration variables, Monte Carlo integration)
            GeneratePhotonSampling();
        }

        //calculate Baier-Katkov integral after this
        //small piece of trajectory (fNSmallTrajectorySteps elements)
        fTotalRadiationProbability = RadIntegral(etotal,mass,
                                                fParticleAnglesX,fParticleAnglesY,
                                                fScatteringAnglesX,fScatteringAnglesY,
                                                fSteps, fImin0);

        //setting the last element of this small trajectory piece
        fImin0 = imax;
        fImax0.push_back(imax*1.);

        //if the radiation probability is high enough or if we are finishing
        //our trajectory => to simulate the photon emission
        if(fTotalRadiationProbability>fSinglePhotonRadiationProbabilityLimit||
                flagEndTrajectory)
        {
            fItrajectories += 1; //count this trajectory !!!correction 19.07.2023

            flagPhotonProduced = SetPhotonProductionParameters(etotal,mass);

            // correction 19.07.2023 fItrajectories += 1; //count this trajectory

            //reinitialize intermediate integrals fFa, fSs, fSc, fSsx, fSsy, fScx, fScy;
            //reset radiation integral internal variables to defaults;
            //reset the trajectory and radiation probability along the trajectory
            ResetRadIntegral();
        }
    }

    return flagPhotonProduced;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BaierKatkov::GeneratePhoton(G4FastStep &fastStep)
{
    const G4DynamicParticle theGamma = G4DynamicParticle(G4Gamma::Gamma(),
                                                         PhMomentumDirection,fEph0);
    //generation of a secondary photon
    fastStep.CreateSecondaryTrack(theGamma,fNewParticleCoordinateXYZ,fNewGlobalTime,true);
}
