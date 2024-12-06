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
// Author:      Alexei Sytov
// Co-author:   Gianfranco Paterno (testing)
// Using the key points of G4BaierKatkov and developments of V.V. Tikhomirov,
// partially described in L. Bandiera et al. Eur. Phys. J. C 82, 699 (2022)

#include "G4CoherentPairProduction.hh"

#include "Randomize.hh"
#include "G4TouchableHistory.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"

#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4EmProcessSubType.hh"
#include "G4TransportationManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CoherentPairProduction::G4CoherentPairProduction(const G4String& aName,
                                                   G4ProcessType):
    G4VDiscreteProcess(aName)
{
    SetProcessSubType(fCoherentPairProduction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CoherentPairProduction::GetMeanFreePath(const G4Track& aTrack,
                                    G4double,
                                    G4ForceCondition* condition)
{
    //current logical volume
    G4LogicalVolume* crystallogic;

    //momentum direction and coordinates (see comments below)
    G4ThreeVector momentumDirectionGamma,xyzGamma0,xyzGamma;
    //angle of the photon in the local reference system of the volume
    G4double txGamma0 = 0, tyGamma0 = 0;

    *condition = NotForced;

    //model activation
    G4bool modelTrigger = false;

    //photon energy
    G4double eGamma = aTrack.GetTotalEnergy();

    //energy cut, at the beginning, to not check everything else
    if(eGamma > ModelMinPrimaryEnergy())
    {
        //current logical volume
        crystallogic = aTrack.GetVolume()->GetLogicalVolume();

        //the model works only in the G4Region fG4RegionName
        if(crystallogic->GetRegion()->GetName()==fG4RegionName)
        {
            fCrystalData->SetGeometryParameters(crystallogic);

            //the momentum direction of the photon in the local reference system of the volume
            momentumDirectionGamma =
                (aTrack.GetTouchableHandle()->GetHistory()->
                    GetTopTransform().NetRotation().inverse())*aTrack.GetMomentumDirection();

            //the coordinates of the photon in the local reference system of the volume
            xyzGamma0 =
                aTrack.GetTouchableHandle()->GetHistory()->
                    GetTopTransform().TransformPoint(aTrack.GetPosition());

            // the coordinates of the photon in the co-rotating reference system within
            //a channel (elementary periodic cell)
            xyzGamma = fCrystalData->CoordinatesFromBoxToLattice(xyzGamma0);

            //angle of the photon in the local reference system of the volume
            //(!!! ONLY FORWARD DIRECTION, momentumDirectionGamma.getZ()>0,
            txGamma0 = std::atan(momentumDirectionGamma.x()/momentumDirectionGamma.z());
            tyGamma0 = std::atan(momentumDirectionGamma.y()/momentumDirectionGamma.z());

            //recalculate angle into the lattice reference system
            G4double angle = fCrystalData->AngleXFromBoxToLattice(txGamma0,xyzGamma.z());
            if (fCrystalData->GetModel()==2)
            {
                angle = std::sqrt(angle*angle+tyGamma0*tyGamma0);
            }

            //Applies the parameterisation not at the last step, only forward local direction
            //above low energy limit and below angular limit
            modelTrigger = (momentumDirectionGamma.z()>0. &&
                            std::abs(angle) < GetHighAngleLimit());
        }
    }

    if(modelTrigger)
    {
        //execute the model

        G4double x=0.,y=0.,z=0.;// the coordinates of charged particles
            //in the co-rotating reference system within
            //a channel (elementary periodic cell)
        G4double tx0=0.,ty0=0.; // the angles of charged particles
            // in the local reference system of the volume
        G4double txPreStep0=0.,tyPreStep0=0.; // the same as tx0, ty0 before the step
            // in the co-rotating reference system within
            //a channel (elementary periodic cell)

        G4ThreeVector scatteringAnglesAndEnergyLoss;//output of scattering functions

        //coordinates in Runge-Kutta calculations
        G4double x1=0.,x2=0.,x3=0.,x4=0.,y1=0.,y2=0.,y3=0.,y4=0.;
        //angles in Runge-Kutta calculations
        G4double tx1=0.,tx2=0.,tx3=0.,tx4=0.,ty1=0.,ty2=0.,ty3=0.,ty4=0.;
        //variables in Runge-Kutta calculations
        G4double kvx1=0.,kvx2=0.,kvx3=0.,kvx4=0.,kvy1=0.,kvy2=0.,kvy3=0.,kvy4=0.;
        //simulation step along z (internal step of the model) and its parts
        G4double dz=0.,dzd3=0.,dzd8=0.;//dzd3 = dz/3; dzd8 = dz/8;
        //simulation step along the momentum direction
        G4double momentumDirectionStep;
        //effective simulation step (taking into account nuclear density along the trajectory)
        G4double effectiveStep=0.;

        // Baier-Katkov variables
        G4double dzMeV=0.; //step in MeV^-1
        G4double axt=0.,ayt=0.; //charged particle accelerations
        G4double vxin=0.,vyin=0.;//the angles vs the photon (with incoherent scattering)
        G4double vxno=0.,vyno=0.;//the angles vs the photon (without incoherent scattering)

        G4double dzmod=0.;
        G4double fa1=0.,faseBefore=0.,faseBeforedz=0.,faseBeforedzd2=0.;
        G4double faseAfter=0.,fa2dfaseBefore2=0.;

        G4double skJ=0, skIx=0., skIy=0.;
        G4double sinfa1=0.,cosfa1=0.;

        //2-vector is needed for an initial parameter collection of 1 pair
        //vector of 2-vectors is an initial parameter collection of all sampling pair

        //collection of etotal for a single pair
        CLHEP::Hep2Vector twoVectorEtotal(0.,0.);

        //collection of x for a single pair
        CLHEP::Hep2Vector twoVectorX(0.,0.);

        //collection of y for a single pair
        CLHEP::Hep2Vector twoVectorY(0.,0.);

        //collection of tx for a single pair
        CLHEP::Hep2Vector twoVectorTX(0.,0.);

        //collection of tx for a single pair
        CLHEP::Hep2Vector twoVectorTY(0.,0.);

        fullVectorEtotal.clear();
        fullVectorX.clear();
        fullVectorY.clear();
        fullVectorTX.clear();
        fullVectorTY.clear();
        fPairProductionCDFdz.clear();
        fPairProductionCDFdz.push_back(0.);//0th element equal to 0

        const G4double charge[2] = {-1.,1.}; //particle charge
        const G4String particleName[2] = {"e-", "e+"};

        // the coordinates of a charged particle in the reference system within
        //a channel (elementary periodic cell)
        G4ThreeVector xyzparticle = xyzGamma;//changed below

        //the idea of pair production simulation is analogical to radiation in G4BaierKatkov
        //since the matrix element of these processes is the same => we solve inverse problem
        //to radiation: sample the pairs, calculate their trajectories and then calculate the
        //probabilities using Baier-Katkov analogically to radiation

        //cycle by sampling e+- pairs
        for(G4int i=0; i<fNMCPairs;i++)
        {
            //pair energy uniform sampling
            G4double etotal = fMass + fPPKineticEnergyCut +
                              G4UniformRand()*(eGamma-2*(fMass+fPPKineticEnergyCut));//particle
                                                                                     //total energy

            G4double phi = CLHEP::twopi*G4UniformRand();//necessary for pair kinematics

            //the probability of the production of the current pair (will be simulated)
            //per distance
            G4double probabilityPPdz = 0.;

            //cycle e- and e+ within single pair
            for(G4int j=0; j<2;j++)
            {
                if(j==1){etotal=eGamma-etotal;} //2nd particle energy
                twoVectorEtotal[j]=etotal;

                //Baier-Katkov input
                //intermediate variables to reduce calculations (the same names as in G4BaierKatkov)
                G4double e2 = etotal*etotal;
                G4double gammaInverse2 = fMass*fMass/(etotal*etotal);// 1/gamma^2
                //normalization coefficient
                G4double coefNorm = CLHEP::fine_structure_const/(8*(CLHEP::pi2))/(2.*fNMCPairs);
                //G4double phi = CLHEP::twopi*G4UniformRand();//necessary for pair kinematics
                G4double om = eGamma;
                G4double eprime=om-etotal; //E'=omega-E
                G4double eprime2 = eprime*eprime;
                G4double e2pluseprime2 =e2+eprime2;
                G4double omprime=etotal*om/eprime;//om'=E*om/(om-E)
                G4double omprimed2=omprime/2;

                //difference vs G4BaierKatkov: om -> etotal
                G4double coefNorme2deprime2 = coefNorm*e2/eprime2; //e2/om/om;//e2/eprime2;

                G4double gammaInverse2om = gammaInverse2*om*om;

                //initialize intermediate integrals with zeros
                G4double fa=0.,ss=0.,sc=0.,ssx=0.,ssy=0.,scx=0.,scy=0.;

                //End of Baier-Katkov input

                G4bool fbreak = false;//flag of the trajectory cycle break

                //set fCrystalData parameters depending on the particle parameters
                fCrystalData->SetParticleProperties(etotal, fMass,
                                                    charge[j], particleName[j]);

                //needed just to setup the correct value of channel No in the crystal
                //since later it may be changed during the trajectory calculation
                fCrystalData->CoordinatesFromBoxToLattice(xyzGamma0);

                //coordinate sampling: random x and y due to coordinate uncertainty
                //in the interaction point
                if(j==0)
                {
                    x = fCrystalData->GetChannelWidthX()*G4UniformRand();
                    y = fCrystalData->GetChannelWidthY()*G4UniformRand();
                }
                else
                {
                    x=twoVectorX[0];
                    y=twoVectorY[0];
                }
                twoVectorX[j] = x;
                twoVectorY[j] = y;
                //definite z as a coordinate of the photon (uncertainty of the
                //interaction point is taking into account later by simulation
                //of the position of pair production)
                z = xyzGamma.z();

                //angles of the photon in the co-rotating reference system within a channel =>
                //angular distribution center
                G4double tx = fCrystalData->AngleXFromBoxToLattice(txGamma0,z);
                G4double ty = tyGamma0;
                G4double momentumDirectionZGamma = 1./
                                                   std::sqrt(1.+std::pow(std::tan(tx),2)+
                                                             std::pow(std::tan(ty),2));

                //angle sampling: depends on angular range within a particle trajectory
                //defined by the Lindhard angle and on the angle of radiation proportional
                //to 1/gamma

                //range of MC integration on angles
                G4double paramParticleAngle = fChargeParticleAngleFactor*fMass/etotal;

                G4double axangle=0.;
                if (fCrystalData->GetModel()==1)//1D model (only angle vs plane matters)
                {
                    axangle = std::abs(tx);
                }
                else if  (fCrystalData->GetModel()==2)//2D model
                {
                    axangle = std::sqrt(tx*tx+ty*ty);
                }

                if(axangle>fCrystalData->GetLindhardAngle()+DBL_EPSILON)
                {
                    paramParticleAngle+=axangle
                                          -std::sqrt(axangle*axangle
                                                      -fCrystalData->GetLindhardAngle()
                                                            *fCrystalData->GetLindhardAngle());
                }
                else
                {
                    paramParticleAngle+=fCrystalData->GetLindhardAngle();
                }


                //ONLY forward direction
                if (paramParticleAngle>CLHEP::halfpi-DBL_EPSILON){paramParticleAngle=CLHEP::halfpi;}

                G4double rho=1.;
                G4double rhocut=CLHEP::halfpi/paramParticleAngle;//radial angular cut of
                                                                 //the distribution
                G4double norm=std::atan(rhocut*rhocut)*
                                CLHEP::pi*paramParticleAngle*paramParticleAngle;


                //distribution with long tails (useful to not exclude particle angles
                //after a strong single scattering)
                //at ellipsescale < 1 => half of statistics
                do
                {
                    rho = std::sqrt(std::tan(CLHEP::halfpi*G4UniformRand()));
                }
                while (rho>rhocut);

                //normalization coefficient for intergration on angles of charged particles
                G4double angleNormCoef = (1.+rho*rho*rho*rho)*norm;

                tx+=charge[j]*paramParticleAngle*rho*std::cos(phi);
                twoVectorTX[j] = tx;
                ty+=charge[j]*paramParticleAngle*rho*std::sin(phi);
                twoVectorTY[j] = ty;

                G4double zalongGamma = 0;//necessary for renormalization of PP probability
                    //depending on the trajectory length along Gamma direction
                //starting the trajectory
                //here we don't care about the boundaries of the crystal volume
                //the trajectory is very short and the pair production probability obtained
                //in Baier-Katkov will be extrapolated to the real step inside the crystal volume
                for(G4int k=0; k<fNTrajectorySteps;k++)
                {
                    //back to the local reference system of the volume
                    txPreStep0 = fCrystalData->AngleXFromLatticeToBox(tx,z);
                    tyPreStep0 = ty;

                    dz = fCrystalData->GetSimulationStep(tx,ty);
                    dzd3=dz/3;
                    dzd8=dz/8;

                    //trajectory calculation:
                    //Runge-Cutt "3/8"
                    //fCrystalData->GetCurv(z)*fCrystalData->GetCorrectionZ() is due to dependence
                    //of the radius on x; GetCurv gets 1/R for the central ("central plane/axis")

                    //first step
                    kvx1=fCrystalData->Ex(x,y);
                    x1=x+tx*dzd3;
                    tx1=tx+(kvx1-fCrystalData->GetCurv(z)*fCrystalData->GetCorrectionZ())*dzd3;
                    if (fCrystalData->GetModel()==2)
                    {
                        kvy1=fCrystalData->Ey(x,y);
                        y1=y+ty*dzd3;
                        ty1=ty+kvy1*dzd3;
                    }

                    //second step
                    kvx2=fCrystalData->Ex(x1,y1);
                    x2=x-tx*dzd3+tx1*dz;
                    tx2=tx-(kvx1-fCrystalData->GetCurv(z)*fCrystalData->GetCorrectionZ())*dzd3+
                          (kvx2-fCrystalData->GetCurv(z)*fCrystalData->GetCorrectionZ())*dz;
                    if (fCrystalData->GetModel()==2)
                    {
                        kvy2=fCrystalData->Ey(x1,y1);
                        y2=y-ty*dzd3+ty1*dz;
                        ty2=ty-kvy1*dzd3+kvy2*dz;
                    }

                    //third step
                    kvx3=fCrystalData->Ex(x2,y2);
                    x3=x+(tx-tx1+tx2)*dz;
                    tx3=tx+(kvx1-kvx2+kvx3-
                                fCrystalData->GetCurv(z)*fCrystalData->GetCorrectionZ())*dz;
                    if (fCrystalData->GetModel()==2)
                    {
                        kvy3=fCrystalData->Ey(x2,y2);
                        y3=y+(ty-ty1+ty2)*dz;
                        ty3=ty+(kvy1-kvy2+kvy3)*dz;
                    }

                    //fourth step
                    kvx4=fCrystalData->Ex(x3,y3);
                    x4=x+(tx+3.*tx1+3.*tx2+tx3)*dzd8;
                    tx4=tx+(kvx1+3.*kvx2+3.*kvx3+kvx4)*dzd8-
                          fCrystalData->GetCurv(z)*fCrystalData->GetCorrectionZ()*dz;
                    if (fCrystalData->GetModel()==2)
                    {
                        kvy4=fCrystalData->Ey(x3,y3);
                        y4=y+(ty+3.*ty1+3.*ty2+ty3)*dzd8;
                        ty4=ty+(kvy1+3.*kvy2+3.*kvy3+kvy4)*dzd8;
                    }
                    else
                    {
                        y4 =y+ty*dz;
                        ty4=ty;
                    }

                    x=x4;
                    tx=tx4;
                    y=y4;
                    ty=ty4;

                    z+=dz*fCrystalData->GetCorrectionZ();//motion along the z coordinate
                        //("central plane/axis", no current plane/axis)

                    xyzparticle = fCrystalData->ChannelChange(x,y,z);
                    x=xyzparticle.x();
                    y=xyzparticle.y();
                    z=xyzparticle.z();

                    momentumDirectionStep =
                        dz*std::sqrt(1+std::pow(std::tan(tx),2)+std::pow(std::tan(ty),2));
                    zalongGamma += dz/momentumDirectionZGamma;

                    //default scattering and energy loss 0
                    scatteringAnglesAndEnergyLoss.set(0.,0.,0.);

                    if(fIncoherentScattering)
                    {
                        //calculate separately for each element of the crystal
                        for (G4int ii = 0; ii < fCrystalData->GetNelements(); ii++)
                        {
                           //effective step taking into account nuclear density along the trajectory
                            effectiveStep = momentumDirectionStep*
                                            fCrystalData->NuclearDensity(x,y,ii);
                            //Coulomb scattering on screened atomic potential
                            //(both multiple and single)
                            scatteringAnglesAndEnergyLoss +=
                                fCrystalData->CoulombAtomicScattering(effectiveStep,
                                                                      momentumDirectionStep,
                                                                      ii);
                        }
                        //electron scattering and coherent part of ionization energy losses
                        scatteringAnglesAndEnergyLoss += fCrystalData->CoulombElectronScattering(
                            fCrystalData->MinIonizationEnergy(x,y),
                            fCrystalData->ElectronDensity(x,y),
                            momentumDirectionStep);
                        tx += scatteringAnglesAndEnergyLoss.x();
                        ty += scatteringAnglesAndEnergyLoss.y();
                    }

                    //To avoid backward direction
                    if(std::abs(tx)>CLHEP::halfpi-DBL_EPSILON||
                        std::abs(ty)>CLHEP::halfpi-DBL_EPSILON)
                    {
                        G4cout << "Warning: particle angle is beyond +-pi/2 range => "
                                  "skipping the calculation of its probability" << G4endl;
                        fbreak = true;
                        break;
                    }

                    //**********Baier-Katkov start

                    //back to the local reference system of the volume
                    tx0 = fCrystalData->AngleXFromLatticeToBox(tx,z);
                    ty0 = ty;

                    dzMeV=momentumDirectionStep/CLHEP::hbarc;// in MeV^-1

                    // accelerations
                    axt=(tx0-scatteringAnglesAndEnergyLoss.x()-txPreStep0)/dzMeV;
                    ayt=(ty0-scatteringAnglesAndEnergyLoss.y()-tyPreStep0)/dzMeV;

                    //the angles vs the photon (with incoherent scattering)
                    vxin = tx0-txGamma0;
                    vyin = ty0-tyGamma0;
                    //the angles vs the photon (without incoherent scattering)
                    vxno = vxin-scatteringAnglesAndEnergyLoss.x();
                    vyno = vyin-scatteringAnglesAndEnergyLoss.y();

                    //phase difference before scattering
                    faseBefore=omprimed2*(gammaInverse2+vxno*vxno+vyno*vyno);//phi' t<ti//MeV

                    faseBeforedz = faseBefore*dzMeV;
                    faseBeforedzd2 = faseBeforedz/2.;
                    fa+=faseBeforedz; //
                    fa1=fa-faseBeforedzd2;//
                    dzmod=2*std::sin(faseBeforedzd2)/faseBefore;//MeV^-1

                    //phi''/faseBefore^2
                    fa2dfaseBefore2 = omprime*(axt*vxno+ayt*vyno)/(faseBefore*faseBefore);

                    //phase difference after scattering
                    faseAfter=omprimed2*(gammaInverse2+vxin*vxin+vyin*vyin);//phi' ti+O//MeV

                    skJ=1/faseAfter-1/faseBefore-fa2dfaseBefore2*dzmod;//MeV^-1
                    skIx=vxin/faseAfter-vxno/faseBefore+dzmod*(axt/faseBefore-
                                                                           vxno*fa2dfaseBefore2);
                    skIy=vyin/faseAfter-vyno/faseBefore+dzmod*(ayt/faseBefore-
                                                                           vyno*fa2dfaseBefore2);

                    sinfa1 = std::sin(fa1);
                    cosfa1 = std::cos(fa1);

                    ss+=sinfa1*skJ;//sum sin integral J of BK
                    sc+=cosfa1*skJ;//sum cos integral J of BK
                    ssx+=sinfa1*skIx;// sum sin integral Ix of BK
                    ssy+=sinfa1*skIy;// sum sin integral Iy of BK
                    scx+=cosfa1*skIx;// sum cos integral Ix of BK
                    scy+=cosfa1*skIy;// sum cos integral Iy of BK
                }

                //only of the trajectory cycle was not broken
                if(!fbreak)
                {
                    G4double i2=ssx*ssx+scx*scx+ssy*ssy+scy*scy;//MeV^-2
                    G4double j2=ss*ss+sc*sc;//MeV^-2

                    probabilityPPdz += coefNorme2deprime2*angleNormCoef*
                                       (i2*e2pluseprime2+j2*gammaInverse2om)/zalongGamma;
                }
            }

            //filling the CDF of probabilities of the production of sampling pairs
            fPairProductionCDFdz.push_back(fPairProductionCDFdz[i]+probabilityPPdz);
            //**********Baier-Katkov end

            //accumulation of initial parameters of sampling pairs
            fullVectorEtotal.push_back(twoVectorEtotal);
            fullVectorX.push_back(twoVectorX);
            fullVectorY.push_back(twoVectorY);
            fullVectorTX.push_back(twoVectorTX);
            fullVectorTY.push_back(twoVectorTY);
        }

        //photon mean free path
        //fPairProductionCDFdz.back() = full pair production probability
        //simulated for the current photon along photon direction
        G4double lMeanFreePath = 1/fPairProductionCDFdz.back();

        fEffectiveLrad = 7.*lMeanFreePath/9.;//only for scoring purpose

        return lMeanFreePath;
    }
    else
    {
        //dummy process, does not occur
        return DBL_MAX;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4CoherentPairProduction::PostStepDoIt(const G4Track& aTrack,
                                           const G4Step& aStep)
{
	//example with no physical sense
    aParticleChange.Initialize(aTrack);
    //G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();

    const G4ParticleDefinition* chargedParticleDefinition[2] =
        {G4Electron::Electron(),G4Positron::Positron()};

    // the coordinates of the photon in the local reference system of the volume
    G4ThreeVector xyzGamma0 =
        aTrack.GetTouchableHandle()->GetHistory()->
            GetTopTransform().TransformPoint(aTrack.GetPosition());

    // the coordinates of the photon in the co-rotating reference system within
    //a channel (elementary periodic cell)
    G4ThreeVector xyzGamma = fCrystalData->CoordinatesFromBoxToLattice(xyzGamma0);

    //global time
    G4double tGlobalGamma = aTrack.GetGlobalTime();

    G4double ksi1 = G4UniformRand()*fPairProductionCDFdz.back();

    //randomly choosing the pair to be produced from the sampling list
    //according to the probabilities calculated in the Baier-Katkov integral
    G4int ipair = FindVectorIndex(fPairProductionCDFdz,ksi1)-1;//index of
        //a pair produced

    // the coordinates of a charged particle in the reference system within
    //a channel (elementary periodic cell)
    G4ThreeVector xyzparticle;
    //cycle e- and e+ within single pair
    for(G4int j=0; j<2;j++)
    {
        xyzparticle.set(fullVectorX[ipair][j],fullVectorY[ipair][j],xyzGamma.z());

        //in the local reference system of the volume
        G4ThreeVector newParticleCoordinateXYZ =
                fCrystalData->CoordinatesFromLatticeToBox(xyzparticle);
        //the same in the global reference system
        newParticleCoordinateXYZ =
            aTrack.GetTouchableHandle()->GetHistory()->
                GetTopTransform().Inverse().TransformPoint(newParticleCoordinateXYZ);

        //back to the local reference system of the volume
        G4double tx0 = fCrystalData->AngleXFromLatticeToBox(fullVectorTX[ipair][j],xyzGamma.z());
        G4double ty0 = fullVectorTY[ipair][j];

        G4double momentumDirectionZ = 1./
                                      std::sqrt(1.+std::pow(std::tan(tx0),2)+
                                                std::pow(std::tan(ty0),2));

        //momentum direction vector of the charged particle produced
        //in the local reference system of the volume
        G4ThreeVector momentumDirectionParticle = G4ThreeVector(momentumDirectionZ*std::tan(tx0),
                                                                momentumDirectionZ*std::tan(ty0),
                                                                momentumDirectionZ);
        //the same in the global reference system
        momentumDirectionParticle =
            (aTrack.GetTouchableHandle()->GetHistory()->GetTopTransform().NetRotation()) *
                momentumDirectionParticle;

        G4DynamicParticle* chargedParticle =
            new G4DynamicParticle(chargedParticleDefinition[j],
                              momentumDirectionParticle,
                              fullVectorEtotal[ipair][j]-fMass);

        // Create the track for the secondary particle
        G4Track* secondaryTrack = new G4Track(chargedParticle,
                                              tGlobalGamma,
                                              newParticleCoordinateXYZ);
        secondaryTrack->SetTouchableHandle(aStep.GetPostStepPoint()->GetTouchableHandle());
        secondaryTrack->SetParentID(aTrack.GetTrackID());

        //generation of a secondary charged particle
        aParticleChange.AddSecondary(secondaryTrack);
    }

    //killing the photon
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4CoherentPairProduction::FindVectorIndex(std::vector<G4double> &myvector, G4double value)
{
    auto iteratorbegin = myvector.begin();
    auto iteratorend   = myvector.end();

    //vector index (for non precise values lower_bound gives upper value)
    auto loweriterator = std::lower_bound(iteratorbegin, iteratorend, value);
    //return the index of the vector element
    return (G4int)std::distance(iteratorbegin, loweriterator);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoherentPairProduction::Input(const G4Material *crystal,
                                     const G4String &lattice,
                                     const G4String &filePath)
{
    //initializing the class with containing all
    //the crystal material and crystal lattice data and
    //Channeling scattering and ionization processes
    fCrystalData = new G4ChannelingFastSimCrystalData();
    //setting all the crystal material and lattice data
    fCrystalData->SetMaterialProperties(crystal,lattice,filePath);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CoherentPairProduction::Input(const G4ChannelingFastSimCrystalData *crystalData)
{
    //setting the class with containing all
    //the crystal material and crystal lattice data and
    //Channeling scattering and ionization processes
    //fCrystalData = new G4ChannelingFastSimCrystalData();

    fCrystalData = const_cast<G4ChannelingFastSimCrystalData*>(crystalData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CoherentPairProduction::ProcessDescription(std::ostream& out) const
{
    out << "  Coherent pair production";
    G4VDiscreteProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
