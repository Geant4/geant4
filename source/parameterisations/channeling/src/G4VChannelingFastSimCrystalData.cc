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

#include "G4VChannelingFastSimCrystalData.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4VChannelingFastSimCrystalData::G4VChannelingFastSimCrystalData()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VChannelingFastSimCrystalData::~G4VChannelingFastSimCrystalData(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VChannelingFastSimCrystalData::SetGeometryParameters
                                     (const G4LogicalVolume *crystallogic)
{
    G4int crystalID = crystallogic->GetInstanceID();

    //set bending angle if the volume exists in the list, otherwise default = 0
    (fMapBendingAngle.count(crystalID) > 0)
    ? SetBendingAngle(fMapBendingAngle[crystalID],crystallogic)
    : SetBendingAngle(0.,crystallogic);

    //set miscut angle if the volume exists in the list, otherwise default = 0
    (fMapMiscutAngle.count(crystalID) > 0)
    ? SetMiscutAngle(fMapMiscutAngle[crystalID],crystallogic)
    : SetMiscutAngle(0.,crystallogic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VChannelingFastSimCrystalData::SetBendingAngle(G4double tetab,
						      const G4LogicalVolume* crystallogic)
{
    G4int crystalID = crystallogic->GetInstanceID();

    //set the bending angle for this logical volume
    fMapBendingAngle[crystalID]=tetab;

    G4ThreeVector limboxmin;//minimal limits of the box bounding the logical volume
    G4ThreeVector limboxmax;//maximal limits of the box bounding the logical volume
    //save the limits of the box bounding the logical volume
    crystallogic->GetSolid()->BoundingLimits(limboxmin,limboxmax);

    //bounding box half dimensions
    fHalfDimBoundingBox = (limboxmax-limboxmin)/2.;

    G4double lcr = limboxmax.getZ()-limboxmin.getZ();//crystal thickness

    fBendingAngle=std::abs(tetab);
    if (fBendingAngle<0.000001)//no bending less then 1 urad
    {
       fBent=0;
       fBendingAngle=0.;
       fBendingR=0.;//just for convenience (infinity in reality)
       fBending2R=0.;
       fBendingRsquare=0.;
       fCurv=0.;

       G4cout << "Channeling model: volume " << crystallogic->GetName() << G4endl;
       G4cout << "Warning: bending angle is lower than 1 urad => set to 0" << G4endl;
    }
    else
    {
       fBent=1;
       fBendingR=lcr/fBendingAngle;
       fBending2R=2.*fBendingR;
       fBendingRsquare=fBendingR*fBendingR;
       fCurv=1./fBendingR;

       if (tetab<0.)
       {
           G4cout << "Channeling model: volume " << crystallogic->GetName() << G4endl;
           G4cout << "Warning: bending angle is negative => set to be positive" << G4endl;
       }
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VChannelingFastSimCrystalData::SetMiscutAngle(G4double tetam,
						     const G4LogicalVolume *crystallogic)
{
    G4int crystalID = crystallogic->GetInstanceID();

    //set the bending angle for this logical volume
    fMapMiscutAngle[crystalID]=tetam;

    // fMiscutAngle>0: rotation of xz coordinate planes clockwise in the xz plane
    fMiscutAngle=tetam;
    if (std::abs(tetam)>1.*mrad)
    {
        G4cout << "Channeling model: volume " << crystallogic->GetName() << G4endl;
        G4cout << "Warning: miscut angle is higher than 1 mrad => " << G4endl;
        G4cout << "coordinate transformation routines may be unstable" << G4endl;
    }
    fCosMiscutAngle=std::cos(fMiscutAngle);
    fSinMiscutAngle=std::sin(fMiscutAngle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VChannelingFastSimCrystalData::SetParticleProperties(G4double etotal,
                                                         G4double mass,
                                                         G4double charge,
                                                         G4bool ifhadron)
{
    G4double teta1;
    fZ2=charge;
    G4double zz22=fZ2*fZ2;
    fHadron=ifhadron;

//     particle momentum and energy
       G4double t=etotal*etotal-mass*mass; //  economy of operations
       fPz=std::sqrt(t); //  momentum of particle
       fPV=t/etotal;    //  pv
       fBeta=fPz/etotal;   //  velocity/c
       fTetaL = std::sqrt(fVmax2/fPV); //Lindhard angle
       fChannelingStep = fChangeStep/fTetaL; //standard simulation step

//     Energy losses
       fV2 = fBeta*fBeta; // particle (velocity/c)^2
       fGamma = etotal/mass; //  Lorentz factor
       fMe2Gamma = 2*CLHEP::electron_mass_c2*fGamma;
//     max ionization losses
       fTmax = fMe2Gamma*fGamma*fV2/
               (CLHEP::electron_mass_c2/mass*CLHEP::electron_mass_c2/mass +
                1. + fMe2Gamma/mass);

       for(G4int i=0; i<fNelements; i++)
       {

//       minimal scattering angle by coulomb scattering on nuclei
//       defining by shielding by electrons
//       teta1=hdc/(fPz*fRF)*DSQRT(1.13D0+3.76D0*(alpha*fZ1*fZ2/fBeta)**2){ev*cm/(eV*cm)}
         teta1=fTeta10[i]*std::sqrt(1.13+fK40[i]*zz22/fV2); // /fPz later to speed up
                                                            //  the calculations

//       the coefficient for multiple scattering
         fBB[i]=teta1*teta1*fPu11[i];
         fE1XBbb[i]=expint(fBB[i]);
         fBBDEXP[i]=(1.+fBB[i])*std::exp(fBB[i]);
//       necessary for suppression of incoherent scattering
//       by the atomic correlations in crystals for single scattering on nucleus
//       (screened atomic potential): EXP(-(fPz*teta*fU1)**2)=EXP(-fPzu11*teta**2).GE.ksi
//         =>no scattering
         fPzu11[i]=fPu11[i]*fPz*fPz;

         teta1=teta1/fPz; //
         fTeta12[i]=teta1*teta1;
//       maximal scattering angle by coulomb scattering on nuclei
//       defining by nucleus radius
//       tetamax=hc/(fPz*1.D-6*fR0*fAN**(1.D0/3.D0))//  {Mev*fermi/(MeV*fermi)}
         G4double tetamax=fTetamax0[i]/fPz;
         fTetamax2[i]=tetamax*tetamax;
         fTetamax12[i]=fTeta12[i]+fTetamax2[i];

//       a cofficient in a formula for scattering (for high speed of simulation)
//       fK2=(fZ2*alpha*hdc)**2*4.*pi*fN0*(fZ1/fPV)**2
//       fK3=(fZ2*alpha*hdc)**2*4.*pi*fN0/(fPV)**2
         fK2[i]=fK20[i]*zz22/fPV/fPV;
       }

//     nuclear diffractive scattering angle
       //tetaQEL=1./sqrt(2.*(9.26-4.94/sqrt(fPz/GeV)+0.28*log(fPz/GeV)));

       fK3=fK30/fV2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VChannelingFastSimCrystalData::GetLindhardAngle(G4double etotal, G4double mass)
{
    G4double pv0 = etotal-mass*mass/etotal;
    return std::sqrt(2*fVmax/pv0); //Calculate the value of the Lindhard angle
                                   //(!!! the value for a straight crystal)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VChannelingFastSimCrystalData::GetLindhardAngle()
{
    return fTetaL; //return the Lindhard angle value calculated in SetParticleProperties
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VChannelingFastSimCrystalData::GetSimulationStep(G4double tx,G4double ty)
{
    G4double simulationstep;
    //find angle of particle w.r.t. the plane or axis
    G4double angle=0.;
    if (iModel==1)//1D model
    {
        angle = std::abs(tx);
    }
    else if  (iModel==2)//2D model
    {
        angle = std::sqrt(tx*tx+ty*ty);
    }

    //compare this angle with the Lindhard angle
    if (angle<fTetaL)
    {
        simulationstep = fChannelingStep;
    }
    else
    {
        simulationstep = fChangeStep/angle;
    }

    return simulationstep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VChannelingFastSimCrystalData::GetMaxSimulationStep(G4double etotal,
                                                            G4double mass)
{
    //standard value of step for channeling particles which is the maximal possible step
    return fChangeStep/GetLindhardAngle(etotal, mass);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4VChannelingFastSimCrystalData::CoulombAtomicScattering(
                                                           G4double effectiveStep,
                                                           G4double step,
                                                           G4int ielement)
{
        G4double tx = 0.;//horizontal scattering angle
        G4double ty = 0.;//vertical scattering angle

        G4double ksi=0.1;

//      calculation of the teta2-minimal possible angle of a single scattering
        G4double e1=fK2[ielement]*effectiveStep; //for high speed of a program
//      (real formula is (4*pi*fN0*wpl(x)*dz*fZ1*zz2*alpha*hdc/fPV)**2)
        G4double teta122=fTetamax12[ielement]/(ksi*fTetamax12[ielement]/e1+1.);
        // teta122=fTeta12+teta22=teta1^2+teta2^2

        G4double teta22;
        G4double t;
//      if the angle of a single scattering is less teta1 - minimal possible
//      angle of coulomb scattering defining by the electron shielding than
//      multiple scattering by both nuclei and electrons and electrons will not
//      occur => minimal possible angle of a single scattering is equal to teta1
        if (teta122<=fTeta12[ielement]*1.000125)
        {
          teta22=0.;
          teta122=fTeta12[ielement];
        }
        else
        {
          teta22=teta122-fTeta12[ielement];
          G4double aa=teta22/fTeta12[ielement];
          G4double aa1=1.+aa;

//        crystal, with scattering suppression
          G4double tetamsi=e1*(std::log(aa1)+
                           (1.-std::exp(-aa*fBB[ielement]))/aa1+
                           fBBDEXP[ielement]*
                           (expint(fBB[ielement]*aa1)-fE1XBbb[ielement]));

//        sumilation of multiple coulomb scattering by nuclei and electrons
//        for high speed of a program, real formula is
//        4*pi*fN0*wpl(x)*dz*(fZ1*zz2*alpha*hdc/fPV)**2*
//         *(ln(1+a)+(1-exp(-a*b))/(1+a)+(1+b)*exp(b)*(E1XB(b*(1+a))-E1XB(b)))

          ksi=G4UniformRand();
          t=std::sqrt(-tetamsi*std::log(ksi));

          ksi=G4UniformRand();

          tx+=t*std::cos(CLHEP::twopi*ksi);
          ty+=t*std::sin(CLHEP::twopi*ksi);

        }
//      simulation of single coulomb scattering by nuclei (with screened potential)
        G4double zss=0.;
        G4double dzss=step;

//      (calculation of a distance, at which another single scattering can happen)
        ksi=G4UniformRand();

        zss=-std::log(ksi)*step/(e1*(1./teta122-1./fTetamax12[ielement]));
        G4double tt;

//      At some step several single scattering can occur.
//      So,if the distance of the next scattering is less than the step,
//      another scattering can occur. If the distance of the next scattering
//      is less than the difference between the step and the distance of
//      the previous scattering, another scattering can occur. And so on, and so on.
//      In the cycle we simulate each of them. The cycle is finished, when
//      the remaining part of step is less than a distance of the next single scattering.
//********************************************
//      if at a step a single scattering occurs
        while (zss<dzss)
        {

//        simulation by Monte-Carlo of angles of single scattering
          ksi=G4UniformRand();

          tt=fTetamax12[ielement]/(1.+ksi*(fTetamax2[ielement]-teta22)/teta122)-
                  fTeta12[ielement];

          ksi=G4UniformRand();

//        suppression of incoherent scattering by the atomic correlations in crystals
          t=fPzu11[ielement]*tt;
          t=std::exp(-t);

          if (t<ksi) //if scattering takes place
          {
            //scattering angle
            t=std::sqrt(tt);
            ksi=G4UniformRand();

            tx+=t*std::cos(CLHEP::twopi*ksi);
            ty+=t*std::sin(CLHEP::twopi*ksi);
          }

          dzss-=zss;
//        (calculation of a distance, at which another single scattering can happen)
          ksi=G4UniformRand();

          zss=-std::log(ksi)*step/(e1*(1./teta122-1./fTetamax12[ielement]));
        }
//********************************************
        return G4ThreeVector(tx,ty,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4VChannelingFastSimCrystalData::CoulombElectronScattering(
                                                             G4double eMinIonization,
                                                             G4double electronDensity,
                                                             G4double step)
{

    G4double zss=0.;
    G4double dzss=step;
    G4double ksi = 0.;

    G4double tx = 0.;//horizontal scattering angle
    G4double ty = 0.;//vertical scattering angle
    G4double eloss = 0.;//energy loss

    // eMinIonization - minimal energy transfered to electron
    // a cut to reduce the number of calls of electron scattering
    // is needed only at low density regions, in many cases does not do anything at all
    if (eMinIonization<0.5*eV){eMinIonization=0.5*eV;}

    //  single scattering on electrons routine
    if ((eMinIonization<fTmax)&&(electronDensity>DBL_EPSILON))
    {

//    (calculation of a distance, at which another single scattering can happen)
//    simulation of scattering length (by the same way single scattering by nucleus
      ksi=G4UniformRand();

      zss=-1.0*std::log(ksi)/(fK3*electronDensity)/(1./eMinIonization-1./fTmax);

//********************************************
//    if at a step a single scattering occur
      while (zss<dzss)
      {
//      simulation by Monte-Carlo of angles of single scattering
        ksi=G4UniformRand();

//      energy transfered to electron
        G4double e1=eMinIonization/(1.-ksi*(1.-eMinIonization/fTmax));

//      scattering angle
        G4double t=std::sqrt(e1*(e1+2.*CLHEP::electron_mass_c2))/fPz;

        // energy losses
        if (fHadron) {eloss=e1;} // we don't calculate ionization losses for e+-

        ksi=G4UniformRand();

        tx+=t*std::cos(CLHEP::twopi*ksi);
        ty+=t*std::sin(CLHEP::twopi*ksi);

        dzss-=zss;
//      (calculation of a distance, at which another single scattering can happen)
//      simulation of scattering length
//      (by the same way single scattering by nucleus
        ksi=G4UniformRand();

        zss=-1.0*std::log(ksi)/(fK3*electronDensity)/(1./eMinIonization-1./fTmax);
      }
//********************************************
    }
    return G4ThreeVector(tx,ty,eloss);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VChannelingFastSimCrystalData::IonizationLosses(G4double dz,
                                                    G4int ielement)
{
    G4double elosses = 0.;
    if (fHadron) {elosses=fKD[ielement]/fV2*
                (std::log(fMe2Gamma*fV2/fI0[ielement]/fGamma) - fV2)*dz;}
    return elosses;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VChannelingFastSimCrystalData::expint(G4double X)
{
//      ============================================
//      Purpose: Compute exponential integral E1(x)
//      Input :  x  --- Argument of E1(x)
//      Output:  E1 --- E1(x)
//      ============================================

G4double E1, R, T, T0;
G4int M;

if (X==0)
{
    E1=1.e300;
}
else if (X<=1.)
{
   E1=1.;
   R=1.;


   for(int K=1; K<=25; K++)
   {
       R=-R*K*X/std::pow(K+1.,2.);
       E1=E1+R;
       if (std::abs(R)<=std::abs(E1)*1.0e-15) {break;}
   }

   E1=-0.5772156649015328-std::log(X)+X*E1;
}
else
{
   M=20+std::trunc(80.0/X);
   T0=0.;

   for(int K=M; K>=1; K--)
   {
      T0=K/(1.0+K/(X+T0));
   }

   T=1.0/(X+T0);
   E1=std::exp(-X)*T;
}

return E1;
}
