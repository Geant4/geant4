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
/// \file electromagnetic/TestEm7/src/G4ScreenedNuclearRecoil.cc
/// \brief Implementation of the G4ScreenedNuclearRecoil class
//
//
//
// Class Description
// Process for screened electromagnetic nuclear elastic scattering; 
// Physics comes from:
// Marcus H. Mendenhall and Robert A. Weller, 
// "Algorithms  for  the rapid  computation  of  classical  cross  
// sections  for  screened  Coulomb  collisions  "
// Nuclear  Instruments  and  Methods  in  Physics  Research B58 (1991)  11-17  
// The only input required is a screening function phi(r/a) which is the ratio
// of the actual interatomic potential for two atoms with atomic 
// numbers Z1 and Z2,
// to the unscreened potential Z1*Z2*e^2/r where e^2 is elm_coupling in 
// Geant4 units
//
// First version, April 2004, Marcus H. Mendenhall, Vanderbilt University
//
// 5 May, 2004, Marcus Mendenhall
// Added an option for enhancing hard collisions statistically, to allow 
// backscattering calculations to be carried out with much improved event rates,
// without distorting the multiple-scattering broadening too much.
// the method SetCrossSectionHardening(G4double fraction, G4double 
//                                     HardeningFactor)
// sets what fraction of the events will be randomly hardened,
// and the factor by which the impact area is reduced for such selected events.
//
// 21 November, 2004, Marcus Mendenhall
// added static_nucleus to IsApplicable
// 
// 7 December, 2004, Marcus Mendenhall
// changed mean free path of stopping particle from 0.0 to 1.0*nanometer
// to avoid new verbose warning about 0 MFP in 4.6.2p02
// 
// 17 December, 2004, Marcus Mendenhall
// added code to permit screening out overly close collisions which are 
// expected to be hadronic, not Coulombic
//
// 19 December, 2004, Marcus Mendenhall
// massive rewrite to add modular physics stages and plug-in cross section table
// computation.  This allows one to select (e.g.) between the normal external 
// python process and an embedded python interpreter (which is much faster) 
// for generating the tables.
// It also allows one to switch between sub-sampled scattering (event biasing) 
// and normal scattering, and between non-relativistic kinematics and 
// relativistic kinematic approximations, without having a class for every 
// combination. Further, one can add extra stages to the scattering, which can 
// implement various book-keeping processes.
// 
// January 2007, Marcus Mendenhall
// Reorganized heavily for inclusion in Geant4 Core.  All modules merged into 
// one source and header, all historic code removed.
// 
// Class Description - End

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <stdio.h>

#include "globals.hh"

#include "G4ScreenedNuclearRecoil.hh"

const char* G4ScreenedCoulombCrossSectionInfo::CVSFileVers() { return 
 "G4ScreenedNuclearRecoil.cc,v 1.57 2008/05/07 11:51:26 marcus Exp GEANT4 tag ";
}

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4DataVector.hh"
#include "G4Track.hh"
#include "G4Step.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ElementVector.hh"
#include "G4IsotopeVector.hh"

#include "G4EmProcessSubType.hh"

#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ProcessManager.hh"
#include "G4StableIsotopes.hh"
#include "G4LindhardPartition.hh"

#include "Randomize.hh"

#include <iostream>
#include <iomanip>

#include "c2_factory.hh"
static c2_factory<G4double> c2; // this makes a lot of notation shorter
typedef c2_ptr<G4double> c2p;

G4ScreenedCoulombCrossSection::~G4ScreenedCoulombCrossSection()
{
        screeningData.clear();
        MFPTables.clear();        
}

const G4double G4ScreenedCoulombCrossSection::massmap[nMassMapElements+1]={
  0, 1.007940, 4.002602, 6.941000, 9.012182, 10.811000, 12.010700, 
  14.006700, 15.999400, 18.998403, 20.179700, 22.989770, 24.305000, 26.981538, 
  28.085500, 
  30.973761, 32.065000, 35.453000, 39.948000, 39.098300, 40.078000, 44.955910, 
  47.867000, 
  50.941500, 51.996100, 54.938049, 55.845000, 58.933200, 58.693400, 63.546000, 
  65.409000, 
  69.723000, 72.640000, 74.921600, 78.960000, 79.904000, 83.798000, 85.467800, 
  87.620000, 
  88.905850, 91.224000, 92.906380, 95.940000, 98.000000, 101.070000, 102.905500,
  106.420000, 
  107.868200, 112.411000, 114.818000, 118.710000, 121.760000, 127.600000, 
  126.904470, 131.293000, 
  132.905450, 137.327000, 138.905500, 140.116000, 140.907650, 144.240000, 
  145.000000, 150.360000, 
  151.964000, 157.250000, 158.925340, 162.500000, 164.930320, 167.259000, 
  168.934210, 173.040000, 
  174.967000, 178.490000, 180.947900, 183.840000, 186.207000, 190.230000, 
  192.217000, 195.078000, 
  196.966550, 200.590000, 204.383300, 207.200000, 208.980380, 209.000000, 
  210.000000, 222.000000, 
  223.000000, 226.000000, 227.000000, 232.038100, 231.035880, 238.028910, 
  237.000000, 244.000000, 
  243.000000, 247.000000, 247.000000, 251.000000, 252.000000, 257.000000, 
  258.000000, 259.000000, 
  262.000000, 261.000000, 262.000000, 266.000000, 264.000000, 277.000000, 
  268.000000, 281.000000, 
  272.000000, 285.000000, 282.500000, 289.000000, 287.500000, 292.000000};

G4ParticleDefinition* 
G4ScreenedCoulombCrossSection::SelectRandomUnweightedTarget(
  const G4MaterialCutsCouple* couple)
{
  // Select randomly an element within the material, according to number 
  // density only        
        const G4Material* material = couple->GetMaterial();
        G4int nMatElements = material->GetNumberOfElements();
        const G4ElementVector* elementVector = material->GetElementVector();
        const G4Element *element=0;
        G4ParticleDefinition*target=0;
        
        // Special case: the material consists of one element
        if (nMatElements == 1)
    {
                element= (*elementVector)[0];
    }
        else
    {
      // Composite material
      G4double random = G4UniformRand() * material->GetTotNbOfAtomsPerVolume();
      G4double nsum=0.0;
      const G4double *atomDensities=material->GetVecNbOfAtomsPerVolume();
                
                for (G4int k=0 ; k < nMatElements ; k++ )
        {
                        nsum+=atomDensities[k];
                        element= (*elementVector)[k];
                        if (nsum >= random) break;
        }
    }

        G4int N=0;
        G4int Z=element->GetZasInt();
        
        G4int nIsotopes=element->GetNumberOfIsotopes();
        if(0<nIsotopes) {
          if(Z<=92) {
            // we have no detailed material isotopic info available, 
            // so use G4StableIsotopes table up to Z=92
            static G4StableIsotopes theIso; 
            // get a stable isotope table for default results
            nIsotopes=theIso.GetNumberOfIsotopes(Z);
            G4double random = 100.0*G4UniformRand(); 
            // values are expressed as percent, sum is 100
            G4int tablestart=theIso.GetFirstIsotope(Z);
            G4double asum=0.0;
            for(G4int i=0; i<nIsotopes; i++) {
              asum+=theIso.GetAbundance(i+tablestart);
              N=theIso.GetIsotopeNucleonCount(i+tablestart);
              if(asum >= random) break;
            }
          } else {
            // too heavy for stable isotope table, just use mean mass
            N=(G4int)std::floor(element->GetN()+0.5);
          }
        } else {
                G4int i;
                const G4IsotopeVector *isoV=element->GetIsotopeVector();
                G4double random = G4UniformRand();
                G4double *abundance=element->GetRelativeAbundanceVector();
                G4double asum=0.0;
                for(i=0; i<nIsotopes; i++) {
                        asum+=abundance[i];
                        N=(*isoV)[i]->GetN();
                        if(asum >= random) break;
                }
        }
        
        // get the official definition of this nucleus, to get the correct 
        // value of A note that GetIon is very slow, so we will cache ones 
        // we have already found ourselves.
        ParticleCache::iterator p=targetMap.find(Z*1000+N);
        if (p != targetMap.end()) {
                target=(*p).second;
        } else{
                target=G4IonTable::GetIonTable()->GetIon(Z, N, 0.0); 
                targetMap[Z*1000+N]=target;
        }
        return target;
}

void G4ScreenedCoulombCrossSection::BuildMFPTables()
{
        const G4int nmfpvals=200;

        std::vector<G4double> evals(nmfpvals), mfpvals(nmfpvals);

        // sum up inverse MFPs per element for each material        
        const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        if (materialTable == 0) { return; }
        //G4Exception("G4ScreenedCoulombCrossSection::BuildMFPTables 
        //- no MaterialTable found)");
        
        G4int nMaterials = G4Material::GetNumberOfMaterials();

        for (G4int matidx=0; matidx < nMaterials; matidx++) {

          const G4Material* material= (*materialTable)[matidx];
          const G4ElementVector &elementVector = 
            *(material->GetElementVector());
          const G4int nMatElements = material->GetNumberOfElements();

          const G4Element *element=0;
          const G4double *atomDensities=material->GetVecNbOfAtomsPerVolume();
                
          G4double emin=0, emax=0; 
          // find innermost range of cross section functions
          for (G4int kel=0 ; kel < nMatElements ; kel++ )
            {
                        element=elementVector[kel];
                        G4int Z=(G4int)std::floor(element->GetZ()+0.5);
                        const G4_c2_function &ifunc=sigmaMap[Z];
                        if(!kel || ifunc.xmin() > emin) emin=ifunc.xmin();
                        if(!kel || ifunc.xmax() < emax) emax=ifunc.xmax();
            }
                
          G4double logint=std::log(emax/emin) / (nmfpvals-1) ; 
          // logarithmic increment for tables
                
          // compute energy scale for interpolator.  Force exact values at 
          // both ends to avoid range errors
          for (G4int i=1; i<nmfpvals-1; i++) evals[i]=emin*std::exp(logint*i);
          evals.front()=emin;
          evals.back()=emax;
                
          // zero out the inverse mfp sums to start
          for (G4int eidx=0; eidx < nmfpvals; eidx++) mfpvals[eidx] = 0.0; 
                
          // sum inverse mfp for each element in this material and for each 
          // energy
          for (G4int kel=0 ; kel < nMatElements ; kel++ )
            {
              element=elementVector[kel];
              G4int Z=(G4int)std::floor(element->GetZ()+0.5);
              const G4_c2_function &sigma=sigmaMap[Z];
              G4double ndens = atomDensities[kel]; 
              // compute atom fraction for this element in this material
                                                
              for (G4int eidx=0; eidx < nmfpvals; eidx++) {
                mfpvals[eidx] += ndens*sigma(evals[eidx]);
              }
            }
                
          // convert inverse mfp to regular mfp
          for (G4int eidx=0; eidx < nmfpvals; eidx++) {
            mfpvals[eidx] = 1.0/mfpvals[eidx];
          }
          // and make a new interpolating function out of the sum
          MFPTables[matidx] = c2.log_log_interpolating_function().load(evals, 
                              mfpvals,true,0,true,0);
        }
}

G4ScreenedNuclearRecoil::
G4ScreenedNuclearRecoil(const G4String& processName, 
                        const G4String &ScreeningKey,
                        G4bool GenerateRecoils, 
                        G4double RecoilCutoff, G4double PhysicsCutoff) : 
  G4VDiscreteProcess(processName, fElectromagnetic),
        screeningKey(ScreeningKey),
        generateRecoils(GenerateRecoils), avoidReactions(1), 
        recoilCutoff(RecoilCutoff), physicsCutoff(PhysicsCutoff),
        hardeningFraction(0.0), hardeningFactor(1.0),
        externalCrossSectionConstructor(0),
        NIELPartitionFunction(new G4LindhardRobinsonPartition)
{
  // for now, point to class instance of this. Doing it by creating a new 
  // one fails
  // to correctly update NIEL
  // not even this is needed... done in G4VProcess().
  // pParticleChange=&aParticleChange; 
  processMaxEnergy=50000.0*MeV;
  highEnergyLimit=100.0*MeV;
  lowEnergyLimit=physicsCutoff;
  registerDepositedEnergy=1; // by default, don't hide NIEL
  MFPScale=1.0;
  // SetVerboseLevel(2);
  AddStage(new G4ScreenedCoulombClassicalKinematics);
  AddStage(new G4SingleScatter); 
  SetProcessSubType(fCoulombScattering);
}

void G4ScreenedNuclearRecoil::ResetTables()
{
        
  std::map<G4int, G4ScreenedCoulombCrossSection*>::iterator xt=
    crossSectionHandlers.begin();
        for(;xt != crossSectionHandlers.end(); xt++) {
                delete (*xt).second;
        }
        crossSectionHandlers.clear();
}

void G4ScreenedNuclearRecoil::ClearStages()
{
  // I don't think I like deleting the processes here... they are better 
  // abandoned
  // if the creator doesn't get rid of them
  // std::vector<G4ScreenedCollisionStage *>::iterator stage=
  //collisionStages.begin();
  //for(; stage != collisionStages.end(); stage++) delete (*stage);

        collisionStages.clear();
}

void G4ScreenedNuclearRecoil::SetNIELPartitionFunction(
     const G4VNIELPartition *part)
{
        if(NIELPartitionFunction) delete NIELPartitionFunction;
        NIELPartitionFunction=part;
}

void G4ScreenedNuclearRecoil::DepositEnergy(G4int z1, G4double a1, 
  const G4Material *material, G4double energy)
{
        if(!NIELPartitionFunction) {
                IonizingLoss+=energy;
        } else {
          G4double part=NIELPartitionFunction->PartitionNIEL(z1, a1, 
            material, energy);
          IonizingLoss+=energy*(1-part);
          NIEL += energy*part;
        }
}

G4ScreenedNuclearRecoil::~G4ScreenedNuclearRecoil()
{
        ResetTables();
}

// returns true if it appears the nuclei collided, and we are interested 
// in checking
G4bool G4ScreenedNuclearRecoil::CheckNuclearCollision(
                        G4double A, G4double a1, G4double apsis) {
  return avoidReactions && (apsis < (1.1*(std::pow(A,1.0/3.0)+
                                          std::pow(a1,1.0/3.0)) + 1.4)*fermi);
  // nuclei are within 1.4 fm (reduced pion Compton wavelength) of each 
  // other at apsis, 
  // this is hadronic, skip it
}

G4ScreenedCoulombCrossSection 
*G4ScreenedNuclearRecoil::GetNewCrossSectionHandler(void) {
  G4ScreenedCoulombCrossSection *xc;
  if(!externalCrossSectionConstructor) 
    xc=new G4NativeScreenedCoulombCrossSection;
        else xc=externalCrossSectionConstructor->create();
        xc->SetVerbosity(verboseLevel);
        return xc;
}

G4double G4ScreenedNuclearRecoil::GetMeanFreePath(const G4Track& track,
                                                  G4double, 
                                                  G4ForceCondition* cond)
{
        const G4DynamicParticle* incoming = track.GetDynamicParticle();
        G4double energy = incoming->GetKineticEnergy();
        G4double a1=incoming->GetDefinition()->GetPDGMass()/amu_c2;
        
        G4double meanFreePath;
        *cond=NotForced;
        
        if (energy < lowEnergyLimit || energy < recoilCutoff*a1) {
                *cond=Forced;
                return 1.0*nm; 
/* catch and stop slow particles to collect their NIEL! */
        } else if (energy > processMaxEnergy*a1) {
                return DBL_MAX; // infinite mean free path
        } else if (energy > highEnergyLimit*a1) energy=highEnergyLimit*a1; 
/* constant MFP at high energy */

        G4double fz1=incoming->GetDefinition()->GetPDGCharge();
        G4int z1=(G4int)(fz1/eplus + 0.5);

        std::map<G4int, G4ScreenedCoulombCrossSection*>::iterator xh=
                crossSectionHandlers.find(z1);
        G4ScreenedCoulombCrossSection *xs;
        
        if (xh==crossSectionHandlers.end()) {
                xs =crossSectionHandlers[z1]=GetNewCrossSectionHandler();
                xs->LoadData(screeningKey, z1, a1, physicsCutoff);
                xs->BuildMFPTables();
        } else xs=(*xh).second;
        
        const G4MaterialCutsCouple* materialCouple = 
          track.GetMaterialCutsCouple();
        size_t materialIndex = materialCouple->GetMaterial()->GetIndex();

        const G4_c2_function &mfp=*(*xs)[materialIndex];
        
        // make absolutely certain we don't get an out-of-range energy
        meanFreePath = mfp(std::min(std::max(energy, mfp.xmin()), mfp.xmax()));
        
        // G4cout << "MFP: " << meanFreePath << " index " << materialIndex 
        //<< " energy " << energy << " MFPScale " << MFPScale << G4endl;
        
        return meanFreePath*MFPScale;
}

G4VParticleChange* G4ScreenedNuclearRecoil::PostStepDoIt(
  const G4Track& aTrack, const G4Step& aStep)
{
        validCollision=1;
        pParticleChange->Initialize(aTrack);
        NIEL=0.0; // default is no NIEL deposited
        IonizingLoss=0.0;
        
        // do universal setup
        
        const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
        G4ParticleDefinition *baseParticle=aTrack.GetDefinition();
        
        G4double fz1=baseParticle->GetPDGCharge()/eplus;
        G4int z1=(G4int)(fz1+0.5);
        G4double a1=baseParticle->GetPDGMass()/amu_c2;
        G4double incidentEnergy = incidentParticle->GetKineticEnergy();
                
        // Select randomly one element and (possibly) isotope in the 
        // current material.
        const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
        
        const G4Material* mat = couple->GetMaterial();

        G4double P=0.0; // the impact parameter of this collision

        if(incidentEnergy < GetRecoilCutoff()*a1) { 
          // check energy sanity on entry
          DepositEnergy(z1, baseParticle->GetPDGMass()/amu_c2, mat, 
                        incidentEnergy);
          GetParticleChange().ProposeEnergy(0.0);
          // stop the particle and bail out
          validCollision=0;
        } else {
                
                G4double numberDensity=mat->GetTotNbOfAtomsPerVolume();
                G4double lattice=0.5/std::pow(numberDensity,1.0/3.0); 
                // typical lattice half-spacing
                G4double length=GetCurrentInteractionLength();
                G4double sigopi=1.0/(pi*numberDensity*length);  
                // this is sigma0/pi
                
                // compute the impact parameter very early, so if is rejected 
                // as too far away, little effort is wasted
                // this is the TRIM method for determining an impact parameter 
                // based on the flight path
                // this gives a cumulative distribution of 
                // N(P)= 1-exp(-pi P^2 n l)
                // which says the probability of NOT hitting a disk of area 
                // sigma= pi P^2 =exp(-sigma N l) 
                // which may be reasonable
                if(sigopi < lattice*lattice) { 
                        // normal long-flight approximation
                        P = std::sqrt(-std::log(G4UniformRand()) *sigopi);
                } else {
                        // short-flight limit
                        P = std::sqrt(G4UniformRand())*lattice;
                }
                        
                G4double fraction=GetHardeningFraction();
                if(fraction && G4UniformRand() < fraction) { 
                  // pick out some events, and increase the central cross 
                  // section by reducing the impact parameter
                  P /= std::sqrt(GetHardeningFactor());
                }
                
                
                // check if we are far enough away that the energy transfer 
                // must be below cutoff, 
                // and leave everything alone if so, saving a lot of time.
                if(P*P > sigopi) {
                  if(GetVerboseLevel() > 1) 
    printf("ScreenedNuclear impact reject: length=%.3f P=%.4f limit=%.4f\n",
           length/angstrom, P/angstrom,std::sqrt(sigopi)/angstrom);
                  // no collision, don't follow up with anything
                  validCollision=0;
                }
        }
        
        // find out what we hit, and record it in our kinematics block.
        kinematics.targetMaterial=mat;
        kinematics.a1=a1;
        
        if(validCollision) {
                G4ScreenedCoulombCrossSection *xsect=
                  GetCrossSectionHandlers()[z1];
                G4ParticleDefinition *recoilIon=
                        xsect->SelectRandomUnweightedTarget(couple);
                kinematics.crossSection=xsect;
                kinematics.recoilIon=recoilIon;
                kinematics.impactParameter=P;
                kinematics.a2=recoilIon->GetPDGMass()/amu_c2;
        } else {
                kinematics.recoilIon=0;
                kinematics.impactParameter=0;
                kinematics.a2=0;
        }
        
        std::vector<G4ScreenedCollisionStage *>::iterator stage=
          collisionStages.begin();
        
        for(; stage != collisionStages.end(); stage++) 
                (*stage)->DoCollisionStep(this,aTrack, aStep);
        
        if(registerDepositedEnergy) {
                pParticleChange->ProposeLocalEnergyDeposit(IonizingLoss+NIEL);
                pParticleChange->ProposeNonIonizingEnergyDeposit(NIEL);
                //MHM G4cout << "depositing energy, total = " 
                //<< IonizingLoss+NIEL << " NIEL = " << NIEL << G4endl;
        }

        return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

G4ScreenedCoulombClassicalKinematics::G4ScreenedCoulombClassicalKinematics() :
// instantiate all the needed functions statically, so no allocation is 
// done at run time
// we will be solving x^2 - x phi(x*au)/eps - beta^2 == 0.0
// or, for easier scaling, x'^2 - x' au phi(x')/eps - beta^2 au^2
// note that only the last of these gets deleted, since it owns the rest
  phifunc(c2.const_plugin_function()),
  xovereps(c2.linear(0., 0., 0.)), 
  // will fill this in with the right slope at run time
  diff(c2.quadratic(0., 0., 0., 1.)-xovereps*phifunc)
{
}

G4bool G4ScreenedCoulombClassicalKinematics::DoScreeningComputation(
  G4ScreenedNuclearRecoil *master,
  const G4ScreeningTables *screen, G4double eps, G4double beta)
{
  G4double au=screen->au;
  G4CoulombKinematicsInfo &kin=master->GetKinematics();
  G4double A=kin.a2;
  G4double a1=kin.a1;
        
  G4double xx0; // first estimate of closest approach
  if(eps < 5.0) {
    G4double y=std::log(eps);
    G4double mlrho4=((((3.517e-4*y+1.401e-2)*y+2.393e-1)*y+2.734)*y+2.220);
    G4double rho4=std::exp(-mlrho4); // W&M eq. 18
    G4double bb2=0.5*beta*beta;
    xx0=std::sqrt(bb2+std::sqrt(bb2*bb2+rho4)); // W&M eq. 17
  } else {
    G4double ee=1.0/(2.0*eps);
    xx0=ee+std::sqrt(ee*ee+beta*beta); // W&M eq. 15 (Rutherford value)  
    if(master->CheckNuclearCollision(A, a1, xx0*au)) return 0; 
    // nuclei too close
                
  }
        
  // we will be solving x^2 - x phi(x*au)/eps - beta^2 == 0.0
  // or, for easier scaling, x'^2 - x' au phi(x')/eps - beta^2 au^2
  xovereps.reset(0., 0.0, au/eps); // slope of x*au/eps term
  phifunc.set_function(&(screen->EMphiData.get())); 
  // install interpolating table
  G4double xx1, phip, phip2;
  G4int root_error;        
  xx1=diff->find_root(phifunc.xmin(), std::min(10*xx0*au,phifunc.xmax()), 
                      std::min(xx0*au, phifunc.xmax()), beta*beta*au*au, 
                      &root_error, &phip, &phip2)/au; 
        
        if(root_error) {
                G4cout << "Screened Coulomb Root Finder Error" << G4endl;
                G4cout << "au " << au << " A " << A << " a1 " << a1 
                       << " xx1 " << xx1 << " eps " << eps 
                       << " beta " << beta << G4endl;
                G4cout << " xmin " << phifunc.xmin() << " xmax " 
                       << std::min(10*xx0*au,phifunc.xmax()) ;
                G4cout << " f(xmin) " << phifunc(phifunc.xmin()) 
                       <<  " f(xmax) " 
                       << phifunc(std::min(10*xx0*au,phifunc.xmax())) ;
                G4cout << " xstart " << std::min(xx0*au, phifunc.xmax()) 
                       <<  " target " <<  beta*beta*au*au ;
                G4cout << G4endl;
                throw c2_exception("Failed root find");
        }
        
        // phiprime is scaled by one factor of au because phi is evaluated 
        // at (xx0*au),
        G4double phiprime=phip*au;
        
        //lambda0 is from W&M 19
        G4double lambda0=1.0/std::sqrt(0.5+beta*beta/(2.0*xx1*xx1)
                                       -phiprime/(2.0*eps));
        
        // compute the 6-term Lobatto integral alpha (per W&M 21, with 
        // different coefficients)
        // this is probably completely un-needed but gives the highest 
        // quality results,
        G4double alpha=(1.0+ lambda0)/30.0;
        G4double xvals[]={0.98302349, 0.84652241, 0.53235309, 0.18347974};
        G4double weights[]={0.03472124, 0.14769029, 0.23485003, 0.18602489};
        for(G4int k=0; k<4; k++) {
                G4double x, ff;
                x=xx1/xvals[k];
                ff=1.0/std::sqrt(1.0-phifunc(x*au)/(x*eps)-beta*beta/(x*x));
                alpha+=weights[k]*ff;
        }
        
        phifunc.unset_function(); 
        // throws an exception if used without setting again

        G4double thetac1=pi*beta*alpha/xx1; 
        // complement of CM scattering angle
        G4double sintheta=std::sin(thetac1); //note sin(pi-theta)=sin(theta)
        G4double costheta=-std::cos(thetac1); // note cos(pi-theta)=-cos(theta)
        // G4double psi=std::atan2(sintheta, costheta+a1/A); 
        // lab scattering angle (M&T 3rd eq. 8.69)
        
        // numerics note:  because we checked above for reasonable values 
        // of beta which give real recoils,
        // we don't have to look too closely for theta -> 0 here 
        // (which would cause sin(theta)
        // and 1-cos(theta) to both vanish and make the atan2 ill behaved).
        G4double zeta=std::atan2(sintheta, 1-costheta); 
        // lab recoil angle (M&T 3rd eq. 8.73)
        G4double coszeta=std::cos(zeta);
        G4double sinzeta=std::sin(zeta);

        kin.sinTheta=sintheta;
        kin.cosTheta=costheta;
        kin.sinZeta=sinzeta;
        kin.cosZeta=coszeta;
        return 1; // all OK, collision is valid
}

void G4ScreenedCoulombClassicalKinematics::DoCollisionStep(
  G4ScreenedNuclearRecoil *master,
  const G4Track& aTrack, const G4Step&) {
        
        if(!master->GetValidCollision()) return;
        
        G4ParticleChange &aParticleChange=master->GetParticleChange();
        G4CoulombKinematicsInfo &kin=master->GetKinematics();
        
        const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
        G4ParticleDefinition *baseParticle=aTrack.GetDefinition();
        
        G4double incidentEnergy = incidentParticle->GetKineticEnergy();
        
        // this adjustment to a1 gives the right results for soft 
        // (constant gamma)
        // relativistic collisions.  Hard collisions are wrong anyway, since the
        // Coulombic and hadronic terms interfere and cannot be added.
        G4double gamma=(1.0+incidentEnergy/baseParticle->GetPDGMass());
        G4double a1=kin.a1*gamma; // relativistic gamma correction
        
        G4ParticleDefinition *recoilIon=kin.recoilIon;
        G4double A=recoilIon->GetPDGMass()/amu_c2;
        G4int Z=(G4int)((recoilIon->GetPDGCharge()/eplus)+0.5);
        
        G4double Ec = incidentEnergy*(A/(A+a1)); 
        // energy in CM frame (non-relativistic!)
        const G4ScreeningTables *screen=kin.crossSection->GetScreening(Z);
        G4double au=screen->au; // screening length
        
        G4double beta = kin.impactParameter/au; 
        // dimensionless impact parameter
        G4double eps = Ec/(screen->z1*Z*elm_coupling/au); 
        // dimensionless energy
        
        G4bool ok=DoScreeningComputation(master, screen, eps, beta);        
        if(!ok) {
                master->SetValidCollision(0); // flag bad collision
                return; // just bail out without setting valid flag
        }
        
        G4double eRecoil=4*incidentEnergy*a1*A*kin.cosZeta*kin.cosZeta
          /((a1+A)*(a1+A));
        kin.eRecoil=eRecoil;
        
        if(incidentEnergy-eRecoil < master->GetRecoilCutoff()*a1) {
                aParticleChange.ProposeEnergy(0.0);
                master->DepositEnergy(int(screen->z1), a1, kin.targetMaterial, 
                                      incidentEnergy-eRecoil);
        } 
        
        if(master->GetEnableRecoils() && 
           eRecoil > master->GetRecoilCutoff() * kin.a2) {
                kin.recoilIon=recoilIon;                
        } else {
                kin.recoilIon=0; // this flags no recoil to be generated
                master->DepositEnergy(Z, A, kin.targetMaterial, eRecoil) ;
        }        
}

void G4SingleScatter::DoCollisionStep(G4ScreenedNuclearRecoil *master,
                                      const G4Track& aTrack, const G4Step&) {
        
        if(!master->GetValidCollision()) return;

        G4CoulombKinematicsInfo &kin=master->GetKinematics();
        G4ParticleChange &aParticleChange=master->GetParticleChange();

        const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
        G4double incidentEnergy = incidentParticle->GetKineticEnergy();
        G4double eRecoil=kin.eRecoil;
        
        G4double azimuth=G4UniformRand()*(2.0*pi);
        G4double sa=std::sin(azimuth);
        G4double ca=std::cos(azimuth);

        G4ThreeVector recoilMomentumDirection(kin.sinZeta*ca, 
                                              kin.sinZeta*sa, kin.cosZeta);
        G4ParticleMomentum incidentDirection = 
          incidentParticle->GetMomentumDirection();
        recoilMomentumDirection=
          recoilMomentumDirection.rotateUz(incidentDirection);
        G4ThreeVector recoilMomentum=
          recoilMomentumDirection*std::sqrt(2.0*eRecoil*kin.a2*amu_c2);
        
        if(aParticleChange.GetEnergy() != 0.0) { 
          // DoKinematics hasn't stopped it!
          G4ThreeVector beamMomentum=
            incidentParticle->GetMomentum()-recoilMomentum;
          aParticleChange.ProposeMomentumDirection(beamMomentum.unit()) ;
          aParticleChange.ProposeEnergy(incidentEnergy-eRecoil);
        }
        
        if(kin.recoilIon) {
                G4DynamicParticle* recoil = 
                  new G4DynamicParticle (kin.recoilIon,
                                recoilMomentumDirection,eRecoil) ;
                
                aParticleChange.SetNumberOfSecondaries(1);
                aParticleChange.AddSecondary(recoil);   
        }
}

G4bool G4ScreenedNuclearRecoil::
IsApplicable(const G4ParticleDefinition& aParticleType)
{
        return  aParticleType == *(G4Proton::Proton()) ||
        aParticleType.GetParticleType() == "nucleus" ||
        aParticleType.GetParticleType() == "static_nucleus";
}

void 
G4ScreenedNuclearRecoil::
BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  G4String nam = aParticleType.GetParticleName();
  if(nam == "GenericIon" || nam == "proton" 
     || nam == "deuteron" || nam == "triton" 
     || nam == "alpha" || nam == "He3") {
    G4cout << G4endl << GetProcessName() << ":   for  " << nam
           << "    SubType= " << GetProcessSubType() 
           << "    maxEnergy(MeV)= " << processMaxEnergy/MeV << G4endl;
  }
}

void 
G4ScreenedNuclearRecoil::
DumpPhysicsTable(const G4ParticleDefinition&)
{
}

// This used to be the file mhmScreenedNuclearRecoil_native.cc
// it has been included here to collect this file into a smaller 
// number of packages

#include "G4DataVector.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ElementVector.hh"
#include <vector>

G4_c2_function &ZBLScreening(G4int z1, G4int z2, size_t npoints, 
                             G4double rMax, G4double *auval)
{
        static const size_t ncoef=4;
        static G4double scales[ncoef]={-3.2, -0.9432, -0.4028, -0.2016};
        static G4double coefs[ncoef]={0.1818,0.5099,0.2802,0.0281};
        
        G4double au=
          0.8854*angstrom*0.529/(std::pow(z1, 0.23)+std::pow(z2,0.23));
        std::vector<G4double> r(npoints), phi(npoints);
        
        for(size_t i=0; i<npoints; i++) {
                G4double rr=(float)i/(float)(npoints-1);
                r[i]=rr*rr*rMax; 
                // use quadratic r scale to make sampling fine near the center
                G4double sum=0.0;
                for(size_t j=0; j<ncoef; j++) 
                  sum+=coefs[j]*std::exp(scales[j]*r[i]/au);
                phi[i]=sum;
        }

        // compute the derivative at the origin for the spline
        G4double phiprime0=0.0;
        for(size_t j=0; j<ncoef; j++) 
          phiprime0+=scales[j]*coefs[j]*std::exp(scales[j]*r[0]/au);
        phiprime0*=(1.0/au); // put back in natural units;
        
        *auval=au;
        return c2.lin_log_interpolating_function().load(r, phi, false, 
                                                        phiprime0,true,0);
}

G4_c2_function &MoliereScreening(G4int z1, G4int z2, size_t npoints, 
                                 G4double rMax, G4double *auval)
{        
        static const size_t ncoef=3;
        static G4double scales[ncoef]={-6.0, -1.2, -0.3};
        static G4double coefs[ncoef]={0.10, 0.55, 0.35};
        
        G4double au=0.8853*0.529*angstrom/std::sqrt(std::pow(z1, 0.6667)
                                                    +std::pow(z2,0.6667));
        std::vector<G4double> r(npoints), phi(npoints);
        
        for(size_t i=0; i<npoints; i++) {
                G4double rr=(float)i/(float)(npoints-1);
                r[i]=rr*rr*rMax; 
                // use quadratic r scale to make sampling fine near the center
                G4double sum=0.0;
                for(size_t j=0; j<ncoef; j++) 
                  sum+=coefs[j]*std::exp(scales[j]*r[i]/au);
                phi[i]=sum;
        }
        
        // compute the derivative at the origin for the spline
        G4double phiprime0=0.0;
        for(size_t j=0; j<ncoef; j++) 
          phiprime0+=scales[j]*coefs[j]*std::exp(scales[j]*r[0]/au);
        phiprime0*=(1.0/au); // put back in natural units;
        
        *auval=au;
        return c2.lin_log_interpolating_function().load(r, phi, false, 
                                                        phiprime0,true,0);
}

G4_c2_function &LJScreening(G4int z1, G4int z2, size_t npoints, 
                            G4double rMax, G4double *auval)
{        
//from Loftager, Besenbacher, Jensen & Sorensen
//PhysRev A20, 1443++, 1979
        G4double au=0.8853*0.529*angstrom/std::sqrt(std::pow(z1, 0.6667)
                                                    +std::pow(z2,0.6667));
        std::vector<G4double> r(npoints), phi(npoints);
        
        for(size_t i=0; i<npoints; i++) {
                G4double rr=(float)i/(float)(npoints-1);
                r[i]=rr*rr*rMax; 
                // use quadratic r scale to make sampling fine near the center

                G4double y=std::sqrt(9.67*r[i]/au);
                G4double ysq=y*y;
                G4double phipoly=1+y+0.3344*ysq+0.0485*y*ysq+0.002647*ysq*ysq;
                phi[i]=phipoly*std::exp(-y);
                // G4cout << r[i] << " " << phi[i] << G4endl;
        }

        // compute the derivative at the origin for the spline
        G4double logphiprime0=(9.67/2.0)*(2*0.3344-1.0); 
        // #avoid 0/0 on first element
        logphiprime0 *= (1.0/au); // #put back in natural units
        
        *auval=au;
        return c2.lin_log_interpolating_function().load(r, phi, false, 
                                                        logphiprime0*phi[0],
                                                        true,0);
}

G4_c2_function &LJZBLScreening(G4int z1, G4int z2, size_t npoints, 
                               G4double rMax, G4double *auval)
{        
// hybrid of LJ and ZBL, uses LJ if x < 0.25*auniv, ZBL if x > 1.5*auniv, and 
/// connector in between.  These numbers are selected so the switchover
// is very near the point where the functions naturally cross.
        G4double auzbl, aulj;
        
        c2p zbl=ZBLScreening(z1, z2, npoints, rMax, &auzbl);
        c2p lj=LJScreening(z1, z2, npoints, rMax, &aulj);

        G4double au=(auzbl+aulj)*0.5;
        lj->set_domain(lj->xmin(), 0.25*au);
        zbl->set_domain(1.5*au,zbl->xmax());
        
        c2p conn=
          c2.connector_function(lj->xmax(), lj, zbl->xmin(), zbl, true,0);
        c2_piecewise_function_p<G4double> &pw=c2.piecewise_function();
        c2p keepit(pw);
        pw.append_function(lj);
        pw.append_function(conn);
        pw.append_function(zbl);
        
        *auval=au;
        keepit.release_for_return();
        return pw;
}

G4NativeScreenedCoulombCrossSection::~G4NativeScreenedCoulombCrossSection() {
}

G4NativeScreenedCoulombCrossSection::G4NativeScreenedCoulombCrossSection() {
        AddScreeningFunction("zbl", ZBLScreening);
        AddScreeningFunction("lj", LJScreening);        
        AddScreeningFunction("mol", MoliereScreening);        
        AddScreeningFunction("ljzbl", LJZBLScreening);        
}

std::vector<G4String> 
G4NativeScreenedCoulombCrossSection::GetScreeningKeys() const {
  std::vector<G4String> keys;
  // find the available screening keys
  std::map<std::string, ScreeningFunc>::const_iterator sfunciter=phiMap.begin();
  for(; sfunciter != phiMap.end(); sfunciter++) 
    keys.push_back((*sfunciter).first);
  return keys;
}

static inline G4double cm_energy(G4double a1, G4double a2, G4double t0) {
        // "relativistically correct energy in CM frame"
        G4double m1=a1*amu_c2, mass2=a2*amu_c2;
        G4double mc2=(m1+mass2);
        G4double f=2.0*mass2*t0/(mc2*mc2);
        // old way: return (f < 1e-6) ?  0.5*mc2*f : mc2*(std::sqrt(1.0+f)-1.0);
        // formally equivalent to previous, but numerically stable for all 
        // f without conditional
        // uses identity (sqrt(1+x) - 1)(sqrt(1+x) + 1) = x
        return mc2*f/(std::sqrt(1.0+f)+1.0); 
}

static inline G4double  thetac(G4double m1, G4double mass2, G4double eratio) {
        G4double s2th2=eratio*( (m1+mass2)*(m1+mass2)/(4.0*m1*mass2) );
        G4double sth2=std::sqrt(s2th2);
        return 2.0*std::asin(sth2);
}

void G4NativeScreenedCoulombCrossSection::LoadData(G4String screeningKey, 
                                                   G4int z1, G4double a1, 
                                                   G4double recoilCutoff)
{                        
        static const size_t sigLen=200; 
        // since sigma doesn't matter much, a very coarse table will do 
        G4DataVector energies(sigLen);
        G4DataVector data(sigLen);
        
        a1=standardmass(z1); 
        // use standardized values for mass for building tables
        
        const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4int nMaterials = G4Material::GetNumberOfMaterials();
        
        for (G4int im=0; im<nMaterials; im++)
          {
            const G4Material* material= (*materialTable)[im];
            const G4ElementVector* elementVector = material->GetElementVector();
            const G4int nMatElements = material->GetNumberOfElements();
                
            for (G4int iEl=0; iEl<nMatElements; iEl++)
              {
                const G4Element* element = (*elementVector)[iEl];
                G4int Z = element->GetZasInt();
                G4double a2=element->GetA()*(mole/gram);
                        
                if(sigmaMap.find(Z)!=sigmaMap.end()) continue; 
                // we've already got this element
                        
                // find the screening function generator we need
                std::map<std::string, ScreeningFunc>::iterator sfunciter=
                  phiMap.find(screeningKey);
                if(sfunciter==phiMap.end()) {
                  G4ExceptionDescription ed;
                  ed << "No such screening key <" 
                     << screeningKey << ">"; 
                  G4Exception("G4NativeScreenedCoulombCrossSection::LoadData",
                              "em0003",FatalException,ed);
                }
                ScreeningFunc sfunc=(*sfunciter).second;
                        
                G4double au; 
                G4_c2_ptr screen=sfunc(z1, Z, 200, 50.0*angstrom, &au); 
                // generate the screening data
                G4ScreeningTables st;
        
                st.EMphiData=screen; //save our phi table
                st.z1=z1; st.m1=a1; st.z2=Z; st.m2=a2; st.emin=recoilCutoff;
                st.au=au;

                // now comes the hard part... build the total cross section
                // tables from the phi table                        
                // based on (pi-thetac) = pi*beta*alpha/x0, but noting that 
                // alpha is very nearly unity, always
                // so just solve it wth alpha=1, which makes the solution 
                // much easier
                // this function returns an approximation to 
                // (beta/x0)^2=phi(x0)/(eps*x0)-1 ~ ((pi-thetac)/pi)^2
                // Since we don't need exact sigma values, this is good enough 
                // (within a factor of 2 almost always)
                // this rearranges to phi(x0)/(x0*eps) = 
                // 2*theta/pi - theta^2/pi^2
                        
                c2_linear_p<G4double> &c2eps=c2.linear(0.0, 0.0, 1.0); 
                // will store an appropriate eps inside this in loop
                G4_c2_ptr phiau=screen(c2.linear(0.0, 0.0, au));
                G4_c2_ptr x0func(phiau/c2eps); 
                // this will be phi(x)/(x*eps) when c2eps is correctly set
                x0func->set_domain(1e-6*angstrom/au, 0.9999*screen->xmax()/au); 
                // needed for inverse function
                // use the c2_inverse_function interface for the root finder
                // it is more efficient for an ordered 
                // computation of values.
                G4_c2_ptr x0_solution(c2.inverse_function(x0func));
                        
                G4double m1c2=a1*amu_c2;
                G4double escale=z1*Z*elm_coupling/au; 
                // energy at screening distance
                G4double emax=m1c2; 
                // model is doubtful in very relativistic range
                G4double eratkin=0.9999*(4*a1*a2)/((a1+a2)*(a1+a2)); 
                // #maximum kinematic ratio possible at 180 degrees
                G4double cmfact0=st.emin/cm_energy(a1, a2, st.emin);
                G4double l1=std::log(emax);
                G4double l0=std::log(st.emin*cmfact0/eratkin);

                if(verbosity >=1) 
                  G4cout << "Native Screening: " << screeningKey << " " 
                         << z1 << " " << a1 << " " << 
                    Z << " " << a2 << " " << recoilCutoff << G4endl;
                        
                for(size_t idx=0; idx<sigLen; idx++) {
                  G4double ee=std::exp(idx*((l1-l0)/sigLen)+l0);
                  G4double gamma=1.0+ee/m1c2;
                  G4double eratio=(cmfact0*st.emin)/ee; 
                  // factor by which ee needs to be reduced to get emin
                  G4double theta=thetac(gamma*a1, a2, eratio);
                        
                  G4double eps=cm_energy(a1, a2, ee)/escale; 
                  // #make sure lab energy is converted to CM for these 
                  // calculations
                  c2eps.reset(0.0, 0.0, eps); 
                  // set correct slope in this function
                                
                  G4double q=theta/pi;
                  // G4cout << ee << " " << m1c2 << " " << gamma << " " 
                  // << eps << " " << theta << " " << q << G4endl;
                  // old way using root finder
                  // G4double x0= x0func->find_root(1e-6*angstrom/au, 
                  // 0.9999*screen.xmax()/au, 1.0, 2*q-q*q);
                  // new way using c2_inverse_function which caches 
                  // useful information so should be a bit faster
                  // since we are scanning this in strict order.
                  G4double x0=0;
                  try {
                    x0=x0_solution(2*q-q*q);
                  } catch(c2_exception&) {
                    G4Exception("G4ScreenedNuclearRecoil::LoadData",
                      "em0003",FatalException,
                      "failure in inverse solution to generate MFP tables");
                  }
                  G4double betasquared=x0*x0 - x0*phiau(x0)/eps;        
                  G4double sigma=pi*betasquared*au*au;
                  energies[idx]=ee;
                  data[idx]=sigma;
                }
                screeningData[Z]=st;
                sigmaMap[Z] = 
                  c2.log_log_interpolating_function().load(energies, data, 
                                                           true,0,true,0);
              }
          }
}
