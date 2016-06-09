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
// $Id: G4ScreenedNuclearRecoil.cc,v 1.5 2008/01/14 12:11:39 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01-patch-02 $
//
//
// Class Description
// Process for screened electromagnetic nuclear elastic scattering; 
// Physics comes from:
// Marcus H. Mendenhall and Robert A. Weller, 
// "Algorithms  for  the rapid  computation  of  classical  cross  
// sections  for  screened  Coulomb  collisions  "
// Nuclear  Instruments  and  Methods  in  Physics  Research  B58  (1991)  11-17  
// The only input required is a screening function phi(r/a) which is the ratio
// of the actual interatomic potential for two atoms with atomic numbers Z1 and Z2,
// to the unscreened potential Z1*Z2*e^2/r where e^2 is elm_coupling in Geant4 units
//
// First version, April 2004, Marcus H. Mendenhall, Vanderbilt University
//
// 5 May, 2004, Marcus Mendenhall
// Added an option for enhancing hard collisions statistically, to allow 
// backscattering calculations to be carried out with much improved event rates,
// without distorting the multiple-scattering broadening too much.
// the method SetCrossSectionHardening(G4double fraction, G4double HardeningFactor)
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
// computation.  This allows one to select (e.g.) between the normal external python
// process and an embedded python interpreter (which is much faster) for generating
// the tables.
// It also allows one to switch between sub-sampled scattering (event biasing) and
// normal scattering, and between non-relativistic kinematics and relativistic
// kinematic approximations, without having a class for every combination. Further, one can
// add extra stages to the scattering, which can implement various book-keeping processes.
// 
// January 2007, Marcus Mendenhall
// Reorganized heavily for inclusion in Geant4 Core.  All modules merged into 
// one source and header, all historic code removed.
// 
// Class Description - End


#include <stdio.h>

#include "globals.hh"

#include "G4ScreenedNuclearRecoil.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DataVector.hh"
#include "G4Track.hh"
#include "G4Step.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ElementVector.hh"
#include "G4IsotopeVector.hh"

#include "G4RangeTest.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ProcessManager.hh"
#include "G4StableIsotopes.hh"

#include "Randomize.hh"

#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <iomanip>

G4ScreenedCoulombCrossSection::~G4ScreenedCoulombCrossSection()
{
	ScreeningMap::iterator tables=screeningData.begin();
	for (;tables != screeningData.end(); tables++) {
		delete (*tables).second.EMphiData;
	}
	screeningData.clear();
	
	std::map<G4int, c2_function<G4double> *>::iterator mfpit=MFPTables.begin();
	for (;mfpit != MFPTables.end(); mfpit++) {
		delete (*mfpit).second;
	}
	MFPTables.clear();
	
}

const G4double G4ScreenedCoulombCrossSection::massmap[nMassMapElements+1]={
	0, 1.007940, 4.002602, 6.941000, 9.012182, 10.811000, 12.010700, 
	14.006700, 15.999400, 18.998403, 20.179700, 22.989770, 24.305000, 26.981538, 28.085500, 
	30.973761, 32.065000, 35.453000, 39.948000, 39.098300, 40.078000, 44.955910, 47.867000, 
	50.941500, 51.996100, 54.938049, 55.845000, 58.933200, 58.693400, 63.546000, 65.409000, 
	69.723000, 72.640000, 74.921600, 78.960000, 79.904000, 83.798000, 85.467800, 87.620000, 
	88.905850, 91.224000, 92.906380, 95.940000, 98.000000, 101.070000, 102.905500, 106.420000, 
	107.868200, 112.411000, 114.818000, 118.710000, 121.760000, 127.600000, 126.904470, 131.293000, 
	132.905450, 137.327000, 138.905500, 140.116000, 140.907650, 144.240000, 145.000000, 150.360000, 
	151.964000, 157.250000, 158.925340, 162.500000, 164.930320, 167.259000, 168.934210, 173.040000, 
	174.967000, 178.490000, 180.947900, 183.840000, 186.207000, 190.230000, 192.217000, 195.078000, 
	196.966550, 200.590000, 204.383300, 207.200000, 208.980380, 209.000000, 210.000000, 222.000000, 
	223.000000, 226.000000, 227.000000, 232.038100, 231.035880, 238.028910, 237.000000, 244.000000, 
	243.000000, 247.000000, 247.000000, 251.000000, 252.000000, 257.000000, 258.000000, 259.000000, 
	262.000000, 261.000000, 262.000000, 266.000000, 264.000000, 277.000000, 268.000000, 281.000000, 
	272.000000, 285.000000, 282.500000, 289.000000, 287.500000, 292.000000};

G4ParticleDefinition* G4ScreenedCoulombCrossSection::SelectRandomUnweightedTarget(const G4MaterialCutsCouple* couple)
{
	// Select randomly an element within the material, according to number density only	
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
	G4int Z=(G4int)std::floor(element->GetZ()+0.5);
	
	G4int nIsotopes=element->GetNumberOfIsotopes();
	if(!nIsotopes) {
		if(Z<=92) {
			// we have no detailed material isotopic info available, 
			// so use G4StableIsotopes table up to Z=92
			static G4StableIsotopes theIso; // get a stable isotope table for default results
			nIsotopes=theIso.GetNumberOfIsotopes(Z);
			G4double random = 100.0*G4UniformRand(); // values are expressed as percent, sum is 100
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
	
	// get the official definition of this nucleus, to get the correct value of A
	// note that GetIon is very slow, so we will cache ones we have already found ourselves.
	ParticleCache::iterator p=targetMap.find(Z*1000+N);
	if (p != targetMap.end()) {
		target=(*p).second;
	} else{
		target=G4ParticleTable::GetParticleTable()->GetIon(Z, N, 0.0);	
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
	if (materialTable == 0)
		G4Exception("G4ScreenedCoulombCrossSection::BuildMFPTables - no MaterialTable found)");
	
	G4int nMaterials = G4Material::GetNumberOfMaterials();

	for (G4int matidx=0; matidx < nMaterials; matidx++) {

		const G4Material* material= (*materialTable)[matidx];
		const G4ElementVector &elementVector = *(material->GetElementVector());
		const G4int nMatElements = material->GetNumberOfElements();

		const G4Element *element=0;
		const G4double *atomDensities=material->GetVecNbOfAtomsPerVolume();
		
		G4double emin=0, emax=0; // find innermost range of cross section functions
		for (G4int kel=0 ; kel < nMatElements ; kel++ )
        {
			element=elementVector[kel];
			G4int Z=(G4int)std::floor(element->GetZ()+0.5);
			c2_function<G4double> &ifunc=*sigmaMap[Z];
			if(!kel || ifunc.xmin() > emin) emin=ifunc.xmin();
			if(!kel || ifunc.xmax() < emax) emax=ifunc.xmax();			
        }
		
		G4double logint=std::log(emax/emin) / (nmfpvals-1) ; // logarithmic increment for tables
		
		// compute energy scale for interpolator.  Force exact values at both ends to avoid range errors
		for (G4int i=1; i<nmfpvals-1; i++) evals[i]=emin*std::exp(logint*i);
		evals.front()=emin;
		evals.back()=emax;
		
		// zero out the inverse mfp sums to start
		for (G4int eidx=0; eidx < nmfpvals; eidx++) mfpvals[eidx] = 0.0; 
		
		// sum inverse mfp for each element in this material and for each energy
		for (G4int kel=0 ; kel < nMatElements ; kel++ )
        {
			element=elementVector[kel];
			G4int Z=(G4int)std::floor(element->GetZ()+0.5);
			c2_function<G4double> &sigma=*sigmaMap[Z];
			G4double ndens = atomDensities[kel]; // compute atom fraction for this element in this material
						
			for (G4int eidx=0; eidx < nmfpvals; eidx++) {
					mfpvals[eidx] += ndens*sigma(evals[eidx]);
			}
        }
		
		// convert inverse mfp to regular mfp
		for (G4int eidx=0; eidx < nmfpvals; eidx++) {
			mfpvals[eidx] = 1.0/mfpvals[eidx];
		}
		// and make a new interpolating function out of the sum
		MFPTables[matidx] = static_cast<c2_function<G4double> *>(new log_log_interpolating_function<G4double>(
				evals, mfpvals));
    }
	
#ifdef DEBUG	
	for (G4int matidx=0; matidx < nMaterials; matidx++) {
		const G4Material* material= (*materialTable)[matidx];
		G4cout << "***** MFP (1MeV) ***** " << material->GetName() << "  " << (*MFPTables[matidx])(1.0) << G4endl; 
	}
#endif
	
}

G4ScreenedNuclearRecoil::
G4ScreenedNuclearRecoil(const G4String& processName, 
							   const G4String &ScreeningKey,
							   G4bool GenerateRecoils, 
							   G4double RecoilCutoff, G4double PhysicsCutoff) : 
	G4VDiscreteProcess(processName),
	screeningKey(ScreeningKey),
	generateRecoils(GenerateRecoils), avoidReactions(1), 
	recoilCutoff(RecoilCutoff), physicsCutoff(PhysicsCutoff),
	hardeningFraction(0.0), hardeningFactor(1.0),
	externalCrossSectionConstructor(0)
{
		highEnergyLimit=100.0*MeV;
		lowEnergyLimit=physicsCutoff;
		registerDepositedEnergy=1; // by default, don't hide NIEL
		MFPScale=1.0;
		// SetVerboseLevel(2);
		AddStage(new G4ScreenedCoulombClassicalKinematics);
		AddStage(new G4SingleScatter); 
}

void G4ScreenedNuclearRecoil::ResetTables()
{
	std::map<G4int, c2_function<G4double>*>::iterator xh=meanFreePathTables.begin();
	for(;xh != meanFreePathTables.end(); xh++) {
		delete (*xh).second;
	}
	meanFreePathTables.clear();
	
	std::map<G4int, G4ScreenedCoulombCrossSection*>::iterator xt=crossSectionHandlers.begin();
	for(;xt != crossSectionHandlers.end(); xt++) {
		delete (*xt).second;
	}
	crossSectionHandlers.clear();
}

void G4ScreenedNuclearRecoil::ClearStages()
{
	// I don't think I like deleting the processes here... they are better abandoned
	// if the creator doesn't get rid of them
	// std::vector<G4ScreenedCollisionStage *>::iterator stage=collisionStages.begin();
	//for(; stage != collisionStages.end(); stage++) delete (*stage);

	collisionStages.clear();
}

G4ScreenedNuclearRecoil::~G4ScreenedNuclearRecoil()
{
	ResetTables();
}

// returns true if it appears the nuclei collided, and we are interested in checking
G4bool G4ScreenedNuclearRecoil::CheckNuclearCollision(
			G4double A, G4double a1, G4double apsis) {
	return avoidReactions && (apsis < (1.1*(std::pow(A,1.0/3.0)+std::pow(a1,1.0/3.0)) + 1.4)*fermi);
	// nuclei are within 1.4 fm (reduced pion Compton wavelength) of each other at apsis, 
	// this is hadronic, skip it
}

G4ScreenedCoulombCrossSection *G4ScreenedNuclearRecoil::GetNewCrossSectionHandler(void) {
	G4ScreenedCoulombCrossSection *xc;
	if(!externalCrossSectionConstructor) xc=new G4NativeScreenedCoulombCrossSection;
	else xc=externalCrossSectionConstructor->create();
	xc->SetVerbosity(verboseLevel);
	return xc;
}

G4double G4ScreenedNuclearRecoil::GetMeanFreePath(const G4Track& track,
						  G4double, 
						  G4ForceCondition*)
{
	const G4DynamicParticle* incoming = track.GetDynamicParticle();
	G4double energy = incoming->GetKineticEnergy();
	G4double a1=incoming->GetDefinition()->GetPDGMass()/amu_c2;
	
	G4double meanFreePath;
	
	if (energy < lowEnergyLimit || energy < recoilCutoff) return 1.0*nanometer; /* stop slow particles! */
	else if (energy > highEnergyLimit*a1) energy=highEnergyLimit*a1; /* constant MFP at high energy */

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
	
	const G4MaterialCutsCouple* materialCouple = track.GetMaterialCutsCouple();
	size_t materialIndex = materialCouple->GetMaterial()->GetIndex();

	c2_function<G4double> &mfp=*(*xs)[materialIndex];
	
	// make absolutely certain we don't get an out-of-range energy
	meanFreePath = mfp(std::min(std::max(energy, mfp.xmin()), mfp.xmax()));
	
	// G4cout << "MFP: " << meanFreePath << " index " << materialIndex << " energy " << energy << " MFPScale " << MFPScale << G4endl;
	
	return meanFreePath*MFPScale;
}

G4VParticleChange* G4ScreenedNuclearRecoil::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
	validCollision=1;
	aParticleChange.Initialize(aTrack);
	NIEL=0.0; // default is no NIEL deposited
	
	// do universal setup
	
	const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
	G4ParticleDefinition *baseParticle=aTrack.GetDefinition();
	
	G4double fz1=baseParticle->GetPDGCharge()/eplus;
	G4int z1=(G4int)(fz1+0.5);
	G4double incidentEnergy = incidentParticle->GetKineticEnergy();
		
	// Select randomly one element and (possibly) isotope in the current material.
	const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
	
	if(incidentEnergy < GetRecoilCutoff()) { // check energy sanity on entry
		if(!baseParticle->GetProcessManager()->
		   GetAtRestProcessVector()->size())
			aParticleChange.ProposeTrackStatus(fStopAndKill);
		else
			aParticleChange.ProposeTrackStatus(fStopButAlive);

		AddToNIEL(incidentEnergy);
		aParticleChange.ProposeEnergy(0.0);
		// stop the particle and bail out
		validCollision=0;
	} 	
	
	const G4Material* mat = couple->GetMaterial();
	G4double numberDensity=mat->GetTotNbOfAtomsPerVolume();
	G4double lattice=0.5/std::pow(numberDensity,1.0/3.0); // typical lattice half-spacing
	G4double length=GetCurrentInteractionLength();
	G4double sigopi=1.0/(CLHEP::pi*numberDensity*length);  // this is sigma0/pi
	
	// compute the impact parameter very early, so if is rejected as too far away, little effort is wasted
	// this is the TRIM method for determining an impact parameter based on the flight path
	// this gives a cumulative distribution of N(P)= 1-exp(-pi P^2 n l)
	// which says the probability of NOT hitting a disk of area sigma= pi P^2 =exp(-sigma N l) 
	// which may be reasonable
	G4double P;
	if(sigopi < lattice*lattice) { 
		// normal long-flight approximation
		P = std::sqrt(-std::log(G4UniformRand()) *sigopi);
	} else {
		// short-flight limit
		P = std::sqrt(G4UniformRand())*lattice;
	}
		
	G4double fraction=GetHardeningFraction();
	if(fraction && G4UniformRand() < fraction) { 
		// pick out some events, and increase the central cross section
		// by reducing the impact parameter
		P /= std::sqrt(GetHardeningFactor());
	}
	
	
	// check if we are far enough away that the energy transfer must be below cutoff, 
	// and leave everything alone if so, saving a lot of time.
	if(P*P > sigopi) {
		if(GetVerboseLevel() > 1) 
			printf("ScreenedNuclear impact reject: length=%.3f P=%.4f limit=%.4f\n",
				   length/angstrom, P/angstrom,std::sqrt(sigopi)/angstrom);
		// no collision, don't follow up with anything
		validCollision=0;
	}
	
	// find out what we hit, and record it in our kinematics block.
	if(validCollision) {
		G4ScreenedCoulombCrossSection	*xsect=GetCrossSectionHandlers()[z1];
		G4ParticleDefinition *recoilIon=
			xsect->SelectRandomUnweightedTarget(couple);
		kinematics.crossSection=xsect;
		kinematics.recoilIon=recoilIon;
		kinematics.impactParameter=P;
		kinematics.a1=baseParticle->GetPDGMass()/amu_c2;
		kinematics.a2=recoilIon->GetPDGMass()/amu_c2;
	} else {
		kinematics.recoilIon=0;
		kinematics.impactParameter=0;
		kinematics.a1=baseParticle->GetPDGMass()/amu_c2;
		kinematics.a2=0;
	}
	
	std::vector<G4ScreenedCollisionStage *>::iterator stage=collisionStages.begin();
	
	for(; stage != collisionStages.end(); stage++) 
		(*stage)->DoCollisionStep(this,aTrack, aStep);
	
	if(registerDepositedEnergy) {
	  aParticleChange.ProposeLocalEnergyDeposit(NIEL);
	  aParticleChange.ProposeNonIonizingEnergyDeposit(NIEL);
	}
	return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

G4bool G4ScreenedCoulombClassicalKinematics::DoScreeningComputation(G4ScreenedNuclearRecoil *master,
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
		if(master->CheckNuclearCollision(A, a1, xx0*au)) return 0; // nuclei too close
		
	}
	
	c2_function<G4double> &phiData=*(screen->EMphiData);
	// instantiate all the needed functions statically, so no allocation is done at run time
	// we will be solving x^2 - x phi(x*au)/eps - beta^2 == 0.0
	// or, for easier scaling, x'^2 - x' au phi(x')/eps - beta^2 au^2
	static c2_plugin_function<G4double> phifunc;
	static c2_quadratic<G4double> xsq(0., 0., 0., 1.); //  x^2
	static c2_linear<G4double> xovereps(0., 0., 0.); // will fill this in with the right slope at run time
	static c2_function<G4double> &xphi=xovereps*phifunc;
	static c2_function<G4double> &diff=xsq-xphi;
	
	xovereps.reset(0., 0.0, au/eps); // slope of x*au/eps term
	phifunc.set_function(phiData); // install interpolating table
	
	G4double xx1, phip, phip2;
	G4int root_error;
	
	xx1=diff.find_root(phiData.xmin(), std::min(10*xx0*au,phiData.xmax()), 
					   std::min(xx0*au, phiData.xmax()), beta*beta*au*au, &root_error, &phip, &phip2)/au; 
	
	if(root_error) {
		G4cout << "Screened Coulomb Root Finder Error" << G4endl;
		G4cout << "au " << au << " A " << A << " a1 " << a1 << " xx1 " << xx1 << " eps " << eps << " beta " << beta << G4endl;
		G4cout << " xmin " << phiData.xmin() << " xmax " << std::min(10*xx0*au,phiData.xmax()) ;
		G4cout << " f(xmin) " << phifunc(phiData.xmin()) <<  " f(xmax) " << phifunc(std::min(10*xx0*au,phiData.xmax())) ;
		G4cout << " xstart " << std::min(xx0*au, phiData.xmax()) <<  " target " <<  beta*beta*au*au ;
		G4cout << G4endl;
		throw c2_exception("Failed root find");
	}
	
	phifunc.unset_function(); // throws an exception if used without setting again
							  // phiprime is scaled by one factor of au because phi is evaluated at (xx0*au),
	G4double phiprime=phip*au;
	
	//lambda0 is from W&M 19
	G4double lambda0=1.0/std::sqrt(0.5+beta*beta/(2.0*xx1*xx1)-phiprime/(2.0*eps));
	
	//compute the 6-term Lobatto integral alpha (per W&M 21, with different coefficients)
	// this is probably completely un-needed but gives the highest quality results,
	G4double alpha=(1.0+ lambda0)/30.0;
	G4double xvals[]={0.98302349, 0.84652241, 0.53235309, 0.18347974};
	G4double weights[]={0.03472124, 0.14769029, 0.23485003, 0.18602489};
	for(G4int k=0; k<4; k++) {
		G4double x, ff;
		x=xx1/xvals[k];
		ff=1.0/std::sqrt(1.0-phiData(x*au)/(x*eps)-beta*beta/(x*x));
		alpha+=weights[k]*ff;
	}
	
	G4double thetac1=CLHEP::pi*beta*alpha/xx1; // complement of CM scattering angle
	G4double sintheta=std::sin(thetac1); //note sin(pi-theta)=sin(theta)
	G4double costheta=-std::cos(thetac1); // note cos(pi-theta)=-cos(theta)
										  // G4double psi=std::atan2(sintheta, costheta+a1/A); // lab scattering angle (M&T 3rd eq. 8.69)
	
	// numerics note:  because we checked above for reasonable values of beta which give real recoils,
	// we don't have to look too closely for theta -> 0 here (which would cause sin(theta)
	// and 1-cos(theta) to both vanish and make the atan2 ill behaved).
	G4double zeta=std::atan2(sintheta, 1-costheta); // lab recoil angle (M&T 3rd eq. 8.73)
	G4double coszeta=std::cos(zeta);
	G4double sinzeta=std::sin(zeta);

	kin.sinTheta=sintheta;
	kin.cosTheta=costheta;
	kin.sinZeta=sinzeta;
	kin.cosZeta=coszeta;
	return 1; // all OK, collision is valid
}

void G4ScreenedCoulombClassicalKinematics::DoCollisionStep(G4ScreenedNuclearRecoil *master,
		const G4Track& aTrack, const G4Step&) {
	
	if(!master->GetValidCollision()) return;
	
	G4ParticleChange &aParticleChange=master->GetParticleChange();
	G4CoulombKinematicsInfo &kin=master->GetKinematics();
	
	const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
	G4ParticleDefinition *baseParticle=aTrack.GetDefinition();
	
	G4double incidentEnergy = incidentParticle->GetKineticEnergy();
	
	// this adjustment to a1 gives the right results for soft (constant gamma)
	// relativistic collisions.  Hard collisions are wrong anyway, since the
	// Coulombic and hadronic terms interfere and cannot be added.
	G4double gamma=(1.0+incidentEnergy/baseParticle->GetPDGMass());
	G4double a1=kin.a1*gamma; // relativistic gamma correction
	
	G4ParticleDefinition *recoilIon=kin.recoilIon;
	G4double A=recoilIon->GetPDGMass()/amu_c2;
	G4int Z=(G4int)((recoilIon->GetPDGCharge()/eplus)+0.5);
	
	G4double Ec = incidentEnergy*(A/(A+a1)); // energy in CM frame (non-relativistic!)
	const G4ScreeningTables *screen=kin.crossSection->GetScreening(Z);
	G4double au=screen->au; // screening length
	
	G4double beta = kin.impactParameter/au; // dimensionless impact parameter
	G4double eps = Ec/(screen->z1*Z*elm_coupling/au); // dimensionless energy
	
	G4bool ok=DoScreeningComputation(master, screen, eps, beta);	
	if(!ok) {
		master->SetValidCollision(0); // flag bad collision
		return; // just bail out without setting valid flag
	}
	
	G4double eRecoil=4*incidentEnergy*a1*A*kin.cosZeta*kin.cosZeta/((a1+A)*(a1+A));
	kin.eRecoil=eRecoil;
	
	if(incidentEnergy-eRecoil < master->GetRecoilCutoff()) {
		if(!baseParticle->GetProcessManager()->
		   GetAtRestProcessVector()->size())
			aParticleChange.ProposeTrackStatus(fStopAndKill);
		else
			aParticleChange.ProposeTrackStatus(fStopButAlive);
		aParticleChange.ProposeEnergy(0.0);
		master->AddToNIEL(incidentEnergy-eRecoil);
	} 
	
	if(master->GetEnableRecoils() && eRecoil > master->GetRecoilCutoff()) {
		kin.recoilIon=recoilIon;		
	} else {
		kin.recoilIon=0; // this flags no recoil to be generated
		master->AddToNIEL(eRecoil) ;
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
	
	G4double azimuth=G4UniformRand()*(2.0*CLHEP::pi);
	G4double sa=std::sin(azimuth);
	G4double ca=std::cos(azimuth);

	G4ThreeVector recoilMomentumDirection(kin.sinZeta*ca, kin.sinZeta*sa, kin.cosZeta);
	G4ParticleMomentum incidentDirection = incidentParticle->GetMomentumDirection();
	recoilMomentumDirection=recoilMomentumDirection.rotateUz(incidentDirection);
	G4ThreeVector recoilMomentum=recoilMomentumDirection*std::sqrt(2.0*eRecoil*kin.a2*amu_c2);
	
	if(aParticleChange.GetEnergy() != 0.0) { // DoKinematics hasn't stopped it!
		G4ThreeVector beamMomentum=incidentParticle->GetMomentum()-recoilMomentum;
		aParticleChange.ProposeMomentumDirection(beamMomentum.unit()) ;
		aParticleChange.ProposeEnergy(incidentEnergy-eRecoil);
	}
	
	if(kin.recoilIon) {
		G4DynamicParticle* recoil = new G4DynamicParticle (kin.recoilIon,
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
DumpPhysicsTable(const G4ParticleDefinition&)
{
}

// This used to be the file mhmScreenedNuclearRecoil_native.cc
// it has been included here to collect this file into a smaller number of packages

#include "G4DataVector.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ElementVector.hh"
#include <vector>

static c2_function<G4double> &ZBLScreening(G4int z1, G4int z2, size_t npoints, G4double rMax, G4double *auval)
{
	static const size_t ncoef=4;
	static G4double scales[ncoef]={-3.2, -0.9432, -0.4028, -0.2016};
	static G4double coefs[ncoef]={0.1818,0.5099,0.2802,0.0281};
	
	G4double au=0.8854*angstrom*0.529/(std::pow(z1, 0.23)+std::pow(z2,0.23));
	std::vector<G4double> r(npoints), phi(npoints);
	
	for(size_t i=0; i<npoints; i++) {
		G4double rr=(float)i/(float)(npoints-1);
		r[i]=rr*rr*rMax; // use quadratic r scale to make sampling fine near the center
		G4double sum=0.0;
		for(size_t j=0; j<ncoef; j++) sum+=coefs[j]*std::exp(scales[j]*r[i]/au);
		phi[i]=sum;
	}

	// compute the derivative at the origin for the spline
	G4double phiprime0=0.0;
	for(size_t j=0; j<ncoef; j++) phiprime0+=scales[j]*coefs[j]*std::exp(scales[j]*r[0]/au);
	phiprime0*=(1.0/au); // put back in natural units;
	
	*auval=au;
	return *static_cast<c2_function<G4double> *>(new lin_log_interpolating_function<G4double>(r, phi, false, phiprime0));
}

static c2_function<G4double> &MoliereScreening(G4int z1, G4int z2, size_t npoints, G4double rMax, G4double *auval)
{	
	static const size_t ncoef=3;
	static G4double scales[ncoef]={-6.0, -1.2, -0.3};
	static G4double coefs[ncoef]={0.10, 0.55, 0.35};
	
	G4double au=0.8853*0.529*angstrom/std::sqrt(std::pow(z1, 0.6667)+std::pow(z2,0.6667));
	std::vector<G4double> r(npoints), phi(npoints);
	
	for(size_t i=0; i<npoints; i++) {
		G4double rr=(float)i/(float)(npoints-1);
		r[i]=rr*rr*rMax; // use quadratic r scale to make sampling fine near the center
		G4double sum=0.0;
		for(size_t j=0; j<ncoef; j++) sum+=coefs[j]*std::exp(scales[j]*r[i]/au);
		phi[i]=sum;
	}
	
	// compute the derivative at the origin for the spline
	G4double phiprime0=0.0;
	for(size_t j=0; j<ncoef; j++) phiprime0+=scales[j]*coefs[j]*std::exp(scales[j]*r[0]/au);
	phiprime0*=(1.0/au); // put back in natural units;
	
	*auval=au;
	return *static_cast<c2_function<G4double> *>(new lin_log_interpolating_function<G4double>(r, phi, false, phiprime0));
}

static c2_function<G4double> &LJScreening(G4int z1, G4int z2, size_t npoints, G4double rMax, G4double *auval)
{	
//from Loftager, Besenbacher, Jensen & Sorensen
//PhysRev A20, 1443++, 1979
	G4double au=0.8853*0.529*angstrom/std::sqrt(std::pow(z1, 0.6667)+std::pow(z2,0.6667));
	std::vector<G4double> r(npoints), phi(npoints);
	
	for(size_t i=0; i<npoints; i++) {
		G4double rr=(float)i/(float)(npoints-1);
		r[i]=rr*rr*rMax; // use quadratic r scale to make sampling fine near the center

		G4double y=std::sqrt(9.67*r[i]/au);
		G4double ysq=y*y;
		G4double phipoly=1+y+0.3344*ysq+0.0485*y*ysq+0.002647*ysq*ysq;
		phi[i]=phipoly*std::exp(-y);
		// G4cout << r[i] << " " << phi[i] << G4endl;
	}

	// compute the derivative at the origin for the spline
	G4double logphiprime0=(9.67/2.0)*(2*0.3344-1.0); // #avoid 0/0 on first element
	logphiprime0 *= (1.0/au); // #put back in natural units
	
	*auval=au;
	return *static_cast<c2_function<G4double> *>(new lin_log_interpolating_function<G4double>(r, phi, false, logphiprime0*phi[0]));
}

G4NativeScreenedCoulombCrossSection::~G4NativeScreenedCoulombCrossSection() {
}

G4NativeScreenedCoulombCrossSection::G4NativeScreenedCoulombCrossSection() {
	AddScreeningFunction("zbl", ZBLScreening);
	AddScreeningFunction("lj", LJScreening);	
	AddScreeningFunction("mol", MoliereScreening);	
}

std::vector<G4String> G4NativeScreenedCoulombCrossSection::GetScreeningKeys() const {
	std::vector<G4String> keys;
	// find the available screening keys
	std::map<std::string, ScreeningFunc>::const_iterator sfunciter=phiMap.begin();
	for(; sfunciter != phiMap.end(); sfunciter++) keys.push_back((*sfunciter).first);
	return keys;
}

static inline G4double cm_energy(G4double a1, G4double a2, G4double t0) {
	// "relativistically correct energy in CM frame"
	G4double m1=a1*amu_c2, m2=a2*amu_c2;
	G4double mc2=(m1+m2);
	G4double f=2.0*m2*t0/(mc2*mc2);
	// old way: return (f < 1e-6) ?  0.5*mc2*f : mc2*(std::sqrt(1.0+f)-1.0);
	// formally equivalent to previous, but numerically stable for all f without conditional
	// uses identity (sqrt(1+x) - 1)(sqrt(1+x) + 1) = x
	return mc2*f/(std::sqrt(1.0+f)+1.0); 
}

static inline G4double  thetac(G4double m1, G4double m2, G4double eratio) {
	G4double s2th2=eratio*( (m1+m2)*(m1+m2)/(4.0*m1*m2) );
	G4double sth2=std::sqrt(s2th2);
	return 2.0*std::asin(sth2);
}

void G4NativeScreenedCoulombCrossSection::LoadData(G4String screeningKey, G4int z1, G4double a1, G4double recoilCutoff)
{			
	static const size_t sigLen=200; // since sigma doesn't matter much, a very coarse table will do			
	G4DataVector energies(sigLen);
	G4DataVector data(sigLen);
	
	a1=standardmass(z1); // use standardized values for mass for building tables
	
	const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
	if (materialTable == 0)
		G4Exception("mhmNativeCrossSection::LoadData - no MaterialTable found)");
	
	G4int nMaterials = G4Material::GetNumberOfMaterials();
	
	for (G4int m=0; m<nMaterials; m++)
    {
		const G4Material* material= (*materialTable)[m];
		const G4ElementVector* elementVector = material->GetElementVector();
		const G4int nMatElements = material->GetNumberOfElements();
		
		for (G4int iEl=0; iEl<nMatElements; iEl++)
		{
			G4Element* element = (*elementVector)[iEl];
			G4int Z = (G4int) element->GetZ();
			G4double a2=element->GetA()*(mole/gram);
			
			if(sigmaMap.find(Z)!=sigmaMap.end()) continue; // we've already got this element
			
			// find the screening function generator we need
			std::map<std::string, ScreeningFunc>::iterator sfunciter=phiMap.find(screeningKey);
			if(sfunciter==phiMap.end()) {
				G4cout << "no such screening key " << screeningKey << G4endl; // FIXME later
				exit(1);
			}
			ScreeningFunc sfunc=(*sfunciter).second;
			
			G4double au; 
			c2_function<G4double> &screen=sfunc(z1, Z, 200, 50.0*angstrom, &au); // generate the screening data

			G4ScreeningTables st;
			st.EMphiData=&screen; // this is our phi table
			st.z1=z1; st.m1=a1; st.z2=Z; st.m2=a2; st.emin=recoilCutoff;
			st.au=au;

			// now comes the hard part... build the total cross section tables from the phi table			
			//based on (pi-thetac) = pi*beta*alpha/x0, but noting that alpha is very nearly unity, always
			//so just solve it wth alpha=1, which makes the solution much easier
			//this function returns an approximation to (beta/x0)^2=phi(x0)/(eps*x0)-1 ~ ((pi-thetac)/pi)^2
			//Since we don't need exact sigma values, this is good enough (within a factor of 2 almost always)
			//this rearranges to phi(x0)/(x0*eps) = 2*theta/pi - theta^2/pi^2
			
			c2_linear<G4double> c2au(0.0, 0.0, au);
			c2_composed_function<G4double> phiau(screen, c2au); // build phi(x*au) for dimensionless phi
			c2_linear<G4double> c2eps(0.0, 0.0, 0.0); // will store an appropriate eps inside this in loop
			c2_ratio<G4double> x0func(phiau, c2eps); // this will be phi(x)/(x*eps) when c2eps is correctly set
			
			G4double m1c2=a1*amu_c2;
			G4double escale=z1*Z*elm_coupling/au; // energy at screening distance
			G4double emax=m1c2; // model is doubtful in very relativistic range
			G4double eratkin=0.9999*(4*a1*a2)/((a1+a2)*(a1+a2)); // #maximum kinematic ratio possible at 180 degrees
			G4double cmfact0=st.emin/cm_energy(a1, a2, st.emin);
			G4double l1=std::log(emax);
			G4double l0=std::log(st.emin*cmfact0/eratkin);

			if(verbosity >=1) 
				G4cout << "Native Screening: " << screeningKey << " " << z1 << " " << a1 << " " << 
					Z << " " << a2 << " " << recoilCutoff << G4endl;
			
			for(size_t idx=0; idx<sigLen; idx++) {
				G4double ee=std::exp(idx*((l1-l0)/sigLen)+l0);
				G4double gamma=1.0+ee/m1c2;
				G4double eratio=(cmfact0*st.emin)/ee; // factor by which ee needs to be reduced to get emin
				G4double theta=thetac(gamma*a1, a2, eratio);
			
				G4double eps=cm_energy(a1, a2, ee)/escale; // #make sure lab energy is converted to CM for these calculations
				c2eps.reset(0.0, 0.0, eps); // set correct slope in this function
				
				G4double q=theta/pi;
				// G4cout << ee << " " << m1c2 << " " << gamma << " " << eps << " " << theta << " " << q << G4endl;
				G4double x0= x0func.find_root(1e-6*angstrom/au, 0.9999*screen.xmax()/au, 1.0, 2*q-q*q);
				G4double betasquared=x0*x0 - x0*phiau(x0)/eps;	
				G4double sigma=pi*betasquared*au*au;
				energies[idx]=ee;
				data[idx]=sigma;
			}
			
			screeningData[Z]=st;
			sigmaMap[Z] = static_cast<c2_function<G4double> *>(new log_log_interpolating_function<G4double>(
				energies, data));
		}
	}
}
