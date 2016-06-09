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
#include "G4AdjointBremsstrahlungModel.hh"
#include "G4AdjointCSManager.hh"
#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include "G4Timer.hh"

////////////////////////////////////////////////////////////////////////////////
//
G4AdjointBremsstrahlungModel::G4AdjointBremsstrahlungModel():
 G4VEmAdjointModel("AdjointBremModel"),
 probsup(1.0),
 MigdalConstant(classic_electr_radius*electron_Compton_length*electron_Compton_length/pi),
 LPMconstant(fine_structure_const*electron_mass_c2*electron_mass_c2/(4.*pi*hbarc)),
 theLPMflag(true)

{ isElectron= true;
  SetUseMatrix(true);
  SetUseMatrixPerElement(false);
  SetApplyCutInRange(true);
  SetIsIonisation(false);
  highKinEnergy= 100.*TeV;
  lowKinEnergy = 1.0*keV;
  theTimer =new G4Timer();
  
  theTimer->Start();
  InitialiseParameters();
  theTimer->Stop();
  G4cout<<"Time elapsed in second for the initialidation of AdjointBrem "<<theTimer->GetRealElapsed()<<std::endl;
  
  ModeldCS="MODEL1";
  
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointBremsstrahlungModel::~G4AdjointBremsstrahlungModel()
{;}
////////////////////////////////////////////////////////////////////////////////
//
/*G4double G4AdjointBremsstrahlungModel::DiffCrossSectionPerVolumePrimToSecond(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      )
{
 
 static const G4double
     ah10 = 4.67733E+00, ah11 =-6.19012E-01, ah12 = 2.02225E-02,
     ah20 =-7.34101E+00, ah21 = 1.00462E+00, ah22 =-3.20985E-02,
     ah30 = 2.93119E+00, ah31 =-4.03761E-01, ah32 = 1.25153E-02;

  static const G4double
     bh10 = 4.23071E+00, bh11 =-6.10995E-01, bh12 = 1.95531E-02,
     bh20 =-7.12527E+00, bh21 = 9.69160E-01, bh22 =-2.74255E-02,
     bh30 = 2.69925E+00, bh31 =-3.63283E-01, bh32 = 9.55316E-03;

  static const G4double
     al00 =-2.05398E+00, al01 = 2.38815E-02, al02 = 5.25483E-04,
     al10 =-7.69748E-02, al11 =-6.91499E-02, al12 = 2.22453E-03,
     al20 = 4.06463E-02, al21 =-1.01281E-02, al22 = 3.40919E-04;

  static const G4double
     bl00 = 1.04133E+00, bl01 =-9.43291E-03, bl02 =-4.54758E-04,
     bl10 = 1.19253E-01, bl11 = 4.07467E-02, bl12 =-1.30718E-03,
     bl20 =-1.59391E-02, bl21 = 7.27752E-03, bl22 =-1.94405E-04;

  static const G4double tlow = 1.*MeV;
  
  G4double dCrossEprod=0.;
  G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
  G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
 
 
 if (kinEnergyProj>Emin_proj && kinEnergyProj<=Emax_proj){
 	
  	G4double cross = 0.0;
	

	
	G4double E1=kinEnergyProd;
 	G4double E2=kinEnergyProd*1.000000001;
 	G4double dE=(E2-E1);

  	const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  	const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();
  	G4double dum=0.;
  
  	for (size_t i=0; i<aMaterial->GetNumberOfElements(); i++) {

    		G4double fac=
		
		cross += theAtomNumDensityVector[i] * theDirectEMModel->ComputeCrossSectionPerAtom(G4Electron::Electron(),
                      kinEnergyProj, (*theElementVector)[i]->GetZ(), dum,E1);
    	
      		
    		
  	} 
 	dCrossEprod=(cross1-cross2)/dE; //first term
	
	//Now come the correction
	//-----------------------
	
	//First compute fsig for E1
	//-------------------------
	
	
	G4double totalEnergy = kinEnergyProj+electron_mass_c2 ;
  	G4double kp2 = MigdalConstant*totalEnergy*totalEnergy
                                             *(aMaterial->GetElectronDensity());

  	G4double fsig = 0.;
  	G4int nmax = 100;
  	G4double vmin=std::log(E1);
  	G4double vmax=std::log(kinEnergyProj) ;
  	G4int nn = (G4int)(nmax*(vmax-vmin)/(std::log(highKinEnergy)-vmin));
  	G4double u,fac,c,v,dv,y ;
  	if(nn > 0) {

      		dv = (vmax-vmin)/nn ;
      		v  = vmin-dv ;
      		for(G4int n=0; n<=nn; n++) {

        		v += dv;  
        		u = std::exp(v);              
        		fac = SupressionFunction(aMaterial, kinEnergyProj, u);
        		y = u/kinEnergyProj;
        		fac *= (4.-4.*y+3.*y*y)/3.;
        		fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;

        		if ((n==0)||(n==nn)) c=0.5;
        		else    c=1. ;

        		fac  *= c;
        		fsig += fac;
      		}
      		y = E1/kinEnergyProj ;
      		fsig *=dv/(-4.*std::log(y)/3.-4.*(1.-y)/3.+0.5*(1.-y*y));

  	} 
	else {
		fsig = 1.;
  	}
  	if (fsig > 1.) fsig = 1.;
	
	dCrossEprod*=fsig;
	//return dCrossEprod;
	//Now we  compute dfsig 
	//-------------------------
	G4double dfsig = 0.;
  	nn=20;
	vmax=std::log(E2) ;
	dv = (vmax-vmin)/nn ;
      	v  = vmin-dv ;
      	for(G4int n=0; n<=nn; n++) {
		v += dv;  
        	u = std::exp(v);              
        	fac = SupressionFunction(aMaterial, kinEnergyProj, u);
        	y = u/kinEnergyProj;
        	fac *= (4.-4.*y+3.*y*y)/3.;
        	fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;

        	if ((n==0)||(n==nn)) c=0.5;
        	else    c=1. ;

        	fac  *= c;
        	dfsig += fac;
      	}
      	y = E1/kinEnergyProj;
      	dfsig *=dv/(-4.*std::log(y)/3.-4.*(1.-y)/3.+0.5*(1.-y*y));
	
	dCrossEprod+=dfsig*cross1/dE;
	
	
	
	 
	
 }
 return dCrossEprod;
  
} 
*/
G4double G4AdjointBremsstrahlungModel::DiffCrossSectionPerVolumePrimToSecond(const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      )
{if (ModeldCS=="MODEL2") return  DiffCrossSectionPerVolumePrimToSecond2(aMaterial,
                                      			         kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      			         kinEnergyProd);
 if (ModeldCS=="MODEL3") return  DiffCrossSectionPerVolumePrimToSecond3(aMaterial,
                                      			         kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      			         kinEnergyProd);
 return  DiffCrossSectionPerVolumePrimToSecond1(aMaterial,
                                      			         kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      			         kinEnergyProd);								 
}				      
////////////////////////////////////////////////////////////////////////////////
// the one used till now
G4double G4AdjointBremsstrahlungModel::DiffCrossSectionPerVolumePrimToSecond1(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      )
{
 G4double dCrossEprod=0.;
 G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
 G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
 
 
 if (kinEnergyProj>Emin_proj && kinEnergyProj<=Emax_proj){
 	
  	G4double cross1 = 0.0;
	G4double cross2 = 0.0;

	
	G4double E1=kinEnergyProd;
 	G4double E2=kinEnergyProd*1.01;
 	G4double dE=(E2-E1);

  	const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  	const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();
  	G4double dum=0.;
  
  	for (size_t i=0; i<aMaterial->GetNumberOfElements(); i++) {

    		cross1 += theAtomNumDensityVector[i] * theDirectEMModel->ComputeCrossSectionPerAtom(G4Electron::Electron(),
                      kinEnergyProj, (*theElementVector)[i]->GetZ(), dum,E1);
    	
      		cross2 += theAtomNumDensityVector[i] * theDirectEMModel->ComputeCrossSectionPerAtom(G4Electron::Electron(),
                     kinEnergyProj, (*theElementVector)[i]->GetZ(), dum, E2);
    		
  	} 
 	dCrossEprod=(cross1-cross2)/dE; //first term
	
	//Now come the correction
	//-----------------------
	
	//First compute fsig for E1
	//-------------------------
	
	
	G4double totalEnergy = kinEnergyProj+electron_mass_c2 ;
  	G4double kp2 = MigdalConstant*totalEnergy*totalEnergy
                                             *(aMaterial->GetElectronDensity());

  	G4double fsig1 = 0.;
  	G4int nmax = 100;
  	G4double vmin=std::log(E1);
  	G4double vmax=std::log(kinEnergyProj) ;
  	G4int nn = (G4int)(nmax*(vmax-vmin)/(std::log(highKinEnergy)-vmin));
  	G4double u,fac,c,v,dv,y ;
  	if(nn > 0) {

      		dv = (vmax-vmin)/nn ;
      		v  = vmin-dv ;
      		for(G4int n=0; n<=nn; n++) {

        		v += dv;  
        		u = std::exp(v);              
        		fac = SupressionFunction(aMaterial, kinEnergyProj, u);
        		y = u/kinEnergyProj;
        		fac *= (4.-4.*y+3.*y*y)/3.;
        		fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;

        		if ((n==0)||(n==nn)) c=0.5;
        		else    c=1. ;

        		fac  *= c;
        		fsig1 += fac;
      		}
      		y = E1/kinEnergyProj ;
      		fsig1 *=dv/(-4.*std::log(y)/3.-4.*(1.-y)/3.+0.5*(1.-y*y));

  	} 
	else {
		fsig1 = 1.;
  	}
  	if (fsig1 > 1.) fsig1 = 1.;
	
	dCrossEprod*=fsig1;
	
	
	G4double fsig2 = 0.;
  	vmin=std::log(E2);
	nn = (G4int)(nmax*(vmax-vmin)/(std::log(highKinEnergy)-vmin));
  	if(nn > 0) {

      		dv = (vmax-vmin)/nn ;
      		v  = vmin-dv ;
      		for(G4int n=0; n<=nn; n++) {

        		v += dv;  
        		u = std::exp(v);              
        		fac = SupressionFunction(aMaterial, kinEnergyProj, u);
        		y = u/kinEnergyProj;
        		fac *= (4.-4.*y+3.*y*y)/3.;
        		fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;

        		if ((n==0)||(n==nn)) c=0.5;
        		else    c=1. ;

        		fac  *= c;
        		fsig2 += fac;
      		}
      		y = E2/kinEnergyProj ;
      		fsig2 *=dv/(-4.*std::log(y)/3.-4.*(1.-y)/3.+0.5*(1.-y*y));

  	} 
	else {
		fsig2 = 1.;
  	}
  	if (fsig2 > 1.) fsig2 = 1.;
	

	G4double dfsig=(fsig2-fsig1);
	dCrossEprod+=dfsig*cross1/dE;
	
	dCrossEprod=(fsig1*cross1-fsig2*cross2)/dE;
	
	
	
	
	
	/*if (fsig < 1.){
		//Now we  compute dfsig 
		//-------------------------
		G4double dfsig = 0.;
  		nn=20;
		vmax=std::log(E2) ;
		dv = (vmax-vmin)/nn ;
      		v  = vmin-dv ;
      		for(G4int n=0; n<=nn; n++) {
			v += dv;  
        		u = std::exp(v);              
        		fac = SupressionFunction(aMaterial, kinEnergyProj, u);
        		y = u/kinEnergyProj;
        		fac *= (4.-4.*y+3.*y*y)/3.;
        		fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;

        		if ((n==0)||(n==nn)) c=0.5;
        		else    c=1. ;

        		fac  *= c;
        		dfsig += fac;
      		}
      		y = E1/kinEnergyProj;
      		dfsig *=dv/(-4.*std::log(y)/3.-4.*(1.-y)/3.+0.5*(1.-y*y));
		dCrossEprod+=dfsig*cross1/dE;
		
	}	
	*/
	
	
	
	
	
	 
	
 }
 return dCrossEprod;
  
} 


////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointBremsstrahlungModel::DiffCrossSectionPerVolumePrimToSecond2(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      )
{
 G4double dCrossEprod=0.;
 G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
 G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
 
 
 if (kinEnergyProj>Emin_proj && kinEnergyProj<=Emax_proj){
 	
  	G4double dEdX1 = 0.0;
	G4double dEdX2 = 0.0;

	
	G4double E1=kinEnergyProd;
 	G4double E2=kinEnergyProd*1.001;
 	G4double dE=(E2-E1);
  	//G4double dum=0.;
	
	dEdX1 = theDirectEMModel->ComputeDEDXPerVolume(aMaterial,G4Electron::Electron(),kinEnergyProj,E1); 
	dEdX2 = theDirectEMModel->ComputeDEDXPerVolume(aMaterial,G4Electron::Electron(),kinEnergyProj,E2);
	dCrossEprod=(dEdX2-dEdX1)/dE/E1;
	
	 
	
	
	
	 
	
 }
 return dCrossEprod;
  
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointBremsstrahlungModel::DiffCrossSectionPerVolumePrimToSecond3(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      )
{
 
 return G4VEmAdjointModel::DiffCrossSectionPerVolumePrimToSecond(aMaterial,
                                      			         kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      			         kinEnergyProd);
  
}

////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointBremsstrahlungModel::SupressionFunction(const G4Material* material,
                                 G4double kineticEnergy, G4double gammaEnergy)
{
  // supression due to the LPM effect+polarisation of the medium/
  // supression due to the polarisation alone


  G4double totEnergy = kineticEnergy+electron_mass_c2 ;
  G4double totEnergySquare = totEnergy*totEnergy ;

  G4double LPMEnergy = LPMconstant*(material->GetRadlen()) ;

  G4double gammaEnergySquare = gammaEnergy*gammaEnergy ;

  G4double electronDensity = material->GetElectronDensity();

  G4double sp = gammaEnergySquare/
   (gammaEnergySquare+MigdalConstant*totEnergySquare*electronDensity);

  G4double supr = 1.0;

  if (theLPMflag) {

    G4double s2lpm = LPMEnergy*gammaEnergy/totEnergySquare;

    if (s2lpm < 1.) {

      G4double LPMgEnergyLimit = totEnergySquare/LPMEnergy ;
      G4double LPMgEnergyLimit2 = LPMgEnergyLimit*LPMgEnergyLimit;
      G4double splim = LPMgEnergyLimit2/
        (LPMgEnergyLimit2+MigdalConstant*totEnergySquare*electronDensity);
      G4double w = 1.+1./splim ;

      if ((1.-sp) < 1.e-6) w = s2lpm*(3.-sp);
      else                 w = s2lpm*(1.+1./sp);

      supr = (std::sqrt(w*w+4.*s2lpm)-w)/(std::sqrt(w*w+4.)-w) ;
      supr /= sp;    
    } 
    
  } 
  return supr;
}

////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointBremsstrahlungModel::SampleSecondaries(const G4Track& aTrack,
                       G4bool IsScatProjToProjCase,
	               G4ParticleChange* fParticleChange)
{ 

  //G4cout<<"Adjoint Brem"<<std::endl;
  const G4DynamicParticle* theAdjointPrimary =aTrack.GetDynamicParticle();
  
  size_t ind=0;
  
  if (UseMatrixPerElement ) { //Select Material
   	std::vector<double>* CS_Vs_Element = &CS_Vs_ElementForScatProjToProjCase;
  	if ( !IsScatProjToProjCase) CS_Vs_Element = &CS_Vs_ElementForProdToProjCase;
  	G4double rand_var= G4UniformRand();
  	G4double SumCS=0.;
  	for (size_t i=0;i<CS_Vs_Element->size();i++){
 		SumCS+=(*CS_Vs_Element)[i];
		if (rand_var<=SumCS/lastCS){
			ind=i;
			break;
		}
  	}
  }
  else 	{
  	ind = currentMaterialIndex;
  }
 
 
 //Elastic inverse scattering modified compared to general G4VEmAdjointModel
 //---------------------------
 G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
 G4double adjointPrimTotalEnergy = theAdjointPrimary->GetTotalEnergy();
 //G4double adjointPrimP =theAdjointPrimary->GetTotalMomentum();
 if (adjointPrimKinEnergy>HighEnergyLimit*0.999){
 	return;
 }
 
 //Sample secondary energy
 //-----------------------
 
 G4double projectileKinEnergy = SampleAdjSecEnergyFromCSMatrix(ind,
 						   adjointPrimKinEnergy,
						   IsScatProjToProjCase);
				   
 
 
 
 //Weight correction
 //-----------------------					   
 CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(), adjointPrimKinEnergy,projectileKinEnergy);	
 
 
 //Kinematic
 //---------
 
 G4double projectileM0 = electron_mass_c2;
 G4double projectileTotalEnergy = projectileM0+projectileKinEnergy;
 G4double projectileP2 = projectileTotalEnergy*projectileTotalEnergy - projectileM0*projectileM0;	
 G4double projectileP = std::sqrt(projectileP2);
 
 
 //Angle of the gamma direction with the projectile taken from G4eBremsstrahlungModel
 //------------------------------------------------
  G4double u;
  const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  if (9./(9.+d) > G4UniformRand()) u = - std::log(G4UniformRand()*G4UniformRand())/a1;
     else                          u = - std::log(G4UniformRand()*G4UniformRand())/a2;

  G4double theta = u*electron_mass_c2/projectileTotalEnergy;

  G4double sint = std::sin(theta);
  G4double cost = std::cos(theta);

  G4double phi = twopi * G4UniformRand() ;
  
  G4ThreeVector projectileMomentum;
  projectileMomentum=G4ThreeVector(std::cos(phi)*sint,std::sin(phi)*sint,cost)*projectileP; //gamma frame
  if (IsScatProjToProjCase) {//the adjoint primary is the scattered e-
  	G4ThreeVector gammaMomentum = (projectileTotalEnergy-adjointPrimTotalEnergy)*G4ThreeVector(0.,0.,1.);
	G4ThreeVector dirProd=projectileMomentum-gammaMomentum;
	G4double cost1 = std::cos(dirProd.angle(projectileMomentum));
	G4double sint1 =  std::sqrt(1.-cost1*cost1);
	projectileMomentum=G4ThreeVector(std::cos(phi)*sint1,std::sin(phi)*sint1,cost1)*projectileP;
  
  }
  
  projectileMomentum.rotateUz(theAdjointPrimary->GetMomentumDirection());
 
 
 
  if (!IsScatProjToProjCase && CorrectWeightMode){ //kill the primary and add a secondary
 	fParticleChange->ProposeTrackStatus(fStopAndKill);
 	fParticleChange->AddSecondary(new G4DynamicParticle(theAdjEquivOfDirectPrimPartDef,projectileMomentum));
	//G4cout<<"projectileMomentum "<<projectileMomentum<<std::endl;
  }
  else {
 	fParticleChange->ProposeEnergy(projectileKinEnergy);
	fParticleChange->ProposeMomentumDirection(projectileMomentum.unit());
	//G4cout<<"projectileMomentum "<<projectileMomentum<<std::endl;
  }	
} 
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointBremsstrahlungModel::DefineDirectBremModel(G4eBremsstrahlungModel* aModel)
{theDirectBremModel=aModel;
 DefineDirectEMModel(aModel);
} 
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointBremsstrahlungModel::InitialiseParameters()
{  
  static const G4double
     ah10 = 4.67733E+00, ah11 =-6.19012E-01, ah12 = 2.02225E-02,
     ah20 =-7.34101E+00, ah21 = 1.00462E+00, ah22 =-3.20985E-02,
     ah30 = 2.93119E+00, ah31 =-4.03761E-01, ah32 = 1.25153E-02;

  static const G4double
     bh10 = 4.23071E+00, bh11 =-6.10995E-01, bh12 = 1.95531E-02,
     bh20 =-7.12527E+00, bh21 = 9.69160E-01, bh22 =-2.74255E-02,
     bh30 = 2.69925E+00, bh31 =-3.63283E-01, bh32 = 9.55316E-03;

 /* static const G4double
     al00 =-2.05398E+00, al01 = 2.38815E-02, al02 = 5.25483E-04,
     al10 =-7.69748E-02, al11 =-6.91499E-02, al12 = 2.22453E-03,
     al20 = 4.06463E-02, al21 =-1.01281E-02, al22 = 3.40919E-04;

  static const G4double
     bl00 = 1.04133E+00, bl01 =-9.43291E-03, bl02 =-4.54758E-04,
     bl10 = 1.19253E-01, bl11 = 4.07467E-02, bl12 =-1.30718E-03,
     bl20 =-1.59391E-02, bl21 = 7.27752E-03, bl22 =-1.94405E-04;*/
 
  
  const G4ElementTable* theElementTable = G4Element::GetElementTable();
  FZ.clear();
  ah1.clear();
  ah2.clear();
  ah3.clear();
  
  bh1.clear();
  bh2.clear();
  bh3.clear();
  
  al0.clear();
  al1.clear();
  al2.clear();
  
  bl0.clear();
  bl1.clear();
  bl2.clear();
  SigmaPerAtom.clear(); 
  
  for (size_t j=0; j<theElementTable->size();j++){
	
	G4Element* anElement=(*theElementTable)[j]; 
	G4double lnZ = 3.*(anElement->GetIonisation()->GetlogZ3());
  	FZ.push_back(lnZ* (4.- 0.55*lnZ));
  	G4double ZZ = anElement->GetIonisation()->GetZZ3();
	
	ah1.push_back(ah10 + ZZ* (ah11 + ZZ* ah12));
        ah2.push_back(ah20 + ZZ* (ah21 + ZZ* ah22));
        ah3.push_back(ah30 + ZZ* (ah31 + ZZ* ah32));

        bh1.push_back(bh10 + ZZ* (bh11 + ZZ* bh12));
        bh2.push_back(bh20 + ZZ* (bh21 + ZZ* bh22));
        bh3.push_back(bh30 + ZZ* (bh31 + ZZ* bh32));
	/*SigmaPerAtom.push_back(theDirectEMModel->ComputeCrossSectionPerAtom(
					theDirectPrimaryPartDef,GetHighEnergyLimit()/2., 
					anElement->GetZ(),1.,GetLowEnergyLimit(),1.e20));*/
	
	
  	
  }	
}
