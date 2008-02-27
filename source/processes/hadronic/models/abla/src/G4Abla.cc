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
// $Id: G4Abla.cc,v 1.1 2008-02-27 18:31:11 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include <time.h>

#include "G4Abla.hh"
#include "G4InclAblaDataFile.hh"
#include "Randomize.hh"
#include <assert.h>

G4Abla::G4Abla()
{
  ilast = 0;
}

G4Abla::G4Abla(G4Hazard *hazard, G4Volant *volant)
{
  verboseLevel = 0;
  ilast = 0;
  volant = volant; // ABLA internal particle data
  volant->iv = 0;
  hazard = hazard; // Random seeds

  varntp = new G4VarNtp();
  pace = new G4Pace();
  ald = new G4Ald();
  ablamain = new G4Ablamain();
  emdpar = new G4Emdpar();
  eenuc = new G4Eenuc();
  ec2sub = new G4Ec2sub();
  ecld = new G4Ecld();
  fb = new G4Fb();
  fiss = new G4Fiss();
  opt = new G4Opt();
}

G4Abla::G4Abla(G4Hazard *aHazard, G4Volant *aVolant, G4VarNtp *aVarntp)
{
  verboseLevel = 0;
  ilast = 0;
  volant = aVolant; // ABLA internal particle data
  volant->iv = 0;
  hazard = aHazard; // Random seeds
  varntp = aVarntp; // Output data structure
  varntp->ntrack = 0;
 
  pace = new G4Pace();
  ald = new G4Ald();
  ablamain = new G4Ablamain();
  emdpar = new G4Emdpar();
  eenuc = new G4Eenuc();
  ec2sub = new G4Ec2sub();
  ecld = new G4Ecld();
  fb = new G4Fb();
  fiss = new G4Fiss();
  opt = new G4Opt();
}

G4Abla::~G4Abla()
{
  delete pace;
  delete ald;
  delete ablamain;
  delete emdpar;
  delete eenuc;
  delete ec2sub;
  delete ecld;
  delete fb;
  delete fiss;
  delete opt;
}

// Main interface to the evaporation

// Possible problem with generic Geant4 interface: ABLA evaporation
// needs angular momentum information (calculated by INCL) to
// work. Maybe there is a way to obtain this information from
// G4Fragment?

void G4Abla::breakItUp(G4double nucleusA, G4double nucleusZ, G4double nucleusMass, G4double excitationEnergy,
		       G4double angularMomentum, G4double recoilEnergy, G4double momX, G4double momY, G4double momZ,
		       G4int eventnumber)
{
  const G4double uma = 931.4942;
  const G4double melec = 0.511;
  const G4double fmp = 938.27231;
  const G4double fmn = 939.56563;

  G4double alrem = 0.0, berem = 0.0, garem = 0.0;
  G4double R[4][4]; // Rotation matrix
  G4double csdir1[4];
  G4double csdir2[4];
  G4double csrem[4];
  G4double pfis_rem[4];
  G4double pf1_rem[4];
  for(G4int init_i = 0; init_i < 4; init_i++) {
    csdir1[init_i] = 0.0;
    csdir2[init_i] = 0.0;
    csrem[init_i] = 0.0;
    pfis_rem[init_i] = 0.0;
    pf1_rem[init_i] = 0.0;
    for(G4int init_j = 0; init_j < 4; init_j++) {
      R[init_i][init_j] = 0.0;
    }
  }

  G4double plab1 = 0.0, gam1 = 0.0, eta1 = 0.0;
  G4double plab2 = 0.0, gam2 = 0.0, eta2 = 0.0;

  G4double sitet = 0.0;
  G4double stet1 = 0.0;
  G4double stet2 = 0.0;
  
  G4int nbpevap = 0;
  G4int mempaw = 0, memiv = 0;

  G4double e_evapo = 0.0;
  G4double el = 0.0;
  G4double fmcv = 0.0;

  G4double aff1 = 0.0;
  G4double zff1 = 0.0;
  G4double eff1 = 0.0;
  G4double aff2 = 0.0;
  G4double zff2 = 0.0;
  G4double eff2 = 0.0;

  G4double v1 = 0.0, v2 = 0.0;

  G4double t2 = 0.0;
  G4double ctet1 = 0.0;
  G4double ctet2 = 0.0;
  G4double phi1 = 0.0;
  G4double phi2 = 0.0;
  G4double p2 = 0.0;
  G4double epf2_out = 0.0 ;
  G4int lma_pf1 = 0, lmi_pf1 = 0;
  G4int lma_pf2 = 0, lmi_pf2 = 0;
  G4int nopart = 0;

  G4double cst = 0.0, sst = 0.0, csf = 0.0, ssf = 0.0;
  
  G4double zf = 0.0, af = 0.0, mtota = 0.0, pleva = 0.0, pxeva = 0.0, pyeva = 0.0;
  G4int ff = 0;
  G4int inum = eventnumber;
  G4int inttype = 0;
  G4double esrem = excitationEnergy;
  
  G4double aprf = nucleusA;
  G4double zprf = nucleusZ;
  G4double mcorem = nucleusMass;
  G4double ee = excitationEnergy;
  G4double jprf = angularMomentum; // actually root-mean-squared

  G4double erecrem = recoilEnergy;
  G4double trem = 0.0;
  G4double pxrem = momX;
  G4double pyrem = momY;
  G4double pzrem = momZ;

  G4double remmass = 0.0;
  
  varntp->ntrack = 0;
  //  volant->iv = 0;
  volant->iv = 1;
  
  G4double pcorem = std::sqrt(erecrem*(erecrem +2.*938.2796*nucleusA));
    // G4double pcorem = std::sqrt(std::pow(momX,2) + std::pow(momY,2) + std::pow(momZ,2));
    // assert(isnan(pcorem) == false);
  if(esrem >= 1.0e-3) {
    evapora(zprf,aprf,ee,jprf, &zf, &af, &mtota, &pleva, &pxeva, &pyeva, &ff, &inttype, &inum);
  }
  else {
    ff = 0; 
    zf = zprf;
    af = aprf;
    pxeva = pxrem;
    pyeva = pyrem;
    pleva = pzrem;
  }
  // assert(isnan(zf) == false);
  // assert(isnan(af) == false);
  // assert(isnan(ee) == false);

  if (ff == 1) {
    // Fission:
    // variable ee: Energy of fissioning nucleus above the fission barrier.          
    // Calcul des impulsions des particules evaporees (avant fission) 
    // dans le systeme labo.

    trem = double(erecrem);
    remmass = pace2(aprf,zprf) + aprf*uma - zprf*melec; // canonic
//     remmass = mcorem  + double(esrem);			// ok
//     remmass = mcorem;					//cugnon
    varntp->kfis = 1;
    G4double gamrem = (remmass + trem)/remmass;
    G4double etrem = std::sqrt(trem*(trem + 2.0*remmass))/remmass;
    // assert(isnan(etrem) == false);

    //  This is not treated as accurately as for the non fission case for which
    //  the remnant mass is computed to satisfy the energy conservation 
    //  of evaporated particles. But it is not bad and more canonical!      
    remmass = pace2(aprf,zprf) + aprf*uma - zprf*melec+double(esrem); // !canonic
    //  Essais avec la masse de KHS (9/2002):
    el = 0.0;
    mglms(aprf,zprf,0,&el);
    remmass = zprf*fmp + (aprf-zprf)*fmn + el + double(esrem);

    gamrem = std::sqrt(std::pow(pcorem,2) + std::pow(remmass,2))/remmass;
    // assert(isnan(gamrem) == false);
    etrem = pcorem/remmass;
    // assert(isnan(etrem) == false);
    
    alrem = pxrem/pcorem;
    // assert(isnan(alrem) == false);
    berem = pyrem/pcorem;
    // assert(isnan(berem) == false);
    garem = pzrem/pcorem;
    // assert(isnan(garem) == false);

    csrem[0] = 0.0; // Should not be used.
    csrem[1] = alrem;
    csrem[2] = berem;
    csrem[3] = garem;
      
    // C Pour Vérif Remnant = evapo(Pre fission) + Noyau_fissionant (système  Remnant)
    G4double bil_e = 0.0;
    G4double bil_px = 0.0;
    G4double bil_py = 0.0;
    G4double bil_pz = 0.0;
    G4double masse = 0.0;
    
    for(G4int iloc = 1; iloc <= volant->iv; iloc++) { //DO iloc=1,iv
      // assert(isnan(volant->zpcv[iloc]) == false);
      //      assert(volant->acv[iloc] != 0);
      //      assert(volant->zpcv[iloc] != 0);
      mglms(double(volant->acv[iloc]),double(volant->zpcv[iloc]),0,&el);
      // assert(isnan(el) == false);
      masse = volant->zpcv[iloc]*fmp + (volant->acv[iloc] - volant->zpcv[iloc])*fmn + el;
      // assert(isnan(masse) == false);
      bil_e = bil_e + std::sqrt(std::pow(volant->pcv[iloc],2) + std::pow(masse,2));
      // assert(isnan(bil_e) == false); 
      bil_px = bil_px + volant->pcv[iloc]*(volant->xcv[iloc]);
      bil_py = bil_py + volant->pcv[iloc]*(volant->ycv[iloc]);
      bil_pz = bil_pz + volant->pcv[iloc]*(volant->zcv[iloc]);
    } //  enddo
    // C Ce bilan (impulsion nulle) est parfait. (Bil_Px=Bil_Px+PXEVA....)

    G4int ndec = 1;
    
    if(volant->iv != 0) { //then
      if(verboseLevel > 2) {
	G4cout <<"varntp->ntrack = " << varntp->ntrack << G4endl;
	G4cout <<"1st Translab: Adding indices from " << ndec << " to " << volant->iv << G4endl;
      }
      nopart = varntp->ntrack - 1;
      translab(gamrem,etrem,csrem,nopart,ndec);
      if(verboseLevel > 2) {
	G4cout <<"Translab complete!" << G4endl;
	G4cout <<"varntp->ntrack = " << varntp->ntrack << G4endl;
      }
    }
    nbpevap = volant->iv;	// nombre de particules d'evaporation traitees

    // C                                                                       
    // C Now calculation of the fission fragment distribution including                  
    // C evaporation from the fragments.                                           
    // C                                   

    // C Distribution of the fission fragments:
                                                                       
    //   void fissionDistri(G4double a,G4double z,G4double e,
    // 			   G4double &a1,G4double &z1,G4double &e1,G4double &v1,
    // 			     G4double &a2,G4double &z2,G4double &e2,G4double &v2);

    fissionDistri(af,zf,ee,aff1,zff1,eff1,v1,aff2,zff2,eff2,v2);

    // C verif des A et Z decimaux:
    G4int na_f = int(std::floor(af + 0.5));
    G4int nz_f = int(std::floor(zf + 0.5));
    varntp->izfis = nz_f;   // copie dans le ntuple
    varntp->iafis = na_f;
    G4int na_pf1 = int(std::floor(aff1 + 0.5));
    G4int nz_pf1 = int(std::floor(zff1 + 0.5));  
    G4int na_pf2 = int(std::floor(aff2 + 0.5));
    G4int nz_pf2 = int(std::floor(zff2 + 0.5));

    if((na_f != (na_pf1+na_pf2)) || (nz_f != (nz_pf1+nz_pf2))) {
      if(verboseLevel > 2) {
	G4cout <<"problemes arrondis dans la fission " << G4endl;
	G4cout << "af,zf,aff1,zff1,aff2,zff2" << G4endl;
	G4cout << af <<" , " << zf <<" , " << aff1 <<" , " << zff1 <<" , " << aff2 <<" , " << zff2 << G4endl;
	G4cout << "a,z,a1,z1,a2,z2 integer" << G4endl;
	G4cout << na_f <<" , " << nz_f <<" , " << na_pf1 <<" , " << nz_pf1 <<" , " << na_pf2 <<" , " << nz_pf2 << G4endl; 
      }
    }
    
    //  Calcul de l'impulsion des PF dans le syteme noyau de fission:
    G4int kboud = idnint(zf);                                                  
    G4int jboud = idnint(af-zf);                                             
    //G4double ef = fb->efa[kboud][jboud]; // barriere de fission
    G4double ef = fb->efa[jboud][kboud]; // barriere de fission
    // assert(isnan(ef) == false);
    varntp->estfis = ee + ef;   	// copie dans le ntuple   
     
    // C           MASSEF = pace2(AF,ZF)
    // C      	   MASSEF = MASSEF + AF*UMA - ZF*MELEC + EE + EF
    // C           MASSE1 = pace2(DBLE(AFF1),DBLE(ZFF1))
    // C      	   MASSE1 = MASSE1 + AFF1*UMA - ZFF1*MELEC + EFF1
    // C           MASSE2 = pace2(DBLE(AFF2),DBLE(ZFF2))
    // C      	   MASSE2 = MASSE2 + AFF2*UMA - ZFF2*MELEC + EFF2
    // C        WRITE(6,*) 'MASSEF,MASSE1,MASSE2',MASSEF,MASSE1,MASSE2
    // C MGLMS est la fonction de masse cohérente avec KHS evapo-fis.
    // C   Attention aux parametres, ici 0=OPTSHP, NO microscopic correct. 
    mglms(af,zf,0,&el);
    // assert(isnan(el) == false);
    G4double massef = zf*fmp + (af - zf)*fmn + el + ee + ef;
    // assert(isnan(massef) == false);
    mglms(double(aff1),double(zff1),0,&el);
    // assert(isnan(el) == false);
    G4double masse1 = zff1*fmp + (aff1-zff1)*fmn + el + eff1;
    // assert(isnan(masse1) == false);
    mglms(aff2,zff2,0,&el);
    // assert(isnan(el) == false);
    G4double masse2 = zff2*fmp + (aff2 - zff2)*fmn + el + eff2;
    // assert(isnan(masse2) == false);
    // C        WRITE(6,*) 'MASSEF,MASSE1,MASSE2',MASSEF,MASSE1,MASSE2	   
    G4double b = massef - masse1 - masse2;
    if(b < 0.0) { //then
      b=0.0;
      if(verboseLevel > 2) {
	G4cout <<"anomalie dans la fission: " << G4endl; 
	G4cout << inum<< " , " << af<< " , " <<zf<< " , " <<massef<< " , " <<aff1<< " , " <<zff1<< " , " <<masse1<< " , " <<aff2<< " , " <<zff2<< " , " << masse2 << G4endl;
      }
    } //endif
    G4double t1 = b*(b + 2.0*masse2)/(2.0*massef);
    // assert(isnan(t1) == false);
    G4double p1 = std::sqrt(t1*(t1 + 2.0*masse1));
    // assert(isnan(p1) == false);
    
    G4double rndm;
    standardRandom(&rndm, &(hazard->igraine[13]));
    ctet1 = 2.0*rndm - 1.0;
    standardRandom(&rndm,&(hazard->igraine[9]));
    phi1 = rndm*2.0*3.141592654;
           
    // C ----Coefs de la transformation de Lorentz (noyau de fission -> Remnant) 
    G4double peva = std::pow(pxeva,2) + std::pow(pyeva,2) + std::pow(pleva,2);
    G4double gamfis = std::sqrt(std::pow(massef,2) + peva)/massef;
    // assert(isnan(gamfis) == false);
    peva = std::sqrt(peva);
    // assert(isnan(peva) == false);
    G4double etfis = peva/massef;
      
    G4double epf1_in = 0.0;
    G4double epf1_out = 0.0;

    // C ----Matrice de rotation (noyau de fission -> Remnant)
    if(peva >= 1.0e-4) {
      sitet = std::sqrt(std::pow(pxeva,2)+std::pow(pyeva,2))/peva;
      // assert(isnan(sitet) == false);
    }
    if(sitet > 1.0e-5) { //then
      G4double cstet = pleva/peva;
      G4double siphi = pyeva/(sitet*peva);
      G4double csphi = pxeva/(sitet*peva);
	
      R[1][1] = cstet*csphi;
      R[1][2] = -siphi;
      R[1][3] = sitet*csphi;
      R[2][1] = cstet*siphi;
      R[2][2] = csphi;
      R[2][3] = sitet*siphi;
      R[3][1] = -sitet;
      R[3][2] = 0.0;
      R[3][3] = cstet;
    }
    else {
      R[1][1] = 1.0;
      R[1][2] = 0.0;
      R[1][3] = 0.0;
      R[2][1] = 0.0;
      R[2][2] = 1.0;
      R[2][3] = 0.0;
      R[3][1] = 0.0;
      R[3][2] = 0.0;
      R[3][3] = 1.0;
    } //  endif
    // c test de verif:                                      

    if((zff1 <= 0.0) || (aff1 <= 0.0) || (aff1 < zff1)) { //then   
      if(verboseLevel > 2) {
	G4cout <<"zf = " <<  zf <<" af = " << af <<"ee = " << ee <<"zff1 = " << zff1 <<"aff1 = " << aff1 << G4endl;
      }
    }
    else {
      // C ---------------------- PF1 will evaporate 
      epf1_in = double(eff1);
      epf1_out = epf1_in;
      //   void evapora(G4double zprf, G4double aprf, G4double ee, G4double jprf, 
      // 	       G4double *zf_par, G4double *af_par, G4double *mtota_par,
      // 	       G4double *pleva_par, G4double *pxeva_par, G4double *pyeva_par,
      // 	       G4double *ff_par, G4int *inttype_par, G4int *inum_par);
      G4double zf1, af1, malpha1, ffpleva1, ffpxeva1, ffpyeva1;
      G4int ff1, ftype1;
      evapora(zff1, aff1, epf1_out, 0.0, &zf1, &af1, &malpha1, &ffpleva1,
	      &ffpxeva1, &ffpyeva1, &ff1, &ftype1, &inum);
      // C On ajoute le fragment:
      // assert(af1 > 0);
      volant->iv = volant->iv + 1;
      // assert(af1 != 0);
      // assert(zf1 != 0);
      volant->acv[volant->iv] = af1;
      volant->zpcv[volant->iv] = zf1;
      if(verboseLevel > 2) {
	G4cout <<"Added fission fragment: a = " << volant->acv[volant->iv] << " z = " << volant->zpcv[volant->iv] << G4endl;
      }
      peva = std::sqrt(std::pow(ffpxeva1,2) + std::pow(ffpyeva1,2) + std::pow(ffpleva1,2));
      // assert(isnan(peva) == false);
      volant->pcv[volant->iv] = peva;
      if(peva > 0.001) { // then
	volant->xcv[volant->iv] = ffpxeva1/peva;
	volant->ycv[volant->iv] = ffpyeva1/peva;
	volant->zcv[volant->iv] = ffpleva1/peva;
      }
      else {
	volant->xcv[volant->iv] = 1.0;
	volant->ycv[volant->iv] = 0.0;
	volant->zcv[volant->iv] = 0.0;
      } // end if
	        
      // C Pour Vérif evapo de PF1 dans le systeme du Noyau Fissionant
      G4double bil1_e = 0.0;
      G4double bil1_px = 0.0;
      G4double bil1_py=0.0;
      G4double bil1_pz=0.0;
      for(G4int iloc = nbpevap + 1; iloc <= volant->iv; iloc++) { //do iloc=nbpevap+1,iv
	//      for(G4int iloc = nbpevap + 1; iloc <= volant->iv + 1; iloc++) { //do iloc=nbpevap+1,iv
	mglms(volant->acv[iloc], volant->zpcv[iloc],0,&el);
       	masse = volant->zpcv[iloc]*fmp + (volant->acv[iloc] - volant->zpcv[iloc])*fmn + el; 
	// assert(isnan(masse) == false);
 	bil1_e = bil1_e + std::sqrt(std::pow(volant->pcv[iloc],2) + std::pow(masse,2));
	// assert(isnan(bil1_e) == false);
 	bil1_px = bil1_px + volant->pcv[iloc]*(volant->xcv[iloc]);
 	bil1_py = bil1_py + volant->pcv[iloc]*(volant->ycv[iloc]);
 	bil1_pz = bil1_pz + volant->pcv[iloc]*(volant->zcv[iloc]);
      } // enddo

      // Calcul des cosinus directeurs de PF1 dans le Remnant et calcul
      // des coefs pour la transformation de Lorentz Systeme PF --> Systeme Remnant
      translabpf(masse1,t1,p1,ctet1,phi1,gamfis,etfis,R,&plab1,&gam1,&eta1,csdir1);

      // calcul des impulsions des particules evaporees dans le systeme Remnant:
      if(verboseLevel > 2) {
	G4cout <<"2nd Translab (pf1 evap): Adding indices from " << nbpevap+1 << " to " << volant->iv << G4endl;
      }
      nopart = varntp->ntrack - 1;
      translab(gam1,eta1,csdir1,nopart,nbpevap+1);
      if(verboseLevel > 2) {
	G4cout <<"After translab call... varntp->ntrack = " << varntp->ntrack << G4endl;
      }
      memiv = nbpevap + 1;	  // memoires pour la future transformation
      mempaw = nopart;	  // remnant->labo pour pf1 et pf2.
      lmi_pf1 = nopart + nbpevap + 1;  // indices min et max dans /var_ntp/
      lma_pf1 = nopart + volant->iv;   // des particules issues de pf1
      nbpevap = volant->iv;	// nombre de particules d'evaporation traitees
    } // end if
    // C --------------------- End of PF1 calculation
  
    // c test de verif:                                                                                                         
    if((zff2 <= 0.0) || (aff2 <= 0.0) || (aff2 <= zff2)) { //then   
      if(verboseLevel > 2) {
	G4cout << zf << " " << af << " " << ee  << " " << zff2 << " " << aff2 << G4endl;                                
      }
    }
    else {                                                          
      // C ---------------------- PF2 will evaporate 
      G4double epf2_in = double(eff2);
      G4double epf2_out = epf2_in;
      //   void evapora(G4double zprf, G4double aprf, G4double ee, G4double jprf, 
      // 	       G4double *zf_par, G4double *af_par, G4double *mtota_par,
      // 	       G4double *pleva_par, G4double *pxeva_par, G4double *pyeva_par,
      // 	       G4double *ff_par, G4int *inttype_par, G4int *inum_par);
      G4double zf2, af2, malpha2, ffpleva2, ffpxeva2, ffpyeva2;
      G4int ff2, ftype2;
      evapora(zff2,aff2,epf2_out,0.0,&zf2,&af2,&malpha2,&ffpleva2,
	      &ffpxeva2,&ffpyeva2,&ff2,&ftype2,&inum);
      // C On ajoute le fragment:
      volant->iv = volant->iv + 1;
      volant->acv[volant->iv] = af2;
      volant->zpcv[volant->iv] = zf2; 
      if(verboseLevel > 2) {
	G4cout <<"Added fission fragment: a = " << volant->acv[volant->iv] << " z = " << volant->zpcv[volant->iv] << G4endl;
      }
      peva = std::sqrt(std::pow(ffpxeva2,2) + std::pow(ffpyeva2,2) + std::pow(ffpleva2,2));
      // assert(isnan(peva) == false);
      volant->pcv[volant->iv] = peva;
      //      exit(0);
      if(peva > 0.001) { //then
	volant->xcv[volant->iv] = ffpxeva2/peva;
	volant->ycv[volant->iv] = ffpyeva2/peva;
	volant->zcv[volant->iv] = ffpleva2/peva;
      }
      else {
	volant->xcv[volant->iv] = 1.0;
	volant->ycv[volant->iv] = 0.0;
	volant->zcv[volant->iv] = 0.0;
      } //end if        
      // C Pour Vérif evapo de PF1 dans le systeme du Noyau Fissionant
      G4double bil2_e = 0.0;
      G4double bil2_px = 0.0;
      G4double bil2_py = 0.0;
      G4double bil2_pz = 0.0;
      //      for(G4int iloc = nbpevap + 1; iloc <= volant->iv; iloc++) { //do iloc=nbpevap+1,iv
      for(G4int iloc = nbpevap + 1; iloc <= volant->iv; iloc++) { //do iloc=nbpevap+1,iv
	mglms(volant->acv[iloc],volant->zpcv[iloc],0,&el);
      	masse = volant->zpcv[iloc]*fmp + (volant->acv[iloc] - volant->zpcv[iloc])*fmn + el; 
	bil2_e = bil2_e + std::sqrt(std::pow(volant->pcv[iloc],2) + std::pow(masse,2));
	// assert(isnan(bil2_e) == false);
	bil2_px = bil2_px + volant->pcv[iloc]*(volant->xcv[iloc]);
	bil2_py = bil2_py + volant->pcv[iloc]*(volant->ycv[iloc]);
	bil2_pz = bil2_pz + volant->pcv[iloc]*(volant->zcv[iloc]);
      } //enddo

      // C ----Calcul des cosinus directeurs de PF2 dans le Remnant et calcul
      // c des coefs pour la transformation de Lorentz Systeme PF --> Systeme Remnant
      G4double t2 = b - t1;
      //      G4double ctet2 = -ctet1;
      ctet2 = -1.0*ctet1;
      assert(std::fabs(ctet2) <= 1.0);
      // assert(isnan(ctet2) == false);
      phi2 = dmod(phi1+3.141592654,6.283185308);
      // assert(isnan(phi2) == false);
      G4double p2 = std::sqrt(t2*(t2+2.0*masse2));
      // assert(isnan(p2) == false);
      
      //   void translabpf(G4double masse1, G4double t1, G4double p1, G4double ctet1,
      // 		  G4double phi1, G4double gamrem, G4double etrem, G4double R[][4],
      // 		  G4double *plab1, G4double *gam1, G4double *eta1, G4double csdir[]);
      translabpf(masse2,t2,p2,ctet2,phi2,gamfis,etfis,R,&plab2,&gam2,&eta2,csdir2);
      // C
      // C  calcul des impulsions des particules evaporees dans le systeme Remnant:
      // c
      if(verboseLevel > 2) {
	G4cout <<"3rd Translab (pf2 evap): Adding indices from " << nbpevap+1 << " to " << volant->iv << G4endl;
      }
      nopart = varntp->ntrack - 1;
      translab(gam2,eta2,csdir2,nopart,nbpevap+1);
      lmi_pf2 = nopart + nbpevap + 1;	// indices min et max dans /var_ntp/
      lma_pf2 = nopart + volant->iv;		// des particules issues de pf2
    } // end if
    // C --------------------- End of PF2 calculation

    // C Pour vérifications: calculs du noyau fissionant et des PF dans 
    // C    le systeme du remnant.
    for(G4int iloc = 1; iloc <= 3; iloc++) { // do iloc=1,3
      pfis_rem[iloc] = 0.0;
    } // enddo 
    G4double efis_rem, pfis_trav[4];
    lorab(gamfis,etfis,massef,pfis_rem,&efis_rem,pfis_trav);
    rotab(R,pfis_trav,pfis_rem);
      
    stet1 = std::sqrt(1.0 - std::pow(ctet1,2));
    // assert(isnan(stet1) == false);
    pf1_rem[1] = p1*stet1*std::cos(phi1);
    pf1_rem[2] = p1*stet1*std::sin(phi1);
    pf1_rem[3] = p1*ctet1;
    G4double e1_rem;
    lorab(gamfis,etfis,masse1+t1,pf1_rem,&e1_rem,pfis_trav);
    rotab(R,pfis_trav,pf1_rem);

    stet2 = std::sqrt(1.0 - std::pow(ctet2,2));
    assert(std::pow(ctet2,2) >= 0.0);
    assert(std::pow(ctet2,2) <= 1.0);
    // assert(isnan(stet2) == false);
    
    G4double pf2_rem[4];
    G4double e2_rem;
    pf2_rem[1] = p2*stet2*std::cos(phi2);
    pf2_rem[2] = p2*stet2*std::sin(phi2);
    pf2_rem[3] = p2*ctet2;
    lorab(gamfis,etfis,masse2+t2,pf2_rem,&e2_rem,pfis_trav);
    rotab(R,pfis_trav,pf2_rem);
    // C Verif 0: Remnant = evapo_pre_fission + Noyau Fissionant
    bil_e = remmass - efis_rem - bil_e;
    bil_px = bil_px + pfis_rem[1];
    bil_py = bil_py + pfis_rem[2];  
    bil_pz = bil_pz + pfis_rem[3];  
    // C Verif 1: noyau fissionant = PF1 + PF2 dans le systeme remnant
    //    G4double bilan_e = efis_rem - e1_rem - e2_rem;
    //    G4double bilan_px = pfis_rem[1] - pf1_rem[1] - pf2_rem[1];
    //    G4double bilan_py = pfis_rem[2] - pf1_rem[2] - pf2_rem[2];
    //    G4double bilan_pz = pfis_rem[3] - pf1_rem[3] - pf2_rem[3];
    // C Verif 2: PF1 et PF2 egaux a toutes leurs particules evaporees
    // C   (Systeme remnant)
    if((lma_pf1-lmi_pf1) != 0) { //then
      G4double bil_e_pf1 = e1_rem - epf1_out;
      G4double bil_px_pf1 = pf1_rem[1];
      G4double bil_py_pf1 = pf1_rem[2]; 
      G4double bil_pz_pf1 = pf1_rem[3];
      for(G4int ipf1 = lmi_pf1; ipf1 <= lma_pf1; ipf1++) { //do ipf1=lmi_pf1,lma_pf1
	bil_e_pf1 = bil_e_pf1 - (std::pow(varntp->plab[ipf1],2) + std::pow(varntp->enerj[ipf1],2))/(2.0*(varntp->enerj[ipf1]));
	cst = std::cos(varntp->tetlab[ipf1]/57.2957795);
	sst = std::sin(varntp->tetlab[ipf1]/57.2957795);
	csf = std::cos(varntp->philab[ipf1]/57.2957795);
	ssf = std::sin(varntp->philab[ipf1]/57.2957795);
	bil_px_pf1 = bil_px_pf1 - varntp->plab[ipf1]*sst*csf;
	bil_py_pf1 = bil_py_pf1 - varntp->plab[ipf1]*sst*ssf;
	bil_pz_pf1 = bil_pz_pf1 - varntp->plab[ipf1]*cst;		 
      } // enddo
    } //endif
	 
    if((lma_pf2-lmi_pf2) != 0) { //then
      G4double bil_e_pf2 =  e2_rem - epf2_out;
      G4double bil_px_pf2 = pf2_rem[1];
      G4double bil_py_pf2 = pf2_rem[2]; 
      G4double bil_pz_pf2 = pf2_rem[3];
      for(G4int ipf2 = lmi_pf2; ipf2 <= lma_pf2; ipf2++) { //do ipf2=lmi_pf2,lma_pf2
	bil_e_pf2 = bil_e_pf2 - (std::pow(varntp->plab[ipf2],2) + std::pow(varntp->enerj[ipf2],2))/(2.0*(varntp->enerj[ipf2]));
	G4double cst = std::cos(varntp->tetlab[ipf2]/57.2957795);
	G4double sst = std::sin(varntp->tetlab[ipf2]/57.2957795);
	G4double csf = std::cos(varntp->philab[ipf2]/57.2957795);
	G4double ssf = std::sin(varntp->philab[ipf2]/57.2957795);
	bil_px_pf2 = bil_px_pf2 - varntp->plab[ipf2]*sst*csf;
	bil_py_pf2 = bil_py_pf2 - varntp->plab[ipf2]*sst*ssf;
	bil_pz_pf2 = bil_pz_pf2 - varntp->plab[ipf2]*cst;		 
      } // enddo
    } //endif 
    // C
    // C ---- Transformation systeme Remnant -> systeme labo. (evapo de PF1 ET PF2)
    // C
    //    G4double mempaw, memiv;
    if(verboseLevel > 2) {
      G4cout <<"4th Translab: Adding indices from " << memiv << " to " << volant->iv << G4endl;
    }
    translab(gamrem,etrem,csrem,mempaw,memiv);
    // C *******************  END of fission calculations ************************
  }
  else {
    // C ************************ Evapo sans fission *****************************
    // C Here, FF=0, --> Evapo sans fission, on ajoute le fragment:
    // C *************************************************************************
    varntp->kfis = 0;
    if(verboseLevel > 2) {
      G4cout <<"Evaporation without fission" << G4endl;
    }
    volant->iv = volant->iv + 1;
    volant->acv[volant->iv] = af;
    volant->zpcv[volant->iv] = zf;
    G4double peva = std::sqrt(std::pow(pxeva,2)+std::pow(pyeva,2)+std::pow(pleva,2));
    // assert(isnan(peva) == false);
    volant->pcv[volant->iv] = peva;
    if(peva > 0.001) { //then
      volant->xcv[volant->iv] = pxeva/peva;
      volant->ycv[volant->iv] = pyeva/peva;
      volant->zcv[volant->iv] = pleva/peva;        
    }
    else {
      volant->xcv[volant->iv] = 1.0;
      volant->ycv[volant->iv] = 0.0;
      volant->zcv[volant->iv] = 0.0;
    } // end if        
	
    // C
    // C  calcul des impulsions des particules evaporees dans le systeme labo:
    // c
    trem = double(erecrem);
    // C      REMMASS = pace2(APRF,ZPRF) + APRF*UMA - ZPRF*MELEC	!Canonic
    // C      REMMASS = MCOREM  + DBLE(ESREM)				!OK
    remmass = mcorem;						//Cugnon
    // C      GAMREM = (REMMASS + TREM)/REMMASS			!OK
    // C      ETREM = DSQRT(TREM*(TREM + 2.*REMMASS))/REMMASS		!OK
    csrem[0] = 0.0; // Should not be used.
    csrem[1] = alrem;
    csrem[2] = berem;
    csrem[3] = garem;

    //    for(G4int j = 1; j <= volant->iv; j++) { //do j=1,iv
    for(G4int j = 1; j <= volant->iv; j++) { //do j=1,iv
      if(volant->acv[j] == 0) {
	if(verboseLevel > 2) {
	  G4cout <<"volant->acv[" << j << "] = 0" << G4endl;
	  G4cout <<"volant->iv = " << volant->iv << G4endl;
	}
      }
      if(volant->acv[j] > 0) {
	assert(volant->acv[j] != 0);
	//      assert(volant->zpcv[j] != 0);
	mglms(volant->acv[j],volant->zpcv[j],0,&el);
	fmcv = volant->zpcv[j]*fmp + (volant->acv[j] - volant->zpcv[j])*fmn + el;
	e_evapo = e_evapo + std::sqrt(std::pow(volant->pcv[j],2) + std::pow(fmcv,2));
	// assert(isnan(e_evapo) == false);
      }
    } // enddo

    // C Redefinition pour conservation d'impulsion!!!
    // C   this mass obtained by energy balance is very close to the
    // C   mass of the remnant computed by pace2 + excitation energy (EE). (OK)      
    remmass = e_evapo;
      
    G4double gamrem = std::sqrt(std::pow(pcorem,2)+std::pow(remmass,2))/remmass;
    // assert(isnan(gamrem) == false);
    G4double etrem = pcorem/remmass;

    if(verboseLevel > 2) {
      G4cout <<"5th Translab (no fission): Adding indices from " << 1 << " to " << volant->iv << G4endl;
    }
    nopart = varntp->ntrack - 1;
    translab(gamrem,etrem,csrem,nopart,1);
                  
    // C End of the (FISSION - NO FISSION) condition (FF=1 or 0)                                          
  } //end if 
  // C *********************** FIN de l'EVAPO KHS ******************** 
}

// Evaporation code
void G4Abla::initEvapora()
{
  //     37	C     PROJECTILE AND TARGET PARAMETERS + CROSS SECTIONS                 
  //     38	C     COMMON /ABLAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,       
  //     39	C                       R_0,R_P,R_T, IMAX,IRNDM,PI,                     
  //     40	C                       BFPRO,SNPRO,SPPRO,SHELL                         
  //     41	C                                                                       
  //     42	C     AP,ZP,AT,ZT   - PROJECTILE AND TARGET MASSES                      
  //     43	C     EAP,BETA      - BEAM ENERGY PER NUCLEON, V/C                      
  //     44	C     BMAXNUC       - MAX. IMPACT PARAMETER FOR NUCL. REAC.             
  //     45	C     CRTOT,CRNUC   - TOTAL AND NUCLEAR REACTION CROSS SECTION          
  //     46	C     R_0,R_P,R_T,  - RADIUS PARAMETER, PROJECTILE+ TARGET RADII        
  //     47	C     IMAX,IRNDM,PI - MAXIMUM NUMBER OF EVENTS, DUMMY, 3.141...         
  //     48	C     BFPRO         - FISSION BARRIER OF THE PROJECTILE                 
  //     49	C     SNPRO         - NEUTRON SEPARATION ENERGY OF THE PROJECTILE       
  //     50	C     SPPRO         - PROTON    "           "   "    "   "              
  //     51	C     SHELL         - GROUND STATE SHELL CORRECTION                     
  //     52	C---------------------------------------------------------------------  
  //     53	C                                                                       
  //     54	C     ENERGIES WIDTHS AND CROSS SECTIONS FOR EM EXCITATION              
  //     55	C     COMMON /EMDPAR/ EGDR,EGQR,FWHMGDR,FWHMGQR,CREMDE1,CREMDE2,        
  //     56	C                     AE1,BE1,CE1,AE2,BE2,CE2,SR1,SR2,XR                
  //     57	C                                                                       
  //     58	C     EGDR,EGQR       - MEAN ENERGY OF GDR AND GQR                      
  //     59	C     FWHMGDR,FWHMGQR - FWHM OF GDR, GQR                                
  //     60	C     CREMDE1,CREMDE2 - EM CROSS SECTION FOR E1 AND E2                  
  //     61	C     AE1,BE1,CE1     - ARRAYS TO CALCULATE                             
  //     62	C     AE2,BE2,CE2     - THE EXCITATION ENERGY AFTER E.M. EXC.           
  //     63	C     SR1,SR2,XR      - WITH MONTE CARLO                                
  //     64	C---------------------------------------------------------------------  
  //     65	C                                                                       
  //     66	C     DEFORMATIONS AND G.S. SHELL EFFECTS                               
  //     67	C     COMMON /ECLD/   ECGNZ,ECFNZ,VGSLD,ALPHA                           
  //     68	C                                                                       
  //     69	C     ECGNZ - GROUND STATE SHELL CORR. FRLDM FOR A SPHERICAL G.S.       
  //     70	C     ECFNZ - SHELL CORRECTION FOR THE SADDLE POINT (NOW: == 0)         
  //     71	C     VGSLD - DIFFERENCE BETWEEN DEFORMED G.S. AND LDM VALUE            
  //     72	C     ALPHA - ALPHA GROUND STATE DEFORMATION (THIS IS NOT BETA2!)       
  //     73	C             BETA2 = SQRT(5/(4PI)) * ALPHA                             
  //     74	C---------------------------------------------------------------------  
  //     75	C                                                                       
  //     76	C     ARRAYS FOR EXCITATION ENERGY BY STATISTICAL HOLE ENERY MODEL      
  //     77	C     COMMON /EENUC/  SHE, XHE                                          
  //     78	C                                                                       
  //     79	C     SHE, XHE - ARRAYS TO CALCULATE THE EXC. ENERGY AFTER              
  //     80	C                ABRASION BY THE STATISTICAL HOLE ENERGY MODEL          
  //     81	C---------------------------------------------------------------------  
  //     82	C                                                                       
  //     83	C     G.S. SHELL EFFECT                                                 
  //     84	C     COMMON /EC2SUB/ ECNZ                                              
  //     85	C                                                                       
  //     86	C     ECNZ G.S. SHELL EFFECT FOR THE MASSES (IDENTICAL TO ECGNZ)        
  //     87	C---------------------------------------------------------------------  
  //     88	C                                                                       
  //     89	C     OPTIONS AND PARAMETERS FOR FISSION CHANNEL                        
  //     90	C     COMMON /FISS/    AKAP,BET,HOMEGA,KOEFF,IFIS,                       
  //     91	C                            OPTSHP,OPTXFIS,OPTLES,OPTCOL               
  //     92	C                                                                       
  //     93	C     AKAP   - HBAR**2/(2* MN * R_0**2) = 10 MEV                        
  //     94	C     BET    - REDUCED NUCLEAR FRICTION COEFFICIENT IN (10**21 S**-1)   
  //     95	C     HOMEGA - CURVATURE OF THE FISSION BARRIER = 1 MEV                 
  //     96	C     KOEFF  - COEFFICIENT FOR THE LD FISSION BARRIER == 1.0            
  //     97	C     IFIS   - 0/1 FISSION CHANNEL OFF/ON                               
  //     98	C     OPTSHP - INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
  //     99	C            = 0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY        
  //    100	C            = 1 SHELL ,  NO PAIRING                                    
  //    101	C            = 2 PAIRING, NO SHELL                                      
  //    102	C            = 3 SHELL AND PAIRING                                      
  //    103	C     OPTCOL - 0/1 COLLECTIVE ENHANCEMENT SWITCHED ON/OFF               
  //    104	C     OPTXFIS- 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV         
  //    105	C              FISSILITY PARAMETER.                                     
  //    106	C     OPTLES - CONSTANT TEMPERATURE LEVEL DENSITY FOR A,Z > TH-224      
  //    107	C     OPTCOL - 0/1 COLLECTIVE ENHANCEMENT OFF/ON                        
  //    108	C---------------------------------------------------------------------  
  //    109	C                                                                       
  //    110	C     OPTIONS                                                           
  //    111	C     COMMON /OPT/    OPTEMD,OPTCHA,EEFAC                               
  //    112	C                                                                       
  //    113	C     OPTEMD - 0/1  NO EMD / INCL. EMD                                  
  //    114	C     OPTCHA - 0/1  0 GDR / 1 HYPERGEOMETRICAL PREFRAGMENT-CHARGE-DIST. 
  //    115	C              ***  RECOMMENDED IS OPTCHA = 1 ***                       
  //    116	C     EEFAC  - EXCITATION ENERGY FACTOR, 2.0 RECOMMENDED                
  //    117	C---------------------------------------------------------------------  
  //    118	C                                                                       
  //    119	C     FISSION BARRIERS                                                  
  //    120	C     COMMON /FB/     EFA                                               
  //    121	C     EFA    - ARRAY OF FISSION BARRIERS                                
  //    122	C---------------------------------------------------------------------  
  //    123	C                                                                       
  //    124	C p    LEVEL DENSITY PARAMETERS                                          
  //    125	C     COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
  //    126	C     AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
  //    127	C                LEVEL DENSITY PARAMETER                                
  //    128	C     OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
  //    129	C               RECOMMENDED IS OPTAFAN = 0                              
  //    130	C---------------------------------------------------------------------  
  //    131	C   ____________________________________________________________________
  //    132	C  /                                                                    
  //    133	C  /  INITIALIZES PARAMETERS IN COMMON /ABRAMAIN/, /EMDPAR/, /ECLD/ ... 
  //    134	C  /  PROJECTILE PARAMETERS, EMD PARAMETERS, SHELL CORRECTION TABLES.   
  //    135	C  /  CALCULATES MAXIMUM IMPACT PARAMETER FOR NUCLEAR COLLISIONS AND    
  //    136	C  /  TOTAL GEOMETRICAL CROSS SECTION + EMD CROSS SECTIONS              
  //    137	C   ____________________________________________________________________
  //    138	C                                                                       
  //    139	C                                                                       
  //    201	C                                                                       
  //    202	C---------- SET INPUT VALUES                                            
  //    203	C                                                                       
  //    204	C *** INPUT FROM UNIT 10 IN THE FOLLOWING SEQUENCE !                    
  //    205	C     AP1 =    INTEGER  !                                               
  //    206	C     ZP1 =    INTEGER  !                                               
  //    207	C     AT1 =    INTEGER  !                                               
  //    208	C     ZT1 =    INTEGER  !                                               
  //    209	C     EAP =    REAL     !                                               
  //    210	C     IMAX =   INTEGER  !                                               
  //    211	C     IFIS =   INTEGER SWITCH FOR FISSION                               
  //    212	C     OPTSHP = INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
  //    213	C            =0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY         
  //    214	C            =1 SHELL , NO PAIRING CORRECTION                           
  //    215	C            =2 PAIRING, NO SHELL CORRECTION                            
  //    216	C            =3 SHELL AND PAIRING CORRECTION IN MASSES AND ENERGY       
  //    217	C     OPTEMD =0,1  0 NO EMD, 1 INCL. EMD                                
  //    218	C               ELECTROMAGNETIC DISSOZIATION IS CALCULATED AS WELL.     
  //    219	C     OPTCHA =0,1  0 GDR- , 1 HYPERGEOMETRICAL PREFRAGMENT-CHARGE-DIST. 
  //    220	C               RECOMMENDED IS OPTCHA=1                                 
  //    221	C     OPTCOL =0,1 COLLECTIVE ENHANCEMENT SWITCHED ON 1 OR OFF 0 IN DENSN
  //    222	C     OPTAFAN=0,1 SWITCH FOR AF/AN = 1 IN DENSNIV 0 AF/AN>1 1 AF/AN=1   
  //    223	C     AKAP =  REAL    ALWAYS EQUALS 10                                  
  //    224	C     BET  =  REAL    REDUCED FRICTION COEFFICIENT / 10**(+21) S**(-1)  
  //    225	C     HOMEGA = REAL   CURVATURE / MEV RECOMMENDED = 1. MEV              
  //    226	C     KOEFF  = REAL   COEFFICIENT FOR FISSION BARRIER                   
  //    227	C     OPTXFIS= INTEGER 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV 
  //    228	C              FISSILITY PARAMETER.                                     
  //    229	C     EEFAC  = REAL EMPIRICAL FACTOR FOR THE EXCITATION ENERGY          
  //    230	C                   RECOMMENDED 2.D0, STATISTICAL ABRASION MODELL 1.D0  
  //    231	C     AV     = REAL KOEFFICIENTS FOR CALCULATION OF A(TILDE)            
  //    232	C     AS     = REAL LEVEL DENSITY PARAMETER                             
  //    233	C     AK     = REAL                                                     
  //    234	C                                                                       
  //    235	C This following inputs will be initialized in the main through the 
  //    236	C         common /ABLAMAIN/  (A.B.)
  //    237	

  // switch-fission.1=on.0=off
  fiss->ifis = 1;

  // shell+pairing.0-1-2-3
  fiss->optshp = 0;

  // optemd =0,1  0 no emd, 1 incl. emd                                
  opt->optemd = 1;
  // read(10,*,iostat=io) dum(10),optcha                               
  opt->optcha = 1;

  // not.to.be.changed.(akap)
  fiss->akap = 10.0;

  // nuclear.viscosity.(beta)
  fiss->bet = 1.5;

  // potential-curvature
  fiss->homega = 1.0;

  // fission-barrier-coefficient
  fiss->koeff = 1.;

  //collective enhancement switched on 1 or off 0 in densn (qr=val or =1.)
  fiss->optcol = 0;

  // switch-for-low-energy-sys
  fiss->optles = 0;

  opt->eefac = 2.;

  ald->optafan = 0;

  ald->av = 0.073e0;
  ald->as = 0.095e0;
  ald->ak = 0.0e0;

  if(verboseLevel > 3) {
    G4cout <<"ifis " << fiss->ifis << G4endl;
    G4cout <<"optshp " << fiss->optshp << G4endl;
    G4cout <<"optemd " << opt->optemd << G4endl;
    G4cout <<"optcha " << opt->optcha << G4endl;
    G4cout <<"akap " << fiss->akap << G4endl;
    G4cout <<"bet " << fiss->bet << G4endl;
    G4cout <<"homega " << fiss->homega << G4endl;
    G4cout <<"koeff " << fiss->koeff << G4endl;
    G4cout <<"optcol " << fiss->optcol << G4endl;
    G4cout <<"optles " << fiss->optles << G4endl;
    G4cout <<"eefac " << opt->eefac << G4endl;
    G4cout <<"optafan " << ald->optafan << G4endl;
    G4cout <<"av " << ald->av << G4endl;
    G4cout <<"as " << ald->as << G4endl;
    G4cout <<"ak " << ald->ak << G4endl;
  }
  fiss->optxfis = 1;

  G4InclAblaDataFile *dataInterface = new G4InclAblaDataFile();
  if(dataInterface->readData() == true) {
    if(verboseLevel > 0) {
      G4cout <<"G4Abla: Datafiles read successfully." << G4endl;
    }
  }
  else {
    G4Exception("ERROR: Failed to read datafiles.");
  }
  
  for(int z = 0; z < 98; z++) { //do 30  z = 0,98,1                                                 
    for(int n = 0; n < 154; n++) { //do 31  n = 0,153,1                                              
      ecld->ecfnz[n][z] = 0.e0;
      ec2sub->ecnz[n][z] = dataInterface->getEcnz(n,z);
      ecld->ecgnz[n][z] = dataInterface->getEcnz(n,z);
      ecld->alpha[n][z] = dataInterface->getAlpha(n,z);
      ecld->vgsld[n][z] = dataInterface->getVgsld(n,z);
    }
  }

  for(int z = 0; z < 500; z++) {
    for(int a = 0; a < 500; a++) {
      pace->dm[z][a] = dataInterface->getPace2(z,a);
    }
  }

  delete dataInterface;
}

void G4Abla::qrot(G4double z, G4double a, G4double bet, G4double sig, G4double u, G4double *qr)
{
  G4double ucr = 10.0; // Critical energy for damping.
  G4double dcr = 40.0; // Width of damping.
  G4double ponq = 0.0, dn = 0.0, n = 0.0, dz = 0.0;

  if(((std::fabs(bet)-1.15) < 0) || ((std::fabs(bet)-1.15) == 0)) {
    goto qrot10;
  }

  if((std::fabs(bet)-1.15) > 0) {
    goto qrot11;
  }

 qrot10:
  n = a - z;
  dz = std::fabs(z - 82.0);
  if (n > 104) {
    dn = std::fabs(n-126.e0);
  }
  else {
    dn = std::fabs(n - 82.0);
  }

  bet = 0.022 + 0.003*dn + 0.005*dz;

  sig = 25.0*std::pow(bet,2) * sig;

 qrot11:   
  ponq = (u - ucr)/dcr;

  if (ponq > 700.0) {
    ponq = 700.0;
  }
  if (sig < 1.0) {
    sig = 1.0;
  }
  (*qr) = 1.0/(1.0 + std::exp(ponq)) * (sig - 1.0) + 1.0;

  if ((*qr) < 1.0) {
    (*qr) = 1.0;
  }

  return;
}

void G4Abla::mglw(G4double a, G4double z, G4double *el)
{
  // MODEL DE LA GOUTTE LIQUIDE DE C. F. WEIZSACKER.
  // USUALLY AN OBSOLETE OPTION

  G4int a1 = 0, z1 = 0;
  G4double xv = 0.0, xs = 0.0, xc = 0.0, xa = 0.0;                                   

  a1 = idnint(a);
  z1 = idnint(z);

  if ((a <= 0.01) || (z < 0.01)) {
    (*el) = 1.0e38;
  }
  else {
    xv = -15.56*a;
    xs = 17.23*std::pow(a,(2.0/3.0));

    if (a > 1.0) {
      xc = 0.7*z*(z-1.0)*std::pow((a-1.0),(-1.e0/3.e0));
    }
    else {
      xc = 0.0;
    }
  }

  xa = 23.6*(std::pow((a-2.0*z),2)/a);
  (*el) = xv+xs+xc+xa;
  return;	
}

void G4Abla::mglms(G4double a, G4double z, G4int refopt4, G4double *el)
{
  // USING FUNCTION EFLMAC(IA,IZ,0)                                    
  // 
  // REFOPT4 = 0 : WITHOUT MICROSCOPIC CORRECTIONS                     
  // REFOPT4 = 1 : WITH SHELL CORRECTION                               
  // REFOPT4 = 2 : WITH PAIRING CORRECTION                             
  // REFOPT4 = 3 : WITH SHELL- AND PAIRING CORRECTION                  

  //   1839	C-----------------------------------------------------------------------
  //   1840	C     A1       LOCAL    MASS NUMBER (INTEGER VARIABLE OF A)             
  //   1841	C     Z1       LOCAL    NUCLEAR CHARGE (INTEGER VARIABLE OF Z)          
  //   1842	C     REFOPT4           OPTION, SPECIFYING THE MASS FORMULA (SEE ABOVE) 
  //   1843	C     A                 MASS NUMBER                                     
  //   1844	C     Z                 NUCLEAR CHARGE                                  
  //   1845	C     DEL               PAIRING CORRECTION                              
  //   1846	C     EL                BINDING ENERGY                                  
  //   1847	C     ECNZ( , )         TABLE OF SHELL CORRECTIONS                      
  //   1848	C-----------------------------------------------------------------------
  //   1849	C                                                                       
  G4int a1 = idnint(a);
  G4int z1 = idnint(z);

  if ( (a1 <= 0) || (z1 <= 0) || ((a1-z1) <= 0) )  { //then 
    // modif pour récupérer une masse p et n correcte:
    (*el) = 0.0;
    return;
    //    goto mglms50;
  }
  else {
    // binding energy incl. pairing contr. is calculated from                
    // function eflmac                                                       
    assert(a1 != 0);
    (*el) = eflmac(a1,z1,0,refopt4);
    // assert(isnan((*el)) == false);
    if (refopt4 > 0) {
      if (refopt4 != 2) {
	(*el) = (*el) + ec2sub->ecnz[a1-z1][z1];
	//(*el) = (*el) + ec2sub->ecnz[z1][a1-z1];
	//assert(isnan((*el)) == false);
      }
    }
  }
  return;
}

G4double G4Abla::spdef(G4int a, G4int z, G4int optxfis)
{

  // INPUT:  A,Z,OPTXFIS MASS AND CHARGE OF A NUCLEUS,                     
  // OPTION FOR FISSILITY                                          
  // OUTPUT: SPDEF                                                         

  // ALPHA2 SADDLE POINT DEF. COHEN&SWIATECKI ANN.PHYS. 22 (1963) 406      
  // RANGING FROM FISSILITY X=0.30 TO X=1.00 IN STEPS OF 0.02              

  G4int index = 0;
  G4double x = 0.0, v = 0.0, dx = 0.0;

  const G4int alpha2Size = 37;
  // The value 0.0 at alpha2[0] added by PK.
  G4double alpha2[alpha2Size] = {0.0, 2.5464e0, 2.4944e0, 2.4410e0, 2.3915e0, 2.3482e0,
				 2.3014e0, 2.2479e0, 2.1982e0, 2.1432e0, 2.0807e0, 2.0142e0, 1.9419e0,
				 1.8714e0, 1.8010e0, 1.7272e0, 1.6473e0, 1.5601e0, 1.4526e0, 1.3164e0,
				 1.1391e0, 0.9662e0, 0.8295e0, 0.7231e0, 0.6360e0, 0.5615e0, 0.4953e0,
				 0.4354e0, 0.3799e0, 0.3274e0, 0.2779e0, 0.2298e0, 0.1827e0, 0.1373e0,
				 0.0901e0, 0.0430e0, 0.0000e0};

  dx = 0.02;
  x  = fissility(a,z,optxfis);

  if (x > 1.0) {
    x = 1.0;
  }

  if (x < 0.0) {
    x = 0.0;
  }

  v  = (x - 0.3)/dx + 1.0;
  index = idnint(v);

  if (index < 1) {
    return(alpha2[1]); // alpha2[0] -> alpha2[1]
  }

  if (index == 36) { //then // :::FIXME:: Possible off-by-one bug...                                            
    return(alpha2[36]); 
  }
  else {
    return(alpha2[index] + (alpha2[index+1] - alpha2[index]) / dx * ( x - (0.3e0 + dx*(index-1)))); //:::FIXME::: Possible off-by-one
  }                                                       

  return alpha2[0]; // The algorithm is not supposed to reach this point.
}

G4double G4Abla::fissility(int a,int z, int optxfis)
{
  // CALCULATION OF FISSILITY PARAMETER                                 
  // 
  // INPUT: A,Z INTEGER MASS & CHARGE OF NUCLEUS                        
  // OPTXFIS = 0 : MYERS, SWIATECKI                              
  //           1 : DAHLINGER                                     
  //           2 : ANDREYEV                                      

  G4double aa = 0.0, zz = 0.0, i = 0.0;
  G4double fissilityResult = 0.0;

  aa = double(a);
  zz = double(z);
  i  = double(a-2*z) / aa;

  // myers & swiatecki droplet modell                        
  if (optxfis == 0) { //then                                            
    fissilityResult = std::pow(zz,2) / aa /50.8830e0 / (1.0e0 - 1.7826e0 * std::pow(i,2));
  }

  if (optxfis == 1) {
    // dahlinger fit:                                          
    fissilityResult = std::pow(zz,2) / aa * std::pow((49.22e0*(1.e0 - 0.3803e0*std::pow(i,2) - 20.489e0*std::pow(i,4))),(-1));
  }

  if (optxfis == 2) {
    // dubna fit:                                              
    fissilityResult = std::pow(zz,2) / aa  /(48.e0*(1.e0 - 17.22e0*std::pow(i,4)));
  }

  return fissilityResult;
}

void G4Abla::evapora(G4double zprf, G4double aprf, G4double ee, G4double jprf, 
		     G4double *zf_par, G4double *af_par, G4double *mtota_par,
		     G4double *pleva_par, G4double *pxeva_par, G4double *pyeva_par,
		     G4int *ff_par, G4int *inttype_par, G4int *inum_par)
{
  G4double zf = (*zf_par);
  G4double af = (*af_par);
  G4double mtota = (*mtota_par);
  G4double pleva = (*pleva_par);
  G4double pxeva = (*pxeva_par);
  G4double pyeva = (*pyeva_par);
  G4int ff = (*ff_par);
  G4int inttype = (*inttype_par);
  G4int inum = (*inum_par);

  //    533	C                                                                       
  //    534	C     INPUT:                                                            
  //    535	C                                                                       
  //    536	C     ZPRF, APRF, EE(EE IS MODIFIED!), JPRF                             
  //    537	C                                                                       
  //    538	C     PROJECTILE AND TARGET PARAMETERS + CROSS SECTIONS                 
  //    539	C     COMMON /ABRAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,       
  //    540	C                       R_0,R_P,R_T, IMAX,IRNDM,PI,                     
  //    541	C                       BFPRO,SNPRO,SPPRO,SHELL                         
  //    542	C                                                                       
  //    543	C     AP,ZP,AT,ZT   - PROJECTILE AND TARGET MASSES                      
  //    544	C     EAP,BETA      - BEAM ENERGY PER NUCLEON, V/C                      
  //    545	C     BMAXNUC       - MAX. IMPACT PARAMETER FOR NUCL. REAC.             
  //    546	C     CRTOT,CRNUC   - TOTAL AND NUCLEAR REACTION CROSS SECTION          
  //    547	C     R_0,R_P,R_T,  - RADIUS PARAMETER, PROJECTILE+ TARGET RADII        
  //    548	C     IMAX,IRNDM,PI - MAXIMUM NUMBER OF EVENTS, DUMMY, 3.141...         
  //    549	C     BFPRO         - FISSION BARRIER OF THE PROJECTILE                 
  //    550	C     SNPRO         - NEUTRON SEPARATION ENERGY OF THE PROJECTILE       
  //    551	C     SPPRO         - PROTON    "           "   "    "   "              
  //    552	C     SHELL         - GROUND STATE SHELL CORRECTION                     
  //    553	C                                                                       
  //    554	C---------------------------------------------------------------------  
  //    555	C     FISSION BARRIERS                                                  
  //    556	C     COMMON /FB/     EFA                                               
  //    557	C     EFA    - ARRAY OF FISSION BARRIERS                                
  //    558	C---------------------------------------------------------------------  
  //    559	C     OUTPUT:                                                           
  //    560	C              ZF, AF, MTOTA, PLEVA, PTEVA, FF, INTTYPE, INUM           
  //    561	C                                                                       
  //    562	C     ZF,AF - CHARGE AND MASS OF FINAL FRAGMENT AFTER EVAPORATION       
  //    563	C     MTOTA _ NUMBER OF EVAPORATED ALPHAS                               
  //    564	C     PLEVA,PXEVA,PYEVA - MOMENTUM RECOIL BY EVAPORATION               
  //    565	C     INTTYPE - TYPE OF REACTION 0/1 NUCLEAR OR ELECTROMAGNETIC         
  //    566	C     FF      - 0/1 NO FISSION / FISSION EVENT                          
  //    567	C     INUM    - EVENTNUMBER                                             
  //    568	C   ____________________________________________________________________
  //    569	C  /                                                                    
  //    570	C  /  CALCUL DE LA MASSE ET CHARGE FINALES D'UNE CHAINE D'EVAPORATION   
  //    571	C  /                                                                    
  //    572	C  /  PROCEDURE FOR CALCULATING THE FINAL MASS AND CHARGE VALUES OF A   
  //    573	C  /  SPECIFIC EVAPORATION CHAIN, STARTING POINT DEFINED BY (APRF, ZPRF,
  //    574	C  /  EE)                                                               
  //    575	C  /  On ajoute les 3 composantes de l'impulsion (PXEVA,PYEVA,PLEVA)
  //    576	C  /    (actuellement PTEVA n'est pas correct; mauvaise norme...)                                               
  //    577	C  /____________________________________________________________________
  //    578	C                                                                       
  //    612	C                                                                       
  //    613	C-----------------------------------------------------------------------
  //    614	C     IRNDM             DUMMY ARGUMENT FOR RANDOM-NUMBER FUNCTION       
  //    615	C     SORTIE   LOCAL    HELP VARIABLE TO END THE EVAPORATION CHAIN      
  //    616	C     ZF                NUCLEAR CHARGE OF THE FRAGMENT                  
  //    617	C     ZPRF              NUCLEAR CHARGE OF THE PREFRAGMENT               
  //    618	C     AF                MASS NUMBER OF THE FRAGMENT                     
  //    619	C     APRF              MASS NUMBER OF THE PREFRAGMENT                  
  //    620	C     EPSILN            ENERGY BURNED IN EACH EVAPORATION STEP          
  //    621	C     MALPHA   LOCAL    MASS CONTRIBUTION TO MTOTA IN EACH EVAPORATION  
  //    622	C                        STEP                                           
  //    623	C     EE                EXCITATION ENERGY (VARIABLE)                    
  //    624	C     PROBP             PROTON EMISSION PROBABILITY                     
  //    625	C     PROBN             NEUTRON EMISSION PROBABILITY                    
  //    626	C     PROBA             ALPHA-PARTICLE EMISSION PROBABILITY             
  //    627	C     PTOTL             TOTAL EMISSION PROBABILITY                      
  //    628	C     E                 LOWEST PARTICLE-THRESHOLD ENERGY                
  //    629	C     SN                NEUTRON SEPARATION ENERGY                       
  //    630	C     SBP               PROTON SEPARATION ENERGY PLUS EFFECTIVE COULOMB 
  //    631	C                        BARRIER                                        
  //    632	C     SBA               ALPHA-PARTICLE SEPARATION ENERGY PLUS EFFECTIVE 
  //    633	C                        COULOMB BARRIER                                
  //    634	C     BP                EFFECTIVE PROTON COULOMB BARRIER                
  //    635	C     BA                EFFECTIVE ALPHA COULOMB BARRIER                 
  //    636	C     MTOTA             TOTAL MASS OF THE EVAPORATED ALPHA PARTICLES    
  //    637	C     X                 UNIFORM RANDOM NUMBER FOR NUCLEAR CHARGE        
  //    638	C     AMOINS   LOCAL    MASS NUMBER OF EVAPORATED PARTICLE              
  //    639	C     ZMOINS   LOCAL    NUCLEAR CHARGE OF EVAPORATED PARTICLE           
  //    640	C     ECP               KINETIC ENERGY OF PROTON WITHOUT COULOMB        
  //    641	C                        REPULSION                                      
  //    642	C     ECN               KINETIC ENERGY OF NEUTRON                       
  //    643	C     ECA               KINETIC ENERGY OF ALPHA PARTICLE WITHOUT COULOMB
  //    644	C                        REPULSION                                      
  //    645	C     PLEVA             TRANSVERSAL RECOIL MOMENTUM OF EVAPORATION      
  //    646	C     PTEVA             LONGITUDINAL RECOIL MOMENTUM OF EVAPORATION     
  //    647	C     FF                FISSION FLAG                                    
  //    648	C     INTTYPE           INTERACTION TYPE FLAG                           
  //    649	C     RNDX              RECOIL MOMENTUM IN X-DIRECTION IN A SINGLE STEP 
  //    650	C     RNDY              RECOIL MOMENTUM IN Y-DIRECTION IN A SINGLE STEP 
  //    651	C     RNDZ              RECOIL MOMENTUM IN Z-DIRECTION IN A SINGLE STEP 
  //    652	C     RNDN              NORMALIZATION OF RECOIL MOMENTUM FOR EACH STEP  
  //    653	C-----------------------------------------------------------------------
  //    654	C                                                                       
  //    655	      SAVE                                                              
  // SAVE -> static
	
  static G4int sortie = 0;                            
  static G4double epsiln = 0.0, probp = 0.0, probn = 0.0, proba = 0.0, ptotl = 0.0, e = 0.0;  
  static G4double sn = 0.0, sbp = 0.0, sba = 0.0, x = 0.0, amoins = 0.0, zmoins = 0.0;
  G4double ecn = 0.0, ecp = 0.0,eca = 0.0, bp = 0.0, ba = 0.0;         
  static G4double pteva = 0.0;                       

  static G4int itest = 0;
  static G4double probf = 0.0;

  static G4int k = 0, j = 0, il = 0;

  static G4double ctet1 = 0.0, stet1 = 0.0, phi1 = 0.0;
  static G4double sbfis = 0.0, rnd = 0.0;
  static G4double selmax = 0.0;
  static G4double segs = 0.0;
  static G4double ef = 0.0;
  static G4int irndm = 0;

  static G4double pc = 0.0, malpha = 0.0;

  zf = zprf;
  af = aprf;
  pleva = 0.0;
  pteva = 0.0;
  pxeva = 0.0;
  pyeva = 0.0;

  sortie = 0;
  ff = 0;

  itest = 0;
  if (itest == 1) {
    G4cout << "***************************" << G4endl;
  }

 evapora10:

  if (itest == 1) {
    G4cout <<"------zf,af,ee------" << idnint(zf) << "," << idnint(af) << "," << ee << G4endl;
  }

  // calculation of the probabilities for the different decay channels     
  // plus separation energies and kinetic energies of the particles        
  direct(zf,af,ee,jprf,&probp,&probn,&proba,&probf,&ptotl,
	 &sn,&sbp,&sba,&ecn,&ecp,&eca,&bp,&ba,inttype,inum,itest); //:::FIXME::: Call
  // assert(isnan(proba) == false);
  // assert(isnan(probp) == false);  
  // assert(isnan(probn) == false);
  // assert(isnan(probf) == false);  
  assert((eca+ba) >= 0);
  assert((ecp+bp) >= 0);
  // assert(isnan(ecp) == false);
  // assert(isnan(ecn) == false);
  // assert(isnan(bp) == false);
  // assert(isnan(ba) == false);
  k = idnint(zf);
  j = idnint(af-zf);

  // now ef is calculated from efa that depends on the subroutine
  // barfit which takes into account the modification on the ang. mom.
  // jb mvr 6-aug-1999
  // note *** shell correction! (ecgnz)  jb mvr 20-7-1999
  il  = idnint(jprf);
  barfit(k,k+j,il,&sbfis,&segs,&selmax);
  // assert(isnan(sbfis) == false);
  
  if ((fiss->optshp == 1) || (fiss->optshp == 3)) { //then                     
    //    fb->efa[k][j] = G4double(sbfis) -  ecld->ecgnz[j][k];
    fb->efa[j][k] = G4double(sbfis) -  ecld->ecgnz[j][k];
  }
  else {
    fb->efa[j][k] = G4double(sbfis);
    //  fb->efa[j][k] = G4double(sbfis);
  } //end if 
  ef = fb->efa[j][k];
  //  ef = fb->efa[j][k];
  // assert(isnan(fb->efa[j][k]) == false);
  // here the final steps of the evaporation are calculated                
  if ((sortie == 1) || (ptotl == 0.e0)) {
    e = dmin1(sn,sbp,sba);
    if (e > 1.0e30) {
      if(verboseLevel > 2) {
	G4cout <<"erreur a la sortie evapora,e>1.e30,af=" << af <<" zf=" << zf << G4endl;
      }
    }
    if (zf <= 6.0) {
      goto evapora100;
    }
    if (e < 0.0) {
      if (sn == e) {
	af = af - 1.e0;
      }
      else if (sbp == e) {
	af = af - 1.0;
	zf = zf - 1.0;
      }
      else if (sba == e) {
	af = af - 4.0;
	zf = zf - 2.0;
      }
      if (af < 2.5) {
	goto evapora100;
      }
      goto evapora10;
    }
    goto evapora100;
  }
  irndm = irndm + 1;

  // here the normal evaporation cascade starts                            

  // random number for the evaporation                                     
  //  x = double(Rndm(irndm))*ptotl;
  x = double(haz(1))*ptotl;

//   G4cout <<"proba = " << proba << G4endl;
//   G4cout <<"probp = " << probp << G4endl;
//   G4cout <<"probn = " << probn << G4endl;
//   G4cout <<"probf = " << probf << G4endl;

  itest = 0;
  if (x < proba) {
    // alpha evaporation                                                     
    if (itest == 1) {
      G4cout <<"< alpha evaporation >" << G4endl;
    }
    amoins = 4.0;
    zmoins = 2.0;
    epsiln = sba + eca;
    assert((std::pow((1.0 + (eca+ba)/3.72834e3),2) - 1.0) >= 0);
    pc = std::sqrt(std::pow((1.0 + (eca+ba)/3.72834e3),2) - 1.0) * 3.72834e3;
    // assert(isnan(pc) == false);
    malpha = 4.0;

    // volant:
    volant->iv = volant->iv + 1;
    volant->acv[volant->iv] = 4.;
    volant->zpcv[volant->iv] = 2.;
    volant->pcv[volant->iv] = pc;
  }
  else if (x < proba+probp) {
    // proton evaporation                                                    
    if (itest == 1) {
      G4cout <<"< proton evaporation >" << G4endl;
    }
    amoins = 1.0;
    zmoins = 1.0;
    epsiln = sbp + ecp;
    assert((std::pow((1.0 + (ecp + bp)/9.3827e2),2) - 1.0) >= 0);
    pc = std::sqrt(std::pow((1.0 + (ecp + bp)/9.3827e2),2) - 1.0) * 9.3827e2;
    // assert(isnan(pc) == false);
    malpha = 0.0;
    // volant:
    volant->iv = volant->iv + 1;
    volant->acv[volant->iv] = 1.0;
    volant->zpcv[volant->iv] = 1.;
    volant->pcv[volant->iv] = pc;
  }
  else if (x < proba+probp+probn) {
    // neutron evaporation                                                   
    if (itest == 1) {
      G4cout <<"< neutron evaporation >" << G4endl;
    }
    amoins = 1.0;
    zmoins = 0.0;
    epsiln = sn + ecn;
    assert((std::pow((1.0 + (ecn)/9.3956e2),2) - 1.0) >= 0);
    pc = std::sqrt(std::pow((1.0 + (ecn)/9.3956e2),2) - 1.0) * 9.3956e2;
    // assert(isnan(pc) == false);
    malpha = 0.0;
  
    // volant:
    volant->iv = volant->iv + 1;
    volant->acv[volant->iv] = 1.;
    volant->zpcv[volant->iv] = 0.;
    volant->pcv[volant->iv] = pc;
  }
  else {
    // fission                                                               
    // in case of fission-events the fragment nucleus is the mother nucleus  
    // before fission occurs with excitation energy above the fis.- barrier. 
    // fission fragment mass distribution is calulated in subroutine fisdis  
    if (itest == 1) {
      G4cout <<"< fission >" << G4endl;
    }
    amoins = 0.0;
    zmoins = 0.0;
    epsiln = ef;

    malpha = 0.0;
    pc = 0.0;
    ff = 1;
    // ff = 0; // For testing, allows to disable fission!
  }

  if (itest == 1) {
    G4cout <<"sn,sbp,sba,ef" << sn << "," << sbp << "," << sba <<"," << ef << G4endl;
    G4cout <<"probn,probp,proba,probf,ptotl " <<","<< probn <<","<< probp <<","<< proba <<","<< probf <<","<< ptotl << G4endl;
  }

  // calculation of the daughter nucleus                                   
  af = af - amoins;
  zf = zf - zmoins;
  ee = ee - epsiln;
  if (ee <= 0.01) {
    ee = 0.01;
  }
  mtota = mtota + malpha;

  if(ff == 0) {
    standardRandom(&rnd,&(hazard->igraine[8]));
    ctet1 = 2.0*rnd - 1.0;
    standardRandom(&rnd,&(hazard->igraine[4]));
    phi1 = rnd*2.0*3.141592654;
    stet1 = std::sqrt(1.0 - std::pow(ctet1,2));
    // assert(isnan(stet1) == false);
    volant->xcv[volant->iv] = stet1*std::cos(phi1);
    volant->ycv[volant->iv] = stet1*std::sin(phi1);
    volant->zcv[volant->iv] = ctet1;
    pxeva = pxeva - pc * volant->xcv[volant->iv];
    pyeva = pyeva - pc * volant->ycv[volant->iv];
    pleva = pleva - pc * ctet1;
    // assert(isnan(pleva) == false);
  }

  // condition for end of evaporation                                   
  if ((af < 2.5) || (ff == 1)) {
    goto evapora100;
  }
  goto evapora10;

 evapora100:
  (*zf_par) = zf;
  (*af_par) = af;
  (*mtota_par) = mtota;
  (*pleva_par) = pleva;
  (*pxeva_par) = pxeva;
  (*pyeva_par) = pyeva;
  (*ff_par) = ff;
  (*inttype_par) = inttype;                                          
  (*inum_par) = inum;

  return;
}

void G4Abla::direct(G4double zprf, G4double a, G4double ee, G4double jprf, 
		    G4double *probp_par, G4double *probn_par, G4double *proba_par, 
		    G4double *probf_par, G4double *ptotl_par, G4double *sn_par,
		    G4double *sbp_par, G4double *sba_par, G4double *ecn_par, 
		    G4double *ecp_par,G4double *eca_par, G4double *bp_par,
		    G4double *ba_par, G4int inttype, G4int inum, G4int itest)
{
  G4int dummy0 = 0;
  
  G4double probp = (*probp_par);
  G4double probn = (*probn_par);
  G4double proba = (*proba_par);
  G4double probf = (*probf_par);
  G4double ptotl = (*ptotl_par);
  G4double sn = (*sn_par);
  G4double sbp = (*sbp_par);
  G4double sba = (*sba_par);
  G4double ecn = (*ecn_par);
  G4double ecp = (*ecp_par);
  G4double eca = (*eca_par);
  G4double bp = (*bp_par);
  G4double ba = (*ba_par);

  // CALCULATION OF PARTICLE-EMISSION PROBABILITIES & FISSION     / 
  // BASED ON THE SIMPLIFIED FORMULAS FOR THE DECAY WIDTH BY      / 
  // MORETTO, ROCHESTER MEETING TO AVOID COMPUTING TIME           / 
  // INTENSIVE INTEGRATION OF THE LEVEL DENSITIES                 / 
  // USES EFFECTIVE COULOMB BARRIERS AND AN AVERAGE KINETIC ENERGY/ 
  // OF THE EVAPORATED PARTICLES                                  / 
  // COLLECTIVE ENHANCMENT OF THE LEVEL DENSITY IS INCLUDED       / 
  // DYNAMICAL HINDRANCE OF FISSION IS INCLUDED BY A STEP FUNCTION/ 
  // APPROXIMATION. SEE A.R. JUNGHANS DIPLOMA THESIS              / 
  // SHELL AND PAIRING STRUCTURES IN THE LEVEL DENSITY IS INCLUDED/ 

  // INPUT:                                                            
  // ZPRF,A,EE  CHARGE, MASS, EXCITATION ENERGY OF COMPOUND     
  // NUCLEUS                                         
  // JPRF       ROOT-MEAN-SQUARED ANGULAR MOMENTUM                           

  // DEFORMATIONS AND G.S. SHELL EFFECTS                               
  // COMMON /ECLD/   ECGNZ,ECFNZ,VGSLD,ALPHA                           

  // ECGNZ - GROUND STATE SHELL CORR. FRLDM FOR A SPHERICAL G.S.       
  // ECFNZ - SHELL CORRECTION FOR THE SADDLE POINT (NOW: == 0)         
  // VGSLD - DIFFERENCE BETWEEN DEFORMED G.S. AND LDM VALUE            
  // ALPHA - ALPHA GROUND STATE DEFORMATION (THIS IS NOT BETA2!)       
  // BETA2 = SQRT((4PI)/5) * ALPHA                             

  //OPTIONS AND PARAMETERS FOR FISSION CHANNEL                        
  //COMMON /FISS/    AKAP,BET,HOMEGA,KOEFF,IFIS,                       
  //                 OPTSHP,OPTXFIS,OPTLES,OPTCOL               
  //
  // AKAP   - HBAR**2/(2* MN * R_0**2) = 10 MEV, R_0 = 1.4 FM          
  // BET    - REDUCED NUCLEAR FRICTION COEFFICIENT IN (10**21 S**-1)   
  // HOMEGA - CURVATURE OF THE FISSION BARRIER = 1 MEV                 
  // KOEFF  - COEFFICIENT FOR THE LD FISSION BARRIER == 1.0            
  // IFIS   - 0/1 FISSION CHANNEL OFF/ON                               
  // OPTSHP - INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
  //          = 0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY        
  //          = 1 SHELL ,  NO PAIRING                                    
  //          = 2 PAIRING, NO SHELL                                      
  //          = 3 SHELL AND PAIRING                                      
  // OPTCOL - 0/1 COLLECTIVE ENHANCEMENT SWITCHED ON/OFF               
  // OPTXFIS- 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV         
  //                FISSILITY PARAMETER.                                     
  // OPTLES - CONSTANT TEMPERATURE LEVEL DENSITY FOR A,Z > TH-224      
  // OPTCOL - 0/1 COLLECTIVE ENHANCEMENT OFF/ON                        

  // LEVEL DENSITY PARAMETERS                                          
  // COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
  //                 AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
  //                            LEVEL DENSITY PARAMETER                                
  // OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
  //           RECOMMENDED IS OPTAFAN = 0                              

  // FISSION BARRIERS                                                  
  // COMMON /FB/     EFA                                               
  // EFA    - ARRAY OF FISSION BARRIERS                                


  // OUTPUT: PROBN,PROBP,PROBA,PROBF,PTOTL:                            
  // - EMISSION PROBABILITIES FOR N EUTRON, P ROTON,  A LPHA     
  // PARTICLES, F ISSION AND NORMALISATION                     
  // SN,SBP,SBA: SEPARATION ENERGIES N P A                     
  // INCLUDING EFFECTIVE BARRIERS                              
  // ECN,ECP,ECA,BP,BA                                         
  // - AVERAGE KINETIC ENERGIES (2*T) AND EFFECTIVE BARRIERS     

  static G4double bk = 0.0;
  static G4int afp = 0;
  static G4double at = 0.0;
  static G4double bs = 0.0;
  static G4double bshell = 0.0;
  static G4double cf = 0.0;
  static G4double dconst = 0.0;
  static G4double defbet = 0.0;
  static G4double denomi = 0.0;
  static G4double densa = 0.0;
  static G4double densf = 0.0;
  static G4double densg = 0.0;
  static G4double densn = 0.0;
  static G4double densp = 0.0;
  static G4double edyn = 0.0;
  static G4double eer = 0.0;
  static G4double ef = 0.0;
  static G4double ft = 0.0;
  static G4double ga = 0.0;
  static G4double gf = 0.0;
  static G4double gn = 0.0;
  static G4double gngf = 0.0;
  static G4double gp = 0.0;
  static G4double gsum = 0.0;
  static G4double hbar = 6.582122e-22; // = 0.0;
  static G4double iflag = 0.0;
  static G4int il = 0;
  static G4int imaxwell = 0;
  static G4int in = 0;
  static G4int iz = 0;
  static G4int j = 0;
  static G4int k = 0;
  static G4double ma1z = 0.0;
  static G4double ma1z1 = 0.0;
  static G4double ma4z2 = 0.0;
  static G4double maz = 0.0;
  static G4double nprf = 0.0;
  static G4double nt = 0.0;
  static G4double parc = 0.0;
  static G4double pi = 3.14159265;
  static G4double pt = 0.0;
  static G4double ra = 0.0;
  static G4double rat = 0.0;
  static G4double refmod = 0.0;
  static G4double rf = 0.0;
  static G4double rn = 0.0;
  static G4double rnd = 0.0;
  static G4double rnt = 0.0;
  static G4double rp = 0.0;
  static G4double rpt = 0.0;
  static G4double sa = 0.0;
  static G4double sbf = 0.0;
  static G4double sbfis = 0.0;
  static G4double segs = 0.0;
  static G4double selmax = 0.0;
  static G4double sp = 0.0;
  static G4double tauc = 0.0;
  static G4double tconst = 0.0;
  static G4double temp = 0.0;
  static G4double ts1 = 0.0;
  static G4double tsum = 0.0;
  static G4double wf = 0.0;
  static G4double wfex = 0.0;
  static G4double xx = 0.0;
  static G4double y = 0.0;

  imaxwell = 1;
  inttype = 0;
  
  // limiting of excitation energy where fission occurs                    
  // Note, this is not the dynamical hindrance (see end of routine)      
  edyn = 1000.0;

  // no limit if statistical model is calculated.                         
  if (fiss->bet <= 1.0e-16) {
    edyn = 10000.0;
  }

  // just a change of name until the end of this subroutine                
  eer = ee;
  if (inum == 1) {
    ilast = 1;
  }

  // calculation of masses                                           
  // refmod = 1 ==> myers,swiatecki model                              
  // refmod = 0 ==> weizsaecker model                                  
  refmod = 1;  // Default = 1

  if (refmod == 1) {
    mglms(a,zprf,fiss->optshp,&maz);
    mglms(a-1.0,zprf,fiss->optshp,&ma1z);
    mglms(a-1.0,zprf-1.0,fiss->optshp,&ma1z1);
    mglms(a-4.0,zprf-2.0,fiss->optshp,&ma4z2);
  }
  else {
    mglw(a,zprf,&maz);
    mglw(a-1.0,zprf,&ma1z);
    mglw(a-1.0,zprf-1.0,&ma1z1);
    mglw(a-4.0,zprf-2.0,&ma4z2);
  }
  // assert(isnan(maz) == false);
  // assert(isnan(ma1z) == false);
  // assert(isnan(ma1z1) == false);
  // assert(isnan(ma4z2) == false);
  
  // separation energies and effective barriers                     
  sn = ma1z - maz;
  sp = ma1z1 - maz;
  sa = ma4z2 - maz - 28.29688;
  if (zprf < 1.0e0) {
    sbp = 1.0e75;
    goto direct30;
  }

  // parameterisation gaimard:
  // bp = 1.44*(zprf-1.d0)/(1.22*std::pow((a - 1.0),(1.0/3.0))+5.6)     
  // parameterisation khs (12-99)
  bp = 1.44*(zprf - 1.0)/(2.1*std::pow((a - 1.0),(1.0/3.0)) + 0.0);

  sbp = sp + bp;
  // assert(isnan(sbp) == false);
  // assert(isinf(sbp) == false);
  if (a-4.0 <= 0.0) {
    sba = 1.0e+75;
    goto direct30;
  }

  // new effective barrier for alpha evaporation d=6.1: khs          
  // ba = 2.88d0*(zprf-2.d0)/(1.22d0*(a-4.d0)**(1.d0/3.d0)+6.1d0)
  // parametrisation khs (12-99)
  ba = 2.88*(zprf - 2.0)/(2.2*std::pow((a - 4.0),(1.0/3.0)) + 0.0);

  sba = sa + ba;
  // assert(isnan(sba) == false);
  // assert(isinf(sba) == false);
 direct30:

  // calculation of surface and curvature integrals needed to      
  // to calculate the level density parameter (in densniv)         
  if (fiss->ifis > 0) {
    k = idnint(zprf);
    j = idnint(a - zprf);

    // now ef is calculated from efa that depends on the subroutine
    // barfit which takes into account the modification on the ang. mom.
    // jb mvr 6-aug-1999
    // note *** shell correction! (ecgnz)  jb mvr 20-7-1999
    il = idnint(jprf);
    barfit(k,k+j,il,&sbfis,&segs,&selmax);
    // assert(isnan(sbfis) == false);
    if ((fiss->optshp == 1) || (fiss->optshp == 3)) {
      //      fb->efa[k][j] = G4double(sbfis) -  ecld->ecgnz[j][k];
      //      fb->efa[j][k] = G4double(sbfis) -  ecld->ecgnz[j][k];
      fb->efa[j][k] = double(sbfis) -  ecld->ecgnz[j][k];
    } 
    else {
      //      fb->efa[k][j] = G4double(sbfis);
      fb->efa[j][k] = double(sbfis);
    }
    //    ef = fb->efa[k][j];
    ef = fb->efa[j][k];

    // to avoid negative values for impossible nuclei                        
    // the fission barrier is set to zero if smaller than zero.              
    //     if (fb->efa[k][j] < 0.0) {
    //       fb->efa[k][j] = 0.0;
    //     }
    if (fb->efa[j][k] < 0.0) {
      if(verboseLevel > 2) {
	G4cout <<"Setting fission barrier to 0" << G4endl;
      }
      fb->efa[j][k] = 0.0;
    }
    // assert(isnan(fb->efa[j][k]) == false);
    
    // factor with jprf should be 0.0025d0 - 0.01d0 for                     
    // approximate influence of ang. momentum on bfis  a.j. 22.07.96        
    // 0.0 means no angular momentum                                       

    if (ef < 0.0) {
      ef = 0.0;
    }
    xx = fissility((k+j),k,fiss->optxfis);
    // assert(isnan(xx) == false);
    // assert(isinf(xx) == false);
    
    y = 1.00 - xx;
    if (y < 0.0) {
      y = 0.0;
    }
    if (y > 1.0) {
      y = 1.0;
    }
    bs = bipol(1,y);
    // assert(isnan(bs) == false);
    // assert(isinf(bs) == false);
    bk = bipol(2,y);
    // assert(isnan(bk) == false);
    // assert(isinf(bk) == false);
  }
  else {
    ef = 1.0e40;
    bs = 1.0;
    bk = 1.0;
  }
  sbf = ee - ef;

  afp = idnint(a);
  iz = idnint(zprf);
  in = afp - iz;
  bshell = ecld->ecfnz[in][iz];
  // assert(isnan(bshell) == false);

  // ld saddle point deformation                                          
  // here: beta2 = std::sqrt(5/(4pi)) * alpha2                                  

  // for the ground state def. 1.5d0 should be used                        
  // because this was just the factor to produce the                       
  // alpha-deformation table 1.5d0 should be used                          
  // a.r.j. 6.8.97                                                         
  defbet = 1.58533e0 * spdef(idnint(a),idnint(zprf),fiss->optxfis);
  // assert(isnan(defbet) == false);
  
  // level density and temperature at the saddle point                     
  //   G4cout <<"a = " << a << G4endl;
  //   G4cout <<"zprf = " << zprf << G4endl;
  //   G4cout <<"ee = " << ee << G4endl;
  //   G4cout <<"ef = " << ef << G4endl;
  //   G4cout <<"bshell = " << bshell << G4endl;
  //   G4cout <<"bs = " << bs << G4endl;
  //   G4cout <<"bk = " << bk << G4endl;
  //   G4cout <<"defbet = " << defbet << G4endl;
  densniv(a,zprf,ee,ef,&densf,bshell,bs,bk,&temp,int(fiss->optshp),int(fiss->optcol),defbet);
  //   G4cout <<"densf = " << densf << G4endl;
  //   G4cout <<"temp = " << temp << G4endl;
  // assert(isnan(densf) == false);
  // assert(isnan(temp) == false);
  //  assert(temp != 0);
  ft = temp;
  if (iz >= 2) {
    bshell = ecld->ecgnz[in][iz-1] - ecld->vgsld[in][iz-1];
    defbet = 1.5 * (ecld->alpha[in][iz-1]);

    // level density and temperature in the proton daughter                  
    densniv(a-1.0,zprf-1.0e0,ee,sbp,&densp, bshell,1.e0,1.e0,&temp,int(fiss->optshp),int(fiss->optcol),defbet);
    assert(temp >= 0);
    // assert(isnan(temp) == false);
    pt = temp;
    if (imaxwell == 1) {
      // valentina - random kinetic energy in a maxwelliam distribution
      // modif juin/2002 a.b. c.v. for light targets; limit on the energy
      // from the maxwell distribution.
      rpt = pt;
      ecp = 2.0 * pt;
      if(rpt <= 1.0e-3) {
	goto direct2914;
      }
      iflag = 0;
    direct1914:
      ecp = fmaxhaz(rpt);
      iflag = iflag + 1;
      if(iflag >= 10) {
	standardRandom(&rnd,&(hazard->igraine[5]));
	ecp=std::sqrt(rnd)*(eer-sbp);
	// assert(isnan(ecp) == false);
	goto direct2914;
      }
      if((ecp+sbp) > eer) {
	goto direct1914;
      }
    }
    else {
      ecp = 2.0 * pt;
    }

  direct2914:
    dummy0 = 0;
    //    G4cout <<""<<G4endl;
  }
  else {
    densp = 0.0;
    ecp = 0.0;
    pt = 0.0;
  }

  if (in >= 2) {
    bshell = ecld->ecgnz[in-1][iz] - ecld->vgsld[in-1][iz];
    defbet = 1.5e0 * (ecld->alpha[in-1][iz]);

    // level density and temperature in the neutron daughter                 
    densniv(a-1.0,zprf,ee,sn,&densn,bshell, 1.e0,1.e0,&temp,int(fiss->optshp),int(fiss->optcol),defbet);
    nt = temp;

    if (imaxwell == 1) {
      // valentina - random kinetic energy in a maxwelliam distribution
      // modif juin/2002 a.b. c.v. for light targets; limit on the energy
      // from the maxwell distribution.
      rnt = nt;
      ecn = 2.0 * nt;
      if(rnt <= 1.e-3) {
	goto direct2915;
      }

      iflag=0;

      ecn = fmaxhaz(rnt);
      iflag=iflag+1;
      if(iflag >= 10) {
	standardRandom(&rnd,&(hazard->igraine[6]));
	ecn = std::sqrt(rnd)*(eer-sn);
	// assert(isnan(ecn) == false);
	goto direct2915;
      }
      //       if((ecn+sn) > eer) {
      // 	goto direct1915;
      //       }
      //       else {
      // 	ecn = 2.e0 * nt;
      //       }
      if((ecn + sn) <= eer) {
	ecn = 2.0 * nt;
      }
    direct2915: 
      dummy0 = 0;
      //      G4cout <<"" <<G4endl;
    }
  }
  else {
    densn = 0.0;
    ecn = 0.0;
    nt = 0.0;
  }

  if ((in >= 3) && (iz >= 3)) {
    bshell = ecld->ecgnz[in-2][iz-2] - ecld->vgsld[in-2][iz-2];
    defbet = 1.5 * (ecld->alpha[in-2][iz-2]);

    // level density and temperature in the alpha daughter                   
    densniv(a-4.0,zprf-2.0e0,ee,sba,&densa,bshell,1.e0,1.e0,&temp,int(fiss->optshp),int(fiss->optcol),defbet);

    // valentina - random kinetic energy in a maxwelliam distribution
    at = temp;
    if (imaxwell == 1) {
      // modif juin/2002 a.b. c.v. for light targets; limit on the energy
      // from the maxwell distribution.
      rat = at;
      eca= 2.e0 * at;
      if(rat <= 1.e-3) {
	goto direct2916;
      }
      iflag=0;
    direct1916:
      eca = fmaxhaz(rat);
      iflag=iflag+1;
      if(iflag >= 10) {
	standardRandom(&rnd,&(hazard->igraine[7]));
	eca=std::sqrt(rnd)*(eer-sba);
	// assert(isnan(eca) == false);
	goto direct2916;
      }
      if((eca+sba) > eer) {
	goto direct1916;
      }
      else {
	eca = 2.0 * at;
      }
    direct2916:
      dummy0 = 0;
      //      G4cout <<"" << G4endl;
    }
    else {
      densa = 0.0;
      eca = 0.0;
      at = 0.0;
    }
  } // PK

  // special treatment for unbound nuclei                                                
  if (sn < 0.0) {
    probn = 1.0;
    probp = 0.0;
    proba = 0.0;
    probf = 0.0;
    goto direct70;
  }
  if (sbp < 0.0) {
    probp = 1.0;
    probn = 0.0;
    proba = 0.0;
    probf = 0.0;
    goto direct70;
  }

  if ((a < 50.e0) || (ee > edyn)) { // no fission if e*> edyn or mass < 50
    //    G4cout <<"densf = 0.0" << G4endl;
    densf = 0.e0;
  }

  bshell = ecld->ecgnz[in][iz] - ecld->vgsld[in][iz];
  defbet = 1.5e0 * (ecld->alpha[in][iz]);

  // compound nucleus level density                                        
  densniv(a,zprf,ee,0.0e0,&densg,bshell,1.e0,1.e0,&temp,int(fiss->optshp),int(fiss->optcol),defbet);
  // assert(isnan(densg) == false);
  // assert(isnan(temp) == false);
  
  if ( densg > 0.e0) {
    // calculation of the partial decay width                                
    // used for both the time scale and the evaporation decay width
    gp = (std::pow(a,(2.0/3.0))/fiss->akap)*densp/densg/pi*std::pow(pt,2);
    gn = (std::pow(a,(2.0/3.0))/fiss->akap)*densn/densg/pi*std::pow(nt,2);
    ga = (std::pow(a,(2.0/3.0))/fiss->akap)*densa/densg/pi*2.0*std::pow(at,2);
    gf = densf/densg/pi/2.0*ft;
    // assert(isnan(gf) == false);
    
    //     assert(isnan(gp) == false);
    //     assert(isnan(gn) == false);
    //     assert(isnan(ga) == false);
    //     assert(isnan(ft) == false);
    //    assert(ft != 0);
    //    assert(isnan(gf) == false);
    
    if(itest == 1) {
      G4cout <<"gn,gp,ga,gf " << gn <<","<< gp <<","<< ga <<","<< gf << G4endl;
    }
  }
  else {
    if(verboseLevel > 2) {
      G4cout <<"direct: densg <= 0.e0 " << a <<","<< zprf <<","<< ee << G4endl;
    }
  }

  gsum = ga + gp + gn;
  // assert(isinf(gsum) == false);
  // assert(isnan(gsum) == false);
  if (gsum > 0.0) {
    ts1  = hbar / gsum;
  }
  else {
    ts1  = 1.0e99;
  }

  if (inum > ilast) {  // new event means reset the time scale
    tsum = 0;
  }

  // calculate the relative probabilities for all decay channels        
  if (densf == 0.0) {
    if (densp == 0.0) {
      if (densn == 0.0) {
	if (densa == 0.0) {
	  // no reaction is possible                                               
	  probf = 0.0;
	  probp = 0.0;
	  probn = 0.0;
	  proba = 0.0;
	  goto direct70;
	}

	// alpha evaporation is the only open channel                            
	rf = 0.0;
	rp = 0.0;
	rn = 0.0;
	ra = 1.0;
	goto direct50;
      }

      // alpha emission and neutron emission                                   
      rf = 0.0;
      rp = 0.0;
      rn = 1.0;
      ra = densa*2.0/densn*std::pow((at/nt),2);
      goto direct50;
    }
    // alpha, proton and neutron emission                                    
    rf = 0.0;
    rp = 1.0;
    rn = densn/densp*std::pow((nt/pt),2);
    // assert(isnan(rn) == false);
    ra = densa*2.0/densp*std::pow((at/pt),2);
    // assert(isnan(ra) == false);
    goto direct50;
  }

  // here fission has taken place                                          
  rf = 1.0;

  // cramers and weidenmueller factors for the dynamical hindrances of     
  // fission                                                               
  if (fiss->bet <= 1.0e-16) {
    cf = 1.0;
    wf = 1.0;
  }
  else if (sbf > 0.0e0) {
    cf = cram(fiss->bet,fiss->homega);
    // if fission barrier ef=0.d0 then fission is the only possible      
    // channel. to avoid std::log(0) in function tau                          
    // a.j. 7/28/93                                                      
    if (ef <= 0.0) {
      rp = 0.0;
      rn = 0.0;
      ra = 0.0;
      goto direct50;
    }
    else {
      // transient time tau()                                                  
      tauc = tau(fiss->bet,fiss->homega,ef,ft);
      // assert(isnan(tauc) == false);
    }
    wfex = (tauc - tsum)/ts1;

    if (wfex < 0.0) {
      wf = 1.0;
    }
    else {
      wf = std::exp( -wfex);
    }
  }
  else {
    cf=1.0;
    wf=1.0;
  }

  if(verboseLevel > 2) {
    G4cout <<"tsum,wf,cf " << tsum <<","<< wf <<","<< cf << G4endl;
  }

  tsum = tsum + ts1;

  // change by g.k. and a.h. 5.9.95                                       
  tconst = 0.7;
  dconst = 12.0/std::sqrt(a);
  // assert(isnan(dconst) == false);
  nprf = a - zprf;

  if (fiss->optshp >= 2) { //then                                           
    parite(nprf,&parc);
    // assert(isnan(parc) == false);
    dconst = dconst*parc;
  }
  else {
    dconst= 0.0;
  }
  if ((ee <= 17.e0) && (fiss->optles == 1) && (iz >= 90) && (in >= 134)) { //then                              
    // constant changed to 5.0 accord to moretto & vandenbosch a.j. 19.3.96  
    gngf = std::pow(a,(2.0/3.0))*tconst/10.0*std::exp((ef-sn+dconst)/tconst);

    // if the excitation energy is so low that densn=0 ==> gn = 0           
    // fission remains the only channel.                                    
    // a. j. 10.1.94                                                        
    if (gn == 0.0) {
      rn = 0.0;
      rp = 0.0;
      ra = 0.0;
    }
    else {
      rn=gngf;
      // assert(isnan(rn) == false);
      rp=gngf*gp/gn;
      // assert(isnan(rp) == false);      
      ra=gngf*ga/gn;
      // assert(isnan(ra) == false);      
    }
  } else {
    // assert(isnan(cf) == false);
    // assert(isinf(gn) == false);
    // assert(isinf(gf) == false);
    // assert(isinf(cf) == false);        
    assert(gn > 0 || (gf != 0 && cf != 0));
    rn = gn/(gf*cf);
//     G4cout <<"rn = " << G4endl;
//     G4cout <<"gn = " << gn << " gf = " << gf << " cf = " << cf << G4endl;
    // assert(isnan(rn) == false);
    rp = gp/(gf*cf);
    // assert(isnan(rp) == false);    
    ra = ga/(gf*cf);
    // assert(isnan(ra) == false);    
  }
 direct50:
  // relative decay probabilities                                          
  // assert(isnan(ra) == false);
  // assert(isnan(rp) == false);
  // assert(isnan(rn) == false);
  // assert(isnan(rf) == false);
  
  denomi = rp+rn+ra+rf;
  // assert(isnan(denomi) == false);
  assert(denomi > 0);
  // decay probabilities after transient time
  probf = rf/denomi;
  // assert(isnan(probf) == false);
  probp = rp/denomi;
  // assert(isnan(probp) == false);
  probn = rn/denomi;
  // assert(isnan(probn) == false);
  proba = ra/denomi;
  // assert(isnan(proba) == false);
  // assert(isinf(proba) == false);
  
  // new treatment of grange-weidenmueller factor, 5.1.2000, khs !!!

  // decay probabilites with transient time included
  // assert(isnan(wf) == false);
  assert(std::fabs(probf) <= 1.0);
  probf = probf * wf;
  if(probf == 1.0) {
    probp = 0.0;
    probn = 0.0;
    proba = 0.0;
  }
  else {
    probp = probp * (wf + (1.e0-wf)/(1.e0-probf));
    probn = probn * (wf + (1.e0-wf)/(1.e0-probf));
    proba = proba * (wf + (1.e0-wf)/(1.e0-probf));
  }
 direct70:
  // assert(isnan(probp) == false);
  // assert(isnan(probn) == false);
  // assert(isnan(probf) == false);  
  // assert(isnan(proba) == false);
  ptotl = probp+probn+proba+probf;
  // assert(isnan(ptotl) == false);
  
  ee = eer;
  ilast = inum;

  // Return values:
  // assert(isnan(proba) == false);
  (*probp_par) = probp;
  (*probn_par) = probn;
  (*proba_par) = proba;
  (*probf_par) = probf;
  (*ptotl_par) = ptotl;
  (*sn_par) = sn;
  (*sbp_par) = sbp;
  (*sba_par) = sba;
  (*ecn_par) = ecn;
  (*ecp_par) = ecp;
  (*eca_par) = eca;
  (*bp_par) = bp;
  (*ba_par) = ba;
}

void G4Abla::densniv(G4double a, G4double z, G4double ee, G4double esous, G4double *dens, G4double bshell, G4double bs, G4double bk, 
		     G4double *temp, G4int optshp, G4int optcol, G4double defbet)
{
  //   1498	C                                                                       
  //   1499	C     INPUT:                                                            
  //   1500	C             A,EE,ESOUS,OPTSHP,BS,BK,BSHELL,DEFBET                     
  //   1501	C                                                                       
  //   1502	C     LEVEL DENSITY PARAMETERS                                          
  //   1503	C     COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
  //   1504	C     AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
  //   1505	C                LEVEL DENSITY PARAMETER                                
  //   1506	C     OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
  //   1507	C               RECOMMENDED IS OPTAFAN = 0                              
  //   1508	C---------------------------------------------------------------------  
  //   1509	C     OUTPUT: DENS,TEMP                                                 
  //   1510	C                                                                       
  //   1511	C   ____________________________________________________________________
  //   1512	C  /                                                                    
  //   1513	C  /  PROCEDURE FOR CALCULATING THE STATE DENSITY OF A COMPOUND NUCLEUS 
  //   1514	C  /____________________________________________________________________
  //   1515	C                                                                       
  //   1516	      INTEGER AFP,IZ,OPTSHP,OPTCOL,J,OPTAFAN                            
  //   1517	      REAL*8 A,EE,ESOUS,DENS,E,Y0,Y1,Y2,Y01,Y11,Y21,PA,BS,BK,TEMP       
  //   1518	C=====INSERTED BY KUDYAEV===============================================
  //   1519	      COMMON /ALD/ AV,AS,AK,OPTAFAN                                     
  //   1520	      REAL*8 ECR,ER,DELTAU,Z,DELTPP,PARA,PARZ,FE,HE,ECOR,ECOR1,Pi6      
  //   1521	      REAL*8 BSHELL,DELTA0,AV,AK,AS,PONNIV,PONFE,DEFBET,QR,SIG,FP       
  //   1522	C=======================================================================
  //   1523	C                                                                       
  //   1524	C                                                                       
  //   1525	C-----------------------------------------------------------------------
  //   1526	C     A                 MASS NUMBER OF THE DAUGHTER NUCLEUS             
  //   1527	C     EE                EXCITATION ENERGY OF THE MOTHER NUCLEUS         
  //   1528	C     ESOUS             SEPARATION ENERGY PLUS EFFECTIVE COULOMB BARRIER
  //   1529	C     DENS              STATE DENSITY OF DAUGHTER NUCLEUS AT EE-ESOUS-EC
  //   1530	C     BSHELL            SHELL CORRECTION                                
  //   1531	C     TEMP              NUCLEAR TEMPERATURE                             
  //   1532	C     E        LOCAL    EXCITATION ENERGY OF THE DAUGHTER NUCLEUS       
  //   1533	C     E1       LOCAL    HELP VARIABLE                                   
  //   1534	C     Y0,Y1,Y2,Y01,Y11,Y21                                              
  //   1535	C              LOCAL    HELP VARIABLES                                  
  //   1536	C     PA       LOCAL    STATE-DENSITY PARAMETER                         
  //   1537	C     EC                KINETIC ENERGY OF EMITTED PARTICLE WITHOUT      
  //   1538	C                        COULOMB REPULSION                              
  //   1539	C     IDEN              FAKTOR FOR SUBSTRACTING KINETIC ENERGY IDEN*TEMP
  //   1540	C     DELTA0            PAIRING GAP 12 FOR GROUND STATE                 
  //   1541	C                       14 FOR SADDLE POINT                             
  //   1542	C     EITERA            HELP VARIABLE FOR TEMPERATURE ITERATION         
  //   1543	C-----------------------------------------------------------------------
  //   1544	C                                                                       
  //   1545	C                                                                       
  G4double afp = 0.0;
  G4double delta0 = 0.0;
  G4double deltau = 0.0;
  G4double deltpp = 0.0;
  G4double e = 0.0;
  G4double ecor = 0.0;
  G4double ecor1 = 0.0;
  G4double ecr = 0.0;
  G4double er = 0.0;
  G4double fe = 0.0;
  G4double fp = 0.0;
  G4double he = 0.0;
  G4double iz = 0.0;
  G4double pa = 0.0;
  G4double para = 0.0;
  G4double parz = 0.0;
  G4double ponfe = 0.0;
  G4double ponniv = 0.0;
  G4double qr = 0.0;
  G4double sig = 0.0;
  G4double y01 = 0.0;
  G4double y11 = 0.0;
  G4double y2 = 0.0;
  G4double y21 = 0.0;
  G4double y1 = 0.0;
  G4double y0 = 0.0;

  G4double pi6 = std::pow(3.1415926535,2) / 6.0;
  ecr=10.0;
  er=28.0;
  afp=idnint(a);
  iz=idnint(z);

  // level density parameter                                               
  if((ald->optafan == 1)) {
    pa = (ald->av)*a + (ald->as)*std::pow(a,(2.e0/3.e0)) + (ald->ak)*std::pow(a,(1.e0/3.e0));
  }
  else {
    pa = (ald->av)*a + (ald->as)*bs*std::pow(a,(2.e0/3.e0)) + (ald->ak)*bk*std::pow(a,(1.e0/3.e0));
  }

  fp = 0.01377937231e0 * std::pow(a,(5.e0/3.e0)) * (1.e0 + defbet/3.e0);

  // pairing corrections                                                   
  if (bs > 1.0) {
    delta0 = 14;
  }
  else {
    delta0 = 12;
  }

  if (esous > 1.0e30) {
    (*dens) = 0.0;
    (*temp) = 0.0;
    goto densniv100;                                                       
  }

  e = ee - esous;

  if (e < 0.0) {
    (*dens) = 0.0;
    (*temp) = 0.0;
    goto densniv100;
  }

  // shell corrections                                                     
  if (optshp > 0) {
    deltau = bshell;
    if (optshp == 2) {
      deltau = 0.0;
    }
    if (optshp >= 2) {
      // pairing energy shift with condensation energy a.r.j. 10.03.97        
      //      deltpp = -0.25e0* (delta0/std::pow(std::sqrt(a),2)) * pa /pi6 + 2.e0*delta0/std::sqrt(a);
      deltpp = -0.25e0* std::pow((delta0/std::sqrt(a)),2) * pa /pi6 + 2.e0*delta0/std::sqrt(a);
      // assert(isnan(deltpp) == false);
      
      parite(a,&para);
      if (para < 0.0) {
	e = e - delta0/std::sqrt(a);
	// assert(isnan(e) == false);
      } else {                                                         
	parite(z,&parz);
	if (parz > 0.e0) {
	  e = e - 2.0*delta0/std::sqrt(a);
	  // assert(isnan(e) == false);
	} else {
	  e = e;
	  // assert(isnan(e) == false);
	}
      }
    } else {                                                          
      deltpp = 0.0;
    }
  } else {
    deltau = 0.0;
    deltpp = 0.0;
  }
  if (e < 0.0) {
    e = 0.0;
    (*temp) = 0.0;
  }

  // washing out is made stronger ! g.k. 3.7.96                           
  ponfe = -2.5*pa*e*std::pow(a,(-4.0/3.0));

  if (ponfe < -700.0)  {
    ponfe = -700.0;
  }
  fe = 1.0 - std::exp(ponfe);
  if (e < ecr) {
    // priv. comm. k.-h. schmidt                                         
    he = 1.0 - std::pow((1.0 - e/ecr),2);
  }
  else {
    he = 1.0;
  }

  // Excitation energy corrected for pairing and shell effects             
  // washing out with excitation energy is included.                        
  ecor = e + deltau*fe + deltpp*he;

  if (ecor <= 0.1) {
    ecor = 0.1;
  }

  // statt 170.d0 a.r.j. 8.11.97                                           

  // iterative procedure according to grossjean and feldmeier              
  // to avoid the singularity e = 0                                        
  if (ee < 5.0) {
    y1 = std::sqrt(pa*ecor);
    // assert(isnan(y1) == false);
    for(int j = 0; j < 5; j++) {
      y2 = pa*ecor*(1.e0-std::exp(-y1));
      // assert(isnan(y2) == false);
      y1 = std::sqrt(y2);
      // assert(isnan(y1) == false);
    }
    
    y0 = pa/y1;
    // assert(isnan(y0) == false);
    assert(y0 != 0.0);
    (*temp)=1.0/y0;
    (*dens) = std::exp(y0*ecor)/ (std::pow((std::pow(ecor,3)*y0),0.5)*std::pow((1.0-0.5*y0*ecor*std::exp(-y1)),0.5))* std::exp(y1)*(1.0-std::exp(-y1))*0.1477045;
    if (ecor < 1.0) {
      ecor1=1.0;
      y11 = std::sqrt(pa*ecor1);
      // assert(isnan(y11) == false);
      for(int j = 0; j < 7; j++) {
	y21 = pa*ecor1*(1.0-std::exp(-y11));
	// assert(isnan(21) == false);
	y11 = std::sqrt(y21);
	// assert(isnan(y11) == false);
      }

      y01 = pa/y11;
      // assert(isnan(y01) == false);
      (*dens) = (*dens)*std::pow((y01/y0),1.5);
      (*temp) = (*temp)*std::pow((y01/y0),1.5);
    }
  }
  else {
    ponniv = 2.0*std::sqrt(pa*ecor);
    // assert(isnan(ponniv) == false);
    if (ponniv > 700.0) {
      ponniv = 700.0;
    }

    // fermi gas state density                                               
    (*dens) = std::pow(pa,(-0.25e0))*std::pow(ecor,(-1.25e0))*std::exp(ponniv) * 0.1477045e0;
    // assert(isnan(std::sqrt(ecor/pa)) == false);
    (*temp) = std::sqrt(ecor/pa);
  }
 densniv100:

  // spin cutoff parameter                                                 
  sig = fp * (*temp);

  // collective enhancement                                                
  if (optcol == 1) {
    qrot(z,a,defbet,sig,ecor,&qr);
  }
  else {
    qr   = 1.0;
  }

  (*dens) = (*dens) * qr;
}


G4double G4Abla::bfms67(G4double zms, G4double ams)
{
  // This subroutine calculates the fission barriers                                                                  
  // of the liquid-drop model of Myers and Swiatecki (1967).                                                                 
  // Analytic parameterization of Dahlinger 1982 
  // replaces tables. Barrier heights from Myers and Swiatecki !!!                                                                 

  G4double nms = 0.0, ims = 0.0, ksims = 0.0, xms = 0.0, ums = 0.0;

  nms = ams - zms;
  ims = (nms-zms)/ams;
  ksims= 50.15e0 * (1.- 1.78 * std::pow(ims,2));
  xms = std::pow(zms,2) / (ams * ksims);
  ums = 0.368e0-5.057e0*xms+8.93e0*std::pow(xms,2)-8.71*std::pow(xms,3);
  return(0.7322e0*std::pow(zms,2)/std::pow(ams,(0.333333e0))*std::pow(10.e0,ums));
}

void G4Abla::lpoly(G4double x, G4int n, G4double pl[])
{
  // THIS SUBROUTINE CALCULATES THE ORDINARY LEGENDRE POLYNOMIALS OF   
  // ORDER 0 TO N-1 OF ARGUMENT X AND STORES THEM IN THE VECTOR PL.    
  // THEY ARE CALCULATED BY RECURSION RELATION FROM THE FIRST TWO      
  // POLYNOMIALS.                                                      
  // WRITTEN BY A.J.SIERK  LANL  T-9  FEBRUARY, 1984                   

  // NOTE: PL AND X MUST BE DOUBLE PRECISION ON 32-BIT COMPUTERS!      

  pl[0] = 1.0;
  pl[1] = x;

  for(int i = 2; i < n; i++) {
    pl[i] = ((2*double(i+1) - 3.0)*x*pl[i-1] - (double(i+1) - 2.0)*pl[i-2])/(double(i+1)-1.0);
  }
}

G4double G4Abla::eflmac(G4int ia, G4int iz, G4int flag, G4int optshp)
{
  // CHANGED TO CALCULATE TOTAL BINDING ENERGY INSTEAD OF MASS EXCESS.     
  // SWITCH FOR PAIRING INCLUDED AS WELL.                                  
  // BINDING = EFLMAC(IA,IZ,0,OPTSHP)                                      
  // FORTRAN TRANSCRIPT OF /U/GREWE/LANG/EEX/FRLDM.C                       
  // A.J. 15.07.96                                                         

  // this function will calculate the liquid-drop nuclear mass for spheri
  // configuration according to the preprint NUCLEAR GROUND-STATE        
  // MASSES and DEFORMATIONS by P. M"oller et al. from August 16, 1993 p.
  // All constants are taken from this publication for consistency.      

  // Parameters:                                                         
  // a:    nuclear mass number                                         
  // z:    nuclear charge                                              
  // flag:     0       - return mass excess                            
  //       otherwise   - return pairing (= -1/2 dpn + 1/2 (Dp + Dn))   

  G4double eflmacResult = 0.0;

  G4int in = 0;
  G4double z = 0.0, n = 0.0, a = 0.0, av = 0.0, as = 0.0;
  G4double a0 = 0.0, c1 = 0.0, c4 = 0.0, b1 = 0.0, b3 = 0.0;
  G4double f = 0.0, ca = 0.0, w = 0.0, dp = 0.0, dn = 0.0, dpn = 0.0, efl = 0.0;
  G4double rmac = 0.0, bs = 0.0, h = 0.0, r0 = 0.0, kf = 0.0, ks = 0.0;
  G4double kv = 0.0, rp = 0.0, ay = 0.0, aden = 0.0, x0 = 0.0, y0 = 0.0;
  G4double mh = 0.0, mn = 0.0, esq = 0.0, ael = 0.0, i = 0.0;
  G4double pi = 3.141592653589793238e0;

  // fundamental constants
  // hydrogen-atom mass excess
  mh  = 7.289034;

  // neutron mass excess
  mn  = 8.071431;

  // electronic charge squared
  esq = 1.4399764;

  // constants from considerations other than nucl. masses
  // electronic binding
  ael = 1.433e-5;

  // proton rms radius
  rp  = 0.8;

  // nuclear radius constant
  r0  = 1.16;

  // range of yukawa-plus-expon. potential
  ay  = 0.68;

  // range of yukawa function used to generate                          
  // nuclear charge distribution
  aden= 0.70;

  // constants from considering odd-even mass differences
  // average pairing gap
  rmac= 4.80;

  // neutron-proton interaction
  h   = 6.6;

  // wigner constant
  w   = 30.0;

  // adjusted parameters
  // volume energy
  av  = 16.00126;

  // volume asymmetry
  kv  =  1.92240;

  // surface energy
  as  = 21.18466;

  // surface asymmetry
  ks  =  2.345;
  // a^0 constant
  a0  =  2.615;

  // charge asymmetry
  ca  =  0.10289;

  // we will account for deformation by using the microscopic           
  // corrections tabulated from p. 68ff */                               
  bs = 1.0;

  z   = double(iz);
  a   = double(ia);
  in  = ia - iz;                                                       
  n   = double(in);
  dn  = rmac*bs/std::pow(n,(1.0/3.0));
  dp  = rmac*bs/std::pow(z,(1.0/3.0));
  dpn = h/bs/std::pow(a,(2.0/3.0));
  // assert(isnan(dpn) == false);
  
  c1  = 3.0/5.0*esq/r0;
  // assert(isnan(c1) == false);
  // assert(isinf(c1) == false);
  
  c4  = 5.0/4.0*std::pow((3.0/(2.0*pi)),(2.0/3.0)) * c1;
  // assert(isnan(c4) == false);
  // assert(isinf(c4) == false);
 
  // assert(isnan(pi) == false);
  // assert(isnan(z) == false);
  // assert(isnan(a) == false);
  // assert(isnan(r0) == false);
  kf  = std::pow((9.0*pi*z/(4.0*a)),(1.0/3.0))/r0;
  // assert(isnan(kf) == false);
  // assert(isinf(kf) == false);
  
  f = -1.0/8.0*rp*rp*esq/std::pow(r0,3) * (145.0/48.0 - 327.0/2880.0*std::pow(kf,2) * std::pow(rp,2) + 1527.0/1209600.0*std::pow(kf,4) * std::pow(rp,4));
  i   = (n-z)/a;

  x0  = r0 * std::pow(a,(1.0/3.0)) / ay;
  y0  = r0 * std::pow(a,(1.0/3.0)) / aden;

  b1  = 1.0 - 3.0/(std::pow(x0,2)) + (1.0 + x0) * (2.0 + 3.0/x0 + 3.0/std::pow(x0,2)) * std::exp(-2.0*x0);

  b3  = 1.0 - 5.0/std::pow(y0,2) * (1.0 - 15.0/(8.0*y0) + 21.0/(8.0 * std::pow(y0,3))
			       - 3.0/4.0 * (1.0 + 9.0/(2.0*y0) + 7.0/std::pow(y0,2)
					    + 7.0/(2.0 * std::pow(y0,3))) * std::exp(-2.0*y0));

  // now calulation of total binding energy a.j. 16.7.96                   

  efl = -1.0 * av*(1.0 - kv*i*i)*a + as*(1.0 - ks*i*i)*b1 * std::pow(a,(2.0/3.0)) + a0
    + c1*z*z*b3/std::pow(a,(1.0/3.0)) - c4*std::pow(z,(4.0/3.0))/std::pow(a,(1.e0/3.e0))
    + f*std::pow(z,2)/a -ca*(n-z) - ael * std::pow(z,(2.39e0));

  if ((in == iz) && (mod(in,2) == 1) && (mod(iz,2) == 1)) {
    // n and z odd and equal
    efl = efl + w*(utilabs(i)+1.e0/a);
  }
  else {
    efl= efl + w* utilabs(i);
  }

  // pairing is made optional                                              
  if (optshp >= 2) {
    // average pairing
    if ((mod(in,2) == 1) && (mod(iz,2) == 1)) {
      efl = efl - dpn;
    }
    if (mod(in,2) == 1) {
      efl = efl + dn;
    }
    if (mod(iz,2) == 1) {
      efl    = efl + dp;
    }
    // end if for pairing term                                               
  }

  if (flag != 0) {
    eflmacResult =  (0.5*(dn + dp) - 0.5*dpn);
  }
  else {
    eflmacResult = efl;
  }

  return eflmacResult;
}

void G4Abla::appariem(G4double a, G4double z, G4double *del)
{
  // CALCUL DE LA CORRECTION, DUE A L'APPARIEMENT, DE L'ENERGIE DE     
  // LIAISON D'UN NOYAU                                                
  // PROCEDURE FOR CALCULATING THE PAIRING CORRECTION TO THE BINDING   
  // ENERGY OF A SPECIFIC NUCLEUS                                      

  double para = 0.0, parz = 0.0;
  // A                 MASS NUMBER                                     
  // Z                 NUCLEAR CHARGE                                  
  // PARA              HELP VARIABLE FOR PARITY OF A                   
  // PARZ              HELP VARIABLE FOR PARITY OF Z                   
  // DEL               PAIRING CORRECTION                              

  parite(a, &para);

  if (para < 0.0) {
    (*del) = 0.0;
  }
  else {
    parite(z, &parz);
    if (parz > 0.0) {
      // assert(isnan(std::sqrt(a)) == false);
      (*del) = -12.0/std::sqrt(a);
    }
    else {
      // assert(isnan(std::sqrt(a)) == false);
      (*del) = 12.0/std::sqrt(a);
    }
  }
}

void G4Abla::parite(G4double n, G4double *par)
{
  // CALCUL DE LA PARITE DU NOMBRE N                                   
  //
  // PROCEDURE FOR CALCULATING THE PARITY OF THE NUMBER N.             
  // RETURNS -1 IF N IS ODD AND +1 IF N IS EVEN                        

  G4double n1 = 0.0, n2 = 0.0, n3 = 0.0;

  // N                 NUMBER TO BE TESTED                             
  // N1,N2             HELP VARIABLES                                  
  // PAR               HELP VARIABLE FOR PARITY OF N                   

  n3 = double(idnint(n));
  n1 = n3/2.0;
  n2 = n1 - dint(n1);

  if (n2 > 0.0) {
    (*par) = -1.0;
  }
  else {
    (*par) = 1.0;
  }
}

G4double G4Abla::tau(G4double bet, G4double homega, G4double ef, G4double t)
{
  // INPUT : BET, HOMEGA, EF, T                                          
  // OUTPUT: TAU - RISE TIME IN WHICH THE FISSION WIDTH HAS REACHED      
  //               90 PERCENT OF ITS FINAL VALUE                               
  // 
  // BETA   - NUCLEAR VISCOSITY                                          
  // HOMEGA - CURVATURE OF POTENTIAL                                     
  // EF     - FISSION BARRIER                                            
  // T      - NUCLEAR TEMPERATURE                                        

  G4double tauResult = 0.0;

  G4double tlim = 8.e0 * ef;
  if (t > tlim) {
    t = tlim;
  }

  // modified bj and khs 6.1.2000:
  if (bet/(std::sqrt(2.0)*10.0*(homega/6.582122)) <= 1.0) {
    tauResult = std::log(10.0*ef/t)/(bet*1.0e21);
    //    assert(isnan(tauResult) == false);
  }
  else {
    tauResult = std::log(10.0*ef/t)/ (2.0*std::pow((10.0*homega/6.582122),2))*(bet*1.0e-21);
    // assert(isnan(tauResult) == false);
  } //end if                                                            

  return tauResult;
}

G4double G4Abla::cram(G4double bet, G4double homega)
{
  // INPUT : BET, HOMEGA  NUCLEAR VISCOSITY + CURVATURE OF POTENTIAL      
  // OUTPUT: KRAMERS FAKTOR  - REDUCTION OF THE FISSION PROBABILITY       
  //                           INDEPENDENT OF EXCITATION ENERGY                             

  G4double rel = bet/(20.0*homega/6.582122);
  G4double cramResult = std::sqrt(1.0 + std::pow(rel,2)) - rel;
  // limitation introduced   6.1.2000  by  khs

  if (cramResult > 1.0) {
    cramResult = 1.0;
  }

  // assert(isnan(cramResult) == false);
  return cramResult;
}

G4double G4Abla::bipol(int iflag, G4double y)
{
  // CALCULATION OF THE SURFACE BS OR CURVATURE BK OF A NUCLEUS        
  // RELATIVE TO THE SPHERICAL CONFIGURATION                           
  // BASED ON  MYERS, DROPLET MODEL FOR ARBITRARY SHAPES               

  // INPUT: IFLAG - 0/1 BK/BS CALCULATION                              
  //         Y    - (1 - X) COMPLEMENT OF THE FISSILITY                

  // LINEAR INTERPOLATION OF BS BK TABLE                               

  int i = 0;

  G4double bipolResult = 0.0;

  const int bsbkSize = 54;

  G4double bk[bsbkSize] = {0.0, 1.00000,1.00087,1.00352,1.00799,1.01433,1.02265,1.03306,
			   1.04576,1.06099,1.07910,1.10056,1.12603,1.15651,1.19348,
			   1.23915,1.29590,1.35951,1.41013,1.44103,1.46026,1.47339,
			   1.48308,1.49068,1.49692,1.50226,1.50694,1.51114,1.51502,
			   1.51864,1.52208,1.52539,1.52861,1.53177,1.53490,1.53803,
			   1.54117,1.54473,1.54762,1.55096,1.55440,1.55798,1.56173,
			   1.56567,1.56980,1.57413,1.57860,1.58301,1.58688,1.58688,
			   1.58688,1.58740,1.58740, 0.0}; //Zeroes at bk[0], and at the end added by PK

  G4double bs[bsbkSize] = {0.0, 1.00000,1.00086,1.00338,1.00750,1.01319,
			   1.02044,1.02927,1.03974,
			   1.05195,1.06604,1.08224,1.10085,1.12229,1.14717,1.17623,1.20963,
			   1.24296,1.26532,1.27619,1.28126,1.28362,1.28458,1.28477,1.28450,
			   1.28394,1.28320,1.28235,1.28141,1.28042,1.27941,1.27837,1.27732,
			   1.27627,1.27522,1.27418,1.27314,1.27210,1.27108,1.27006,1.26906,
			   1.26806,1.26707,1.26610,1.26514,1.26418,1.26325,1.26233,1.26147,
			   1.26147,1.26147,1.25992,1.25992, 0.0};

  i = idint(y/(2.0e-02)) + 1;
  assert(i >= 1);
    
  if(i >= bsbkSize) {
    if(verboseLevel > 2) {
      G4cout <<"G4Abla error: index i = " << i << " is greater than array size permits." << G4endl;
    }
    bipolResult = 0.0;
  }
  else {
    if (iflag == 1) {
      bipolResult = bs[i] + (bs[i+1] - bs[i])/2.0e-02 * (y - 2.0e-02*(i - 1));
    }
    else {
      bipolResult = bk[i] + (bk[i+1] - bk[i])/2.0e-02 * (y - 2.0e-02*(i - 1));
    }
  }
  
  // assert(isnan(bipolResult) == false);
  return bipolResult;
}

void G4Abla::barfit(G4int iz, G4int ia, G4int il, G4double *sbfis, G4double *segs, G4double *selmax)
{
  //   2223	C     VERSION FOR 32BIT COMPUTER                                        
  //   2224	C     THIS SUBROUTINE RETURNS THE BARRIER HEIGHT BFIS, THE              
  //   2225	C     GROUND-STATE ENERGY SEGS, IN MEV, AND THE ANGULAR MOMENTUM        
  //   2226	C     AT WHICH THE FISSION BARRIER DISAPPEARS, LMAX, IN UNITS OF        
  //   2227	C     H-BAR, WHEN CALLED WITH INTEGER AGUMENTS IZ, THE ATOMIC           
  //   2228	C     NUMBER, IA, THE ATOMIC MASS NUMBER, AND IL, THE ANGULAR           
  //   2229	C     MOMENTUM IN UNITS OF H-BAR. (PLANCK'S CONSTANT DIVIDED BY         
  //   2230	C     2*PI).                                                            
  //   2231	C                                                                       
  //   2232	C        THE FISSION BARRIER FO IL = 0 IS CALCULATED FROM A 7TH         
  //   2233	C     ORDER FIT IN TWO VARIABLES TO 638 CALCULATED FISSION              
  //   2234	C     BARRIERS FOR Z VALUES FROM 20 TO 110. THESE 638 BARRIERS ARE      
  //   2235	C     FIT WITH AN RMS DEVIATION OF 0.10 MEV BY THIS 49-PARAMETER        
  //   2236	C     FUNCTION.                                                         
  //   2237	C     IF BARFIT IS CALLED WITH (IZ,IA) VALUES OUTSIDE THE RANGE OF      
  //   2238	C     THE BARRIER HEIGHT IS SET TO 0.0, AND A MESSAGE IS PRINTED        
  //   2239	C     ON THE DEFAULT OUTPUT FILE.                                       
  //   2240	C                                                                       
  //   2241	C        FOR IL VALUES NOT EQUAL TO ZERO, THE VALUES OF L AT WHICH      
  //   2242	C     THE BARRIER IS 80% AND 20% OF THE L=0 VALUE ARE RESPECTIVELY      
  //   2243	C     FIT TO 20-PARAMETER FUNCTIONS OF Z AND A, OVER A MORE             
  //   2244	C     RESTRICTED RANGE OF A VALUES, THAN IS THE CASE FOR L = 0.         
  //   2245	C     THE VALUE OF L WHERE THE BARRIER DISAPPEARS, LMAX IS FIT TO       
  //   2246	C     A 24-PARAMETER FUNCTION OF Z AND A, WITH THE SAME RANGE OF        
  //   2247	C     Z AND A VALUES AS L-80 AND L-20.                                  
  //   2248	C        ONCE AGAIN, IF AN (IZ,IA) PAIR IS OUTSIDE OF THE RANGE OF      
  //   2249	C     VALIDITY OF THE FIT, THE BARRIER VALUE IS SET TO 0.0 AND A        
  //   2250	C     MESSAGE IS PRINTED. THESE THREE VALUES (BFIS(L=0),L-80, AND       
  //   2251	C     L-20) AND THE CONSTRINTS OF BFIS = 0 AND D(BFIS)/DL = 0 AT        
  //   2252	C     L = LMAX AND L=0 LEAD TO A FIFTH-ORDER FIT TO BFIS(L) FOR         
  //   2253	C     L>L-20. THE FIRST THREE CONSTRAINTS LEAD TO A THIRD-ORDER FIT     
  //   2254	C     FOR THE REGION L < L-20.                                          
  //   2255	C                                                                       
  //   2256	C        THE GROUND STATE ENERGIES ARE CALCULATED FROM A                
  //   2257	C     120-PARAMETER FIT IN Z, A, AND L TO 214 GROUND-STATE ENERGIES     
  //   2258	C     FOR 36 DIFFERENT Z AND A VALUES.                                  
  //   2259	C     (THE RANGE OF Z AND A IS THE SAME AS FOR L-80, L-20, AND          
  //   2260	C     L-MAX)                                                            
  //   2261	C                                                                       
  //   2262	C        THE CALCULATED BARRIERS FROM WHICH THE FITS WERE MADE WERE     
  //   2263	C     CALCULATED IN 1983-1984 BY A. J. SIERK OF LOS ALAMOS              
  //   2264	C     NATIONAL LABORATORY GROUP T-9, USING YUKAWA-PLUS-EXPONENTIAL      
  //   2265	C     G4DOUBLE FOLDED NUCLEAR ENERGY, EXACT COULOMB DIFFUSENESS           
  //   2266	C     CORRECTIONS, AND DIFFUSE-MATTER MOMENTS OF INERTIA.               
  //   2267	C     THE PARAMETERS OF THE MODEL R-0 = 1.16 FM, AS 21.13 MEV,          
  //   2268	C     KAPPA-S = 2.3, A = 0.68 FM.                                       
  //   2269	C     THE DIFFUSENESS OF THE MATTER AND CHARGE DISTRIBUTIONS USED       
  //   2270	C     CORRESPONDS TO A SURFACE DIFFUSENESS PARAMETER (DEFINED BY        
  //   2271	C     MYERS) OF 0.99 FM. THE CALCULATED BARRIERS FOR L = 0 ARE          
  //   2272	C     ACCURATE TO A LITTLE LESS THAN 0.1 MEV; THE OUTPUT FROM           
  //   2273	C     THIS SUBROUTINE IS A LITTLE LESS ACCURATE. WORST ERRORS MAY BE    
  //   2274	C     AS LARGE AS 0.5 MEV; CHARACTERISTIC UNCERTAINY IS IN THE RANGE    
  //   2275	C     OF 0.1-0.2 MEV. THE RMS DEVIATION OF THE GROUND-STATE FIT         
  //   2276	C     FROM THE 214 INPUT VALUES IS 0.20 MEV. THE MAXIMUM ERROR          
  //   2277	C     OCCURS FOR LIGHT NUCLEI IN THE REGION WHERE THE GROUND STATE      
  //   2278	C     IS PROLATE, AND MAY BE GREATER THAN 1.0 MEV FOR VERY NEUTRON      
  //   2279	C     DEFICIENT NUCLEI, WITH L NEAR LMAX. FOR MOST NUCLEI LIKELY TO     
  //   2280	C     BE ENCOUNTERED IN REAL EXPERIMENTS, THE MAXIMUM ERROR IS          
  //   2281	C     CLOSER TO 0.5 MEV, AGAIN FOR LIGHT NUCLEI AND L NEAR LMAX.        
  //   2282	C                                                                       
  //   2283	C     WRITTEN BY A. J. SIERK, LANL T-9                                  
  //   2284	C     VERSION 1.0 FEBRUARY, 1984                                        
  //   2285	C                                                                       
  //   2286	C     THE FOLLOWING IS NECESSARY FOR 32-BIT MACHINES LIKE DEC VAX,      
  //   2287	C     IBM, ETC                                                          

  G4double pa[7],pz[7],pl[10];
  for(G4int init_i = 0; init_i < 7; init_i++) {
    pa[init_i] = 0.0; 
    pz[init_i] = 0.0; 
  }
  for(G4int init_i = 0; init_i < 10; init_i++) {
    pl[init_i] = 0.0;
  }

  G4double a = 0.0, z = 0.0, amin = 0.0, amax = 0.0, amin2 = 0.0;
  G4double amax2 = 0.0, aa = 0.0, zz = 0.0, bfis = 0.0;
  G4double bfis0 = 0.0, ell = 0.0, el = 0.0, egs = 0.0, el80 = 0.0, el20 = 0.0;
  G4double elmax = 0.0, sel80 = 0.0, sel20 = 0.0, x = 0.0, y = 0.0, q = 0.0, qa = 0.0, qb = 0.0;
  G4double aj = 0.0, ak = 0.0, a1 = 0.0, a2 = 0.0;

  G4int i = 0, j = 0, k = 0, m = 0;
  G4int l = 0;

  G4double emncof[4][5] = {{-9.01100e+2,-1.40818e+3, 2.77000e+3,-7.06695e+2, 8.89867e+2}, 
			   {1.35355e+4,-2.03847e+4, 1.09384e+4,-4.86297e+3,-6.18603e+2},
			   {-3.26367e+3, 1.62447e+3, 1.36856e+3, 1.31731e+3, 1.53372e+2},
			   {7.48863e+3,-1.21581e+4, 5.50281e+3,-1.33630e+3, 5.05367e-2}};

  G4double elmcof[4][5] = {{1.84542e+3,-5.64002e+3, 5.66730e+3,-3.15150e+3, 9.54160e+2},
			   {-2.24577e+3, 8.56133e+3,-9.67348e+3, 5.81744e+3,-1.86997e+3},
			   {2.79772e+3,-8.73073e+3, 9.19706e+3,-4.91900e+3, 1.37283e+3},
			   {-3.01866e+1, 1.41161e+3,-2.85919e+3, 2.13016e+3,-6.49072e+2}};

  G4double emxcof[4][6] = {{9.43596e4,-2.241997e5,2.223237e5,-1.324408e5,4.68922e4,-8.83568e3},
			   {-1.655827e5,4.062365e5,-4.236128e5,2.66837e5,-9.93242e4,1.90644e4},
			   {1.705447e5,-4.032e5,3.970312e5,-2.313704e5,7.81147e4,-1.322775e4},
			   {-9.274555e4,2.278093e5,-2.422225e5,1.55431e5,-5.78742e4,9.97505e3}};

  G4double elzcof[7][7] = {{5.11819909e+5,-1.30303186e+6, 1.90119870e+6,-1.20628242e+6, 5.68208488e+5, 5.48346483e+4,-2.45883052e+4},
			   {-1.13269453e+6, 2.97764590e+6,-4.54326326e+6, 3.00464870e+6, -1.44989274e+6,-1.02026610e+5, 6.27959815e+4},
			   {1.37543304e+6,-3.65808988e+6, 5.47798999e+6,-3.78109283e+6, 1.84131765e+6, 1.53669695e+4,-6.96817834e+4},
			   {-8.56559835e+5, 2.48872266e+6,-4.07349128e+6, 3.12835899e+6, -1.62394090e+6, 1.19797378e+5, 4.25737058e+4},
			   {3.28723311e+5,-1.09892175e+6, 2.03997269e+6,-1.77185718e+6, 9.96051545e+5,-1.53305699e+5,-1.12982954e+4},
			   {4.15850238e+4, 7.29653408e+4,-4.93776346e+5, 6.01254680e+5, -4.01308292e+5, 9.65968391e+4,-3.49596027e+3},
			   {-1.82751044e+5, 3.91386300e+5,-3.03639248e+5, 1.15782417e+5, -4.24399280e+3,-6.11477247e+3, 3.66982647e+2}};

  const G4int sizex = 5;
  const G4int sizey = 6;
  const G4int sizez = 4;

  G4double egscof[sizey][sizey][sizez];

  G4double egs1[sizey][sizex] = {{1.927813e5, 7.666859e5, 6.628436e5, 1.586504e5,-7.786476e3},
				 {-4.499687e5,-1.784644e6,-1.546968e6,-4.020658e5,-3.929522e3},
				 {4.667741e5, 1.849838e6, 1.641313e6, 5.229787e5, 5.928137e4},
				 {-3.017927e5,-1.206483e6,-1.124685e6,-4.478641e5,-8.682323e4},
				 {1.226517e5, 5.015667e5, 5.032605e5, 2.404477e5, 5.603301e4},
				 {-1.752824e4,-7.411621e4,-7.989019e4,-4.175486e4,-1.024194e4}};

  G4double egs2[sizey][sizex] = {{-6.459162e5,-2.903581e6,-3.048551e6,-1.004411e6,-6.558220e4},
				 {1.469853e6, 6.564615e6, 6.843078e6, 2.280839e6, 1.802023e5},
				 {-1.435116e6,-6.322470e6,-6.531834e6,-2.298744e6,-2.639612e5},
				 {8.665296e5, 3.769159e6, 3.899685e6, 1.520520e6, 2.498728e5},      
				 {-3.302885e5,-1.429313e6,-1.512075e6,-6.744828e5,-1.398771e5},
				 {4.958167e4, 2.178202e5, 2.400617e5, 1.167815e5, 2.663901e4}};

  G4double egs3[sizey][sizex] = {{3.117030e5, 1.195474e6, 9.036289e5, 6.876190e4,-6.814556e4},
				 {-7.394913e5,-2.826468e6,-2.152757e6,-2.459553e5, 1.101414e5},
				 {7.918994e5, 3.030439e6, 2.412611e6, 5.228065e5, 8.542465e3},
				 {-5.421004e5,-2.102672e6,-1.813959e6,-6.251700e5,-1.184348e5},
				 {2.370771e5, 9.459043e5, 9.026235e5, 4.116799e5, 1.001348e5},
				 {-4.227664e4,-1.738756e5,-1.795906e5,-9.292141e4,-2.397528e4}};

  G4double egs4[sizey][sizex] = {{-1.072763e5,-5.973532e5,-6.151814e5, 7.371898e4, 1.255490e5},
				 {2.298769e5, 1.265001e6, 1.252798e6,-2.306276e5,-2.845824e5},
				 {-2.093664e5,-1.100874e6,-1.009313e6, 2.705945e5, 2.506562e5},
				 {1.274613e5, 6.190307e5, 5.262822e5,-1.336039e5,-1.115865e5},
				 {-5.715764e4,-2.560989e5,-2.228781e5,-3.222789e3, 1.575670e4},
				 {1.189447e4, 5.161815e4, 4.870290e4, 1.266808e4, 2.069603e3}};

  for(i = 0; i < sizey; i++) {
    for(j = 0; j < sizex; j++) {
      //       egscof[i][j][0] = egs1[i][j];
      //       egscof[i][j][1] = egs2[i][j];
      //       egscof[i][j][2] = egs3[i][j];
      //       egscof[i][j][3] = egs4[i][j];
      egscof[i][j][0] = egs1[i][j];
      egscof[i][j][1] = egs2[i][j];
      egscof[i][j][2] = egs3[i][j];
      egscof[i][j][3] = egs4[i][j];
    }
  }

  // the program starts here                                           
  if (iz < 19  ||  iz > 111) {
    goto barfit900;
  }

  if(iz > 102   &&  il > 0) {
    goto barfit902;
  }

  z=double(iz);
  a=double(ia);
  el=double(il);
  amin= 1.2e0*z + 0.01e0*z*z;
  amax= 5.8e0*z - 0.024e0*z*z;

  if(a  <  amin  ||  a  >  amax) {
    goto barfit910;
  }

  // angul.mom.zero barrier                 
  aa=2.5e-3*a;
  zz=1.0e-2*z;
  ell=1.0e-2*el;
  bfis0 = 0.0;
  lpoly(zz,7,pz);
  lpoly(aa,7,pa);

  for(i = 0; i < 7; i++) { //do 10 i=1,7                                                       
    for(j = 0; j < 7; j++) { //do 10 j=1,7                                                       
      bfis0=bfis0+elzcof[j][i]*pz[i]*pa[j];
      //bfis0=bfis0+elzcof[i][j]*pz[j]*pa[i];
    }
  }

  bfis=bfis0;
  // assert(isnan(bfis) == false);
  
  (*sbfis)=bfis;
  egs=0.0;
  (*segs)=egs;

  // values of l at which the barrier        
  // is 20%(el20) and 80%(el80) of l=0 value    
  amin2 = 1.4e0*z + 0.009e0*z*z;
  amax2 = 20.e0 + 3.0e0*z;

  if((a < amin2-5.e0  ||  a > amax2+10.e0) &&  il > 0) {
    goto barfit920;
  }

  lpoly(zz,5,pz);
  lpoly(aa,4,pa);
  el80=0.0;
  el20=0.0;
  elmax=0.0;

  for(i = 0; i < 4; i++) {
    for(j = 0; j < 5; j++) {
//       el80 = el80 + elmcof[j][i]*pz[j]*pa[i];
//       el20 = el20 + emncof[j][i]*pz[j]*pa[i];
            el80 = el80 + elmcof[i][j]*pz[j]*pa[i];
            el20 = el20 + emncof[i][j]*pz[j]*pa[i];
    }
  }

  sel80 = el80;
  sel20 = el20;

  // value of l (elmax) where barrier disapp.
  lpoly(zz,6,pz);
  lpoly(ell,9,pl);

  for(i = 0; i < 4; i++) { //do 30 i= 1,4                                                      
    for(j = 0; j < 6; j++) { //do 30 j=1,6                                                       
      //elmax = elmax + emxcof[j][i]*pz[j]*pa[i];
      //      elmax = elmax + emxcof[j][i]*pz[i]*pa[j];
      elmax = elmax + emxcof[i][j]*pz[j]*pa[i];
    }
  }

  // assert(isnan(elmax) == false);
  (*selmax)=elmax;

  // value of barrier at ang.mom.  l          
  if(il < 1){
    return;                                                
  }

  x = sel20/(*selmax);
  // assert(isnan(x) == false);
  y = sel80/(*selmax);
  // assert(isnan(y) == false);
  
  if(el <= sel20) {
    // low l              
    q = 0.2/(std::pow(sel20,2)*std::pow(sel80,2)*(sel20-sel80));
    qa = q*(4.0*std::pow(sel80,3) - std::pow(sel20,3));
    qb = -q*(4.0*std::pow(sel80,2) - std::pow(sel20,2));
    bfis = bfis*(1.0 + qa*std::pow(el,2) + qb*std::pow(el,3));
  }
  else {
    // high l             
    aj = (-20.0*std::pow(x,5) + 25.e0*std::pow(x,4) - 4.0)*std::pow((y-1.0),2)*y*y;
    ak = (-20.0*std::pow(y,5) + 25.0*std::pow(y,4) - 1.0) * std::pow((x-1.0),2)*x*x;
    q = 0.2/(std::pow((y-x)*((1.0-x)*(1.0-y)*x*y),2));
    qa = q*(aj*y - ak*x);
    qb = -q*(aj*(2.0*y + 1.0) - ak*(2.0*x + 1.0));
    z = el/(*selmax);
    a1 = 4.0*std::pow(z,5) - 5.0*std::pow(z,4) + 1.0;
    a2 = qa*(2.e0*z + 1.e0);
    bfis=bfis*(a1 + (z - 1.e0)*(a2 + qb*z)*z*z*(z - 1.e0));
  }
  
  if(bfis <= 0.0) {
    bfis=0.0;
  }

  if(el > (*selmax)) {
    bfis=0.0;
  }
  (*sbfis)=bfis;

  // now calculate rotating ground state energy                        
  if(el > (*selmax)) {
    return;                                           
  }

  for(k = 0; k < 4; k++) {
    for(l = 0; l < 6; l++) {
      for(m = 0; m < 5; m++) {
	//egs = egs + egscof[l][m][k]*pz[l]*pa[k]*pl[2*m-1];
	egs = egs + egscof[l][m][k]*pz[l]*pa[k]*pl[2*m];
	// egs = egs + egscof[m][l][k]*pz[l]*pa[k]*pl[2*m-1];
      }
    }
  }

  (*segs)=egs;
  if((*segs) < 0.0) {
    (*segs)=0.0;
  }

  return;                                                            

 barfit900:  //continue                                                          
  (*sbfis)=0.0;
  // for z<19 sbfis set to 1.0e3                                            
  if (iz < 19)  {
    (*sbfis) = 1.0e3;
  }
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                            

 barfit902:
  (*sbfis)=0.0;
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                            

 barfit910:
  (*sbfis)=0.0;
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                           

 barfit920:
  (*sbfis)=0.0;
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                            
}

G4double G4Abla::expohaz(G4int k, G4double T)
{
  // TIRAGE ALEATOIRE DANS UNE EXPONENTIELLLE : Y=EXP(-X/T)

  // assert(isnan((-1*T*std::log(haz(k)))) == false);
  return (-1.0*T*std::log(haz(k)));
}

G4double G4Abla::fd(G4double E)
{
  // DISTRIBUTION DE MAXWELL

  return (E*std::exp(-E));
}

G4double G4Abla::f(G4double E)
{
  // FONCTION INTEGRALE DE FD(E)
  return (1.0 - (E + 1.0) * std::exp(-E));
}

G4double G4Abla::fmaxhaz(G4double T)
{
  // tirage aleatoire dans une maxwellienne
  // t : temperature
  //
  // declaration des variables
  //

  const int pSize = 101;
  static G4double p[pSize];

  // ial generateur pour le cascade (et les iy pour eviter les correlations)
  static G4int i = 0;
  static G4int itest = 0;
  // programme principal

  // calcul des p(i) par approximation de newton
  p[pSize-1] = 8.0;
  G4double x = 0.1;
  G4double x1 = 0.0;
  G4double y = 0.0;

  if (itest == 1) {
    goto fmaxhaz120;
  }

  for(i = 1; i <= 99; i++) {
  fmaxhaz20:
    x1 = x - (f(x) - double(i)/100.0)/fd(x);
    x = x1;
    if (std::fabs(f(x) - double(i)/100.0) < 1e-5) {
      goto fmaxhaz100;
    }
    goto fmaxhaz20;
  fmaxhaz100:
    p[i] = x;
  } //end do

  //  itest = 1;
  itest = 0;
  // tirage aleatoire et calcul du x correspondant 
  // par regression lineaire
 fmaxhaz120:
  standardRandom(&y, &(hazard->igraine[17]));
  i = nint(y*100);

  //   2590	c ici on evite froidement les depassements de tableaux....(a.b. 3/9/99)        
  if(i == 0) {
    goto fmaxhaz120;
  }

  if (i == 1) {
    x = p[i]*y*100;
  }
  else {
    x = (p[i] - p[i-1])*(y*100 - i) + p[i];
  }

  return(x*T);
}

G4double G4Abla::pace2(G4double a, G4double z)
{
  // PACE2
  // Cette fonction retourne le defaut de masse du noyau A,Z en MeV
  // Révisée pour a, z flottants 25/4/2002	                       =

  G4double pace2 = 0.0;

  G4int ii = idint(a+0.5);
  G4int jj = idint(z+0.5);

  if(ii <= 0 || jj < 0) {
    pace2=0.;
    return pace2;
  }

  if(jj > 300) {
    pace2=0.0;
  }
  else {
    pace2=pace->dm[ii][jj];
  }
  pace2=pace2/1000.;

  if(pace->dm[ii][jj] == 0.) {
    if(ii < 12) {
      pace2=-500.;
    }
    else {
      guet(&a, &z, &pace2);
      pace2=pace2-ii*931.5;
      pace2=pace2/1000.;
    }
  }

  return pace2;
}

void G4Abla::guet(G4double *x_par, G4double *z_par, G4double *find_par)
{
  // TABLE DE MASSES ET FORMULE DE MASSE TIRE DU PAPIER DE BRACK-GUET
  // Gives the theoritical value for mass excess...
  // Révisée pour x, z flottants 25/4/2002

  //real*8 x,z
  //	dimension q(0:50,0:70)
  G4double x = (*x_par);
  G4double z = (*z_par);
  G4double find = (*find_par);

  const G4int qrows = 50;
  const G4int qcols = 70;
  G4double q[qrows][qcols];
  for(G4int init_i = 0; init_i < qrows; init_i++) {
    for(G4int init_j = 0; init_j < qcols; init_j++) {
      q[init_i][init_j] = 0.0;
    }
  }

  G4int ix=G4int(std::floor(x+0.5));
  G4int iz=G4int(std::floor(z+0.5));
  G4double zz = iz;
  G4double xx = ix;
  find = 0.0;
  G4double avol = 15.776;
  G4double asur = -17.22;
  G4double ac = -10.24;
  G4double azer = 8.0;
  G4double xjj = -30.03;
  G4double qq = -35.4;
  G4double c1 = -0.737;
  G4double c2 = 1.28;

  if(ix <= 7) {
    q[0][1]=939.50;
    q[1][1]=938.21;
    q[1][2]=1876.1;
    q[1][3]=2809.39;
    q[2][4]=3728.34;
    q[2][3]=2809.4;
    q[2][5]=4668.8;
    q[2][6]=5606.5;
    q[3][5]=4669.1;
    q[3][6]=5602.9;
    q[3][7]=6535.27;
    q[4][6]=5607.3;
    q[4][7]=6536.1;
    q[5][7]=6548.3;
    find=q[iz][ix];
  }
  else {
    G4double xneu=xx-zz;
    G4double si=(xneu-zz)/xx;
    G4double x13=std::pow(xx,.333);
    G4double ee1=c1*zz*zz/x13;
    G4double ee2=c2*zz*zz/xx;
    G4double aux=1.+(9.*xjj/4./qq/x13);
    G4double ee3=xjj*xx*si*si/aux;
    G4double ee4=avol*xx+asur*(std::pow(xx,.666))+ac*x13+azer;
    G4double tota = ee1 + ee2 + ee3 + ee4;
    find = 939.55*xneu+938.77*zz - tota;
  }

  (*x_par) = x;
  (*z_par) = z;
  (*find_par) = find;
}


// Fission code

void G4Abla::even_odd(G4double r_origin,G4double r_even_odd,G4int &i_out)     
{
  // Procedure to calculate I_OUT from R_IN in a way that
  // on the average a flat distribution in R_IN results in a
  // fluctuating distribution in I_OUT with an even-odd effect as
  // given by R_EVEN_ODD

  //     /* ------------------------------------------------------------ */
  //     /* EXAMPLES :                                                   */
  //     /* ------------------------------------------------------------ */
  //     /*    If R_EVEN_ODD = 0 :                                       */
  //     /*           CEIL(R_IN)  ----                                   */
  //     /*                                                              */
  //     /*              R_IN ->                                         */
  //     /*            (somewhere in between CEIL(R_IN) and FLOOR(R_IN)) */                                            */
  //     /*                                                              */
  //     /*           FLOOR(R_IN) ----       --> I_OUT                   */
  //     /* ------------------------------------------------------------ */
  //     /*    If R_EVEN_ODD > 0 :                                       */
  //     /*      The interval for the above treatment is                 */
  //     /*         larger for FLOOR(R_IN) = even and                    */
  //     /*         smaller for FLOOR(R_IN) = odd                        */
  //     /*    For R_EVEN_ODD < 0 : just opposite treatment              */
  //     /* ------------------------------------------------------------ */

  //     /* ------------------------------------------------------------ */
  //     /* On input:   R_ORIGIN    nuclear charge (real number)         */
  //     /*             R_EVEN_ODD  requested even-odd effect            */
  //     /* Intermediate quantity: R_IN = R_ORIGIN + 0.5                 */
  //     /* On output:  I_OUT       nuclear charge (integer)             */
  //     /* ------------------------------------------------------------ */

  //      G4double R_ORIGIN,R_IN,R_EVEN_ODD,R_REST,R_HELP;
  G4double r_in = 0.0, r_rest = 0.0, r_help = 0.0;
  G4double r_floor = 0.0;
  G4double r_middle = 0.0;
  //      G4int I_OUT,N_FLOOR;
  G4int n_floor = 0;

  r_in = r_origin + 0.5;
  r_floor = (float)((int)(r_in));
  if (r_even_odd < 0.001) {
    i_out = (int)(r_floor);
  } 
  else {
    r_rest = r_in - r_floor;
    r_middle = r_floor + 0.5;
    n_floor = (int)(r_floor);
    if (n_floor%2 == 0) {
      // even before modif.
      r_help = r_middle + (r_rest - 0.5) * (1.0 - r_even_odd);
    } 
    else {
      // odd before modification
      r_help = r_middle + (r_rest - 0.5) * (1.0 + r_even_odd);
    }
    i_out = (int)(r_help);
  }
}

G4double G4Abla::umass(G4double z,G4double n,G4double beta)
{
  // liquid-drop mass, Myers & Swiatecki, Lysekil, 1967
  // pure liquid drop, without pairing and shell effects

  // On input:    Z     nuclear charge of nucleus
  //              N     number of neutrons in nucleus
  //              beta  deformation of nucleus
  // On output:   binding energy of nucleus

  G4double a = 0.0, umass = 0.0;
  G4double alpha = 0.0;
  G4double xcom = 0.0, xvs = 0.0, xe = 0.0;
  const G4double pi = 3.1416;
     
  a = n + z;
  alpha = ( std::sqrt(5.0/(4.0*pi)) ) * beta;
  // assert(isnan(alpha) == false);
  
  xcom = 1.0 - 1.7826 * ((a - 2.0*z)/a)*((a - 2.0*z)/a);
  // assert(isnan(xcom) == false);
  // factor for asymmetry dependence of surface and volume term
  xvs = - xcom * ( 15.4941 * a - 
		   17.9439 * std::pow(a,0.66667) * (1.0+0.4*alpha*alpha) );
  // sum of volume and surface energy
  xe = z*z * (0.7053/(std::pow(a,0.33333)) * (1.0-0.2*alpha*alpha) - 1.1529/a);
  // assert(isnan(xe) == false);
  umass = xvs + xe;
  
  return umass;
}

G4double G4Abla::ecoul(G4double z1,G4double n1,G4double beta1,G4double z2,G4double n2,G4double beta2,G4double d)
{
  // Coulomb potential between two nuclei
  // surfaces are in a distance of d
  // in a tip to tip configuration

  // approximate formulation
  // On input: Z1      nuclear charge of first nucleus
  //           N1      number of neutrons in first nucleus
  //           beta1   deformation of first nucleus
  //           Z2      nuclear charge of second nucleus
  //           N2      number of neutrons in second nucleus
  //           beta2   deformation of second nucleus
  //           d       distance of surfaces of the nuclei

  //      G4double Z1,N1,beta1,Z2,N2,beta2,d,ecoul;
  G4double ecoul = 0;
  G4double dtot = 0;
  const G4double r0 = 1.16;

  dtot = r0 * ( std::pow((z1+n1),0.33333) * (1.0+(2.0/3.0)*beta1)
		+ std::pow((z2+n2),0.33333) * (1.0+(2.0/3.0)*beta2) ) + d;
  ecoul = z1 * z2 * 1.44 / dtot;

  // assert(isnan(ecoul) == false);
  return ecoul;
}

void G4Abla::fissionDistri(G4double &a,G4double &z,G4double &e,
			   G4double &a1,G4double &z1,G4double &e1,G4double &v1,
			   G4double &a2,G4double &z2,G4double &e2,G4double &v2)
{
  //  On input: A, Z, E (mass, atomic number and exc. energy of compound nucleus
  //                     before fission)
  //  On output: Ai, Zi, Ei (mass, atomic number and exc. energy of fragment 1 and 2
  //                     after fission)

  //  Additionally calculated but not put in the parameter list:
  //  Kinetic energy of prefragments EkinR1, EkinR2

  //  Translation of SIMFIS18.PLI (KHS, 2.1.2001)

  // This program calculates isotopic distributions of fission fragments
  // with a semiempirical model                      
  // Copy from SIMFIS3, KHS, 8. February 1995        
  // Modifications made by Jose Benlliure and KHS in August 1996
  // Energy counted from lowest barrier (J. Benlliure, KHS 1997)
  // Some bugs corrected (J. Benlliure, KHS 1997)         
  // Version used for thesis S. Steinhaueser (August 1997)
  // (Curvature of LD potential increased by factor of 2!)

  // Weiter veraendert mit der Absicht, eine Version zu erhalten, die
  // derjenigen entspricht, die von J. Benlliure et al.
  // in Nucl. Phys. A 628 (1998) 458 verwendet wurde,
  // allerdings ohne volle Neutronenabdampfung.

  // The excitation energy was calculate now for each fission channel
  // separately. The dissipation from saddle to scission was taken from
  // systematics, the deformation energy at scission considers the shell
  // effects in a simplified way, and the fluctuation is included. 
  // KHS, April 1999 

  // The width in N/Z was carefully adapted to values given by Lang et al.

  // The width and eventually a shift in N/Z (polarization) follows the
  // following rules:                                        

  // The line N/Z following UCD has an angle of std::atan(Zcn/Ncn)
  // to the horizontal axis on a chart of nuclides.
  // (For 238U the angle is 32.2 deg.)

  // The following relations hold: (from Armbruster)
  // 
  // sigma(N) (A=const) = sigma(Z) (A=const)
  // sigma(A) (N=const) = sigma(Z) (N=const)
  // sigma(A) (Z=const) = sigma(N) (Z=const)
  // 
  // From this we get:
  // sigma(Z) (N=const) * N = sigma(N) (Z=const) * Z
  // sigma(A) (Z=const) = sigma(Z) (A=const) * A/Z
  // sigma(N) (Z=const) = sigma(Z) (A=const) * A/Z
  // Z*sigma(N) (Z=const) = N*sigma(Z) (N=const) = A*sigma(Z) (A=const)

  // Excitation energy now calculated above the lowest potential point
  // Inclusion of a distribution of excitation energies             

  // Several modifications, starting from SIMFIS12: KHS November 2000 
  // This version seems to work quite well for 238U.                  
  // The transition from symmetric to asymmetric fission around 226Th 
  // is reasonably well reproduced, although St. I is too strong and St. II
  // is too weak. St. I and St. II are also weakly seen for 208Pb.

  // Extensions for an event generator of fission events (21.11.2000,KHS)

  // Defalt parameters (IPARS) rather carefully adjusted to
  // pre-neutron mass distributions of Vives et al. (238U + n)
  // Die Parameter Fgamma1 und Fgamma2 sind kleiner als die resultierenden
  // Breiten der Massenverteilungen!!!  
  // Fgamma1 und Fgamma2 wurden angepa�, so da�
  // Sigma-A(ST-I) = 3.3, Sigma-A(St-II) = 5.8 (nach Vives) 

  // Parameters of the model carefully adjusted by KHS (2.2.2001) to
  // 238U + 208Pb, 1000 A MeV, Timo Enqvist et al.         


  G4double     n = 0.0;
  G4double     nlight1 = 0.0, nlight2 = 0.0;
  G4double     aheavy1 = 0.0,alight1 = 0.0, aheavy2 = 0.0, alight2 = 0.0;
  G4double     eheavy1 = 0.0, elight1 = 0.0, eheavy2 = 0.0, elight2 = 0.0;
  G4double     zheavy1_shell = 0.0, zheavy2_shell = 0.0;
  G4double     zlight1 = 0.0, zlight2 = 0.0;
  G4double     masscurv = 0.0;
  G4double     sasymm1 = 0.0, sasymm2 = 0.0, ssymm = 0.0, ysum = 0.0, yasymm = 0.0;
  G4double     ssymm_mode1 = 0.0, ssymm_mode2 = 0.0;
  G4double     cz_asymm1_saddle = 0.0, cz_asymm2_saddle = 0.0;
  // Curvature at saddle, modified by ld-potential
  G4double     wzasymm1_saddle, wzasymm2_saddle, wzsymm_saddle  = 0.0;
  G4double     wzasymm1_scission = 0.0, wzasymm2_scission = 0.0, wzsymm_scission = 0.0;
  G4double     wzasymm1 = 0.0, wzasymm2 = 0.0, wzsymm = 0.0;
  G4double     nlight1_eff = 0.0, nlight2_eff = 0.0;
  G4int  imode = 0;
  G4double     rmode = 0.0;
  G4double     z1mean = 0.0, z2mean = 0.0, z1width = 0.0, za1width = 0.0;
  //      G4double     Z1,Z2,N1R,N2R,A1R,A2R,N1,N2,A1,A2;
  G4double     n1r = 0.0, n2r = 0.0, a1r = 0.0, a2r = 0.0, n1 = 0.0, n2 = 0.0;

  G4double     zsymm = 0.0, nsymm = 0.0, asymm = 0.0;
  G4double     n1mean = 0.0, n2mean, n1width;
  G4double     dueff = 0.0;
  // effective shell effect at lowest barrier
  G4double     eld = 0.0;
  // Excitation energy with respect to ld barrier
  G4double     re1 = 0.0, re2 = 0.0, re3 = 0.0;
  G4double     eps1 = 0.0, eps2 = 0.0;
  G4double     n1ucd = 0.0, n2ucd = 0.0, z1ucd = 0.0, z2ucd = 0.0;
  G4double     beta = 0.0, beta1 = 0.0, beta2 = 0.0;

  G4double     dn1_pol = 0.0;
  // shift of most probable neutron number for given Z,
  // according to polarization
  G4int  i_help = 0;

  //   /* Parameters of the semiempirical fission model */
  G4double a_levdens = 0.0;
  //           /* level-density parameter */
  G4double a_levdens_light1 = 0.0, a_levdens_light2 = 0.0;
  G4double a_levdens_heavy1 = 0.0, a_levdens_heavy2 = 0.0;
  const G4double r_null = 1.16;
  //          /* radius parameter */
  G4double epsilon_1_saddle = 0.0, epsilon0_1_saddle = 0.0;
  G4double epsilon_2_saddle = 0.0, epsilon0_2_saddle = 0.0, epsilon_symm_saddle = 0.0;
  G4double epsilon_1_scission = 0.0, epsilon0_1_scission = 0.0;
  G4double epsilon_2_scission = 0.0, epsilon0_2_scission = 0.0;
  G4double epsilon_symm_scission = 0.0;
  //                                   /* modified energy */
  G4double e_eff1_saddle = 0.0, e_eff2_saddle = 0.0;
  G4double epot0_mode1_saddle = 0.0, epot0_mode2_saddle = 0.0, epot0_symm_saddle = 0.0;
  G4double epot_mode1_saddle = 0.0, epot_mode2_saddle = 0.0, epot_symm_saddle = 0.0;
  G4double e_defo = 0.0, e_defo1 = 0.0, e_defo2 = 0.0, e_scission = 0.0, e_asym = 0.0;
  G4double e1exc = 0.0, e2exc = 0.0;
  G4double e1exc_sigma = 0.0, e2exc_sigma = 0.0;
  G4double e1final = 0.0, e2final = 0.0;

  const G4double r0 = 1.16;
  G4double tker = 0.0;
  G4double ekin1 = 0.0, ekin2 = 0.0;
  //      G4double EkinR1,EkinR2,E1,E2,V1,V2;
  G4double ekinr1 = 0.0, ekinr2 = 0.0;
  G4int icz = 0, k = 0;

  //   Input parameters:
  //OMMENT(Nuclear charge number);
  //      G4double Z;
  //OMMENT(Nuclear mass number);
  //      G4double A;
  //OMMENT(Excitation energy above fission barrier);
  //      G4double E;

  //   Model parameters:
  //OMMENT(position of heavy peak valley 1);
  const G4double nheavy1 = 83.0;
  //OMMENT(position of heavy peak valley 2);
  const G4double nheavy2 = 90.0;
  //OMMENT(Shell effect for valley 1);
  const G4double delta_u1_shell = -2.65;
  //        Parameter (Delta_U1_shell = -2)
  //OMMENT(Shell effect for valley 2);
  const G4double delta_u2_shell = -3.8;
  //        Parameter (Delta_U2_shell = -3.2)
  //OMMENT(I: used shell effect);
  G4double delta_u1 = 0.0;
  //omment(I: used shell effect);
  G4double delta_u2 = 0.0;
  //OMMENT(Curvature of asymmetric valley 1);
  const G4double cz_asymm1_shell = 0.7;
  //OMMENT(Curvature of asymmetric valley 2);
  const G4double cz_asymm2_shell = 0.15;
  //OMMENT(Factor for width of distr. valley 1);
  const G4double fwidth_asymm1 = 0.63;
  //OMMENT(Factor for width of distr. valley 2);
  const G4double fwidth_asymm2 = 0.97;
  //       Parameter (CZ_asymm2_scission = 0.12)
  //OMMENT(Parameter x: a = A/x);
  const G4double xlevdens = 12.0;
  //OMMENT(Factor to gamma_heavy1);
  const G4double fgamma1 = 2.0;
  //OMMENT(I: fading of shells (general));
  G4double gamma = 0.0;
  //OMMENT(I: fading of shell 1);
  G4double gamma_heavy1 = 0.0;
  //OMMENT(I: fading of shell 2);
  G4double gamma_heavy2 = 0.0;
  //OMMENT(Zero-point energy at saddle);
  const G4double e_zero_point = 0.5;
  //OMMENT(I: friction from saddle to scission);
  G4double e_saddle_scission = 0.0;
  //OMMENT(Friction factor);
  const G4double friction_factor = 1.0;
  //OMMENT(I: Internal counter for different modes); INIT(0,0,0)
  //      Integer*4 I_MODE(3)
  //OMMENT(I: Yield of symmetric mode);
  G4double ysymm = 0.0;
  //OMMENT(I: Yield of asymmetric mode 1);
  G4double yasymm1 = 0.0;
  //OMMENT(I: Yield of asymmetric mode 2);
  G4double yasymm2 = 0.0;
  //OMMENT(I: Effective position of valley 1);
  G4double nheavy1_eff = 0.0;
  //OMMENT(I: position of heavy peak valley 1);
  G4double zheavy1 = 0.0;
  //omment(I: Effective position of valley 2);
  G4double nheavy2_eff = 0.0;
  //OMMENT(I: position of heavy peak valley 2);
  G4double zheavy2 = 0.0;
  //omment(I: Excitation energy above saddle 1);
  G4double eexc1_saddle = 0.0;
  //omment(I: Excitation energy above saddle 2);
  G4double eexc2_saddle = 0.0;
  //omment(I: Excitation energy above lowest saddle);
  G4double eexc_max = 0.0;
  //omment(I: Effective mass mode 1);
  G4double aheavy1_mean = 0.0;
  //omment(I: Effective mass mode 2);
  G4double aheavy2_mean = 0.0;
  //omment(I: Width of symmetric mode);
  G4double wasymm_saddle = 0.0;
  //OMMENT(I: Width of asymmetric mode 1);
  G4double waheavy1_saddle = 0.0;
  //OMMENT(I: Width of asymmetric mode 2);
  G4double waheavy2_saddle = 0.0;
  //omment(I: Width of symmetric mode);
  G4double wasymm = 0.0;
  //OMMENT(I: Width of asymmetric mode 1);
  G4double waheavy1 = 0.0;
  //OMMENT(I: Width of asymmetric mode 2);
  G4double waheavy2 = 0.0;
  //OMMENT(I: Even-odd effect in Z);
  G4double r_e_o = 0.0, r_e_o_exp = 0.0;
  //OMMENT(I: Curveture of symmetric valley);
  G4double cz_symm = 0.0;
  //OMMENT(I: Curvature of mass distribution for fixed Z);
  G4double cn = 0.0;
  //OMMENT(I: Curvature of Z distribution for fixed A);
  G4double cz = 0.0;
  //OMMENT(Minimum neutron width for constant Z);
  const G4double sigzmin = 1.16;
  //OMMENT(Surface distance of scission configuration);
  const G4double d = 2.0;

  //   /* Charge polarisation from Wagemanns p. 397: */
  //OMMENT(Charge polarisation standard I);
  const G4double cpol1 = 0.65;
  //OMMENT(Charge polarisation standard II);
  const G4double cpol2 = 0.55;
  //OMMENT(=1: Polarisation simult. in N and Z);
  const G4int nzpol = 1;
  //OMMENT(=1: test output, =0: no test output);
  const G4int itest = 0;
      
  //      G4double UMASS, ECOUL, reps1, reps2, rn1_pol;
  G4double reps1 = 0.0, reps2 = 0.0, rn1_pol = 0.0;
  //      Float_t HAZ,GAUSSHAZ;
  G4int kkk = 0;
  //  G4int kkk = 10; // PK
  
  //     I_MODE = 0;

  if(itest == 1) {
    G4cout << " cn mass " << a << G4endl;
    G4cout << " cn charge " << z << G4endl;
    G4cout << " cn energy " << e << G4endl;
  }

  //     /* average Z of asymmetric and symmetric components: */
  n = a - z;  /* neutron number of the fissioning nucleus */

  k = 0;
  icz = 0;
  if ( (std::pow(z,2)/a < 25.0) || (n < nheavy2) || (e > 500.0) ) {
    icz = -1;
    //          GOTO 1002;
    goto milledeux;
  }

  nlight1 = n - nheavy1;
  nlight2 = n - nheavy2;
	
  //    /* Polarisation assumed for standard I and standard II:
  //      Z - Zucd = cpol (for A = const);
  //      from this we get (see Armbruster)
  //      Z - Zucd =  Acn/Ncn * cpol (for N = const)                        */

  zheavy1_shell = ((nheavy1/n) * z) - ((a/n) * cpol1);
  zheavy2_shell = ((nheavy2/n) * z) - ((a/n) * cpol2);

  e_saddle_scission = 
    (-24.0 + 0.02227 * (std::pow(z,2))/(std::pow(a,0.33333)) ) * friction_factor;
    
  //      /* Energy dissipated from saddle to scission                        */
  //      /* F. Rejmund et al., Nucl. Phys. A 678 (2000) 215, fig. 4 b        */
  //      E_saddle_scission = DMAX1(0.,E_saddle_scission);
  if (e_saddle_scission > 0.) {
    e_saddle_scission = e_saddle_scission;
  }
  else {
    e_saddle_scission = 0.;
  }
  //     /* Semiempirical fission model: */

  //    /* Fit to experimental result on curvature of potential at saddle */
  //           /* reference:                                              */
  //    /* IF Z**2/A < 33.15E0 THEN
  //       MassCurv = 30.5438538E0 - 4.00212049E0 * Z**2/A
  //                               + 0.11983384E0 * Z**4 / (A**2) ;
  //     ELSE
  //       MassCurv = 10.E0 ** (7.16993332E0 - 0.26602401E0 * Z**2/A
  //                               + 0.00283802E0 * Z**4 / (A**2)) ;  */
  //  /* New parametrization of T. Enqvist according to Mulgin et al. 1998 */
  if ( (std::pow(z,2))/a < 34.0) {
    masscurv =  std::pow( 10.0,(-1.093364 + 0.082933 * (std::pow(z,2)/a)
			   - 0.0002602 * (std::pow(z,4)/std::pow(a,2))) );
  } else {
    masscurv = std::pow( 10.0,(3.053536 - 0.056477 * (std::pow(z,2)/a)
			  + 0.0002454 * (std::pow(z,4)/std::pow(a,2))) );
  }

  cz_symm = (8.0/std::pow(z,2)) * masscurv;

  if(itest == 1) {
    G4cout << "cz_symmetry= " << cz_symm << G4endl;
  }

  if (cz_symm < 0) {
    icz = -1;
    //          GOTO 1002;
    goto milledeux;
  }

  //  /* proton number in symmetric fission (centre) */
  zsymm  = z/2.0;
  nsymm  = n/2.0;
  asymm = nsymm + zsymm;

  zheavy1 = (cz_symm*zsymm + cz_asymm1_shell*zheavy1_shell)/(cz_symm + cz_asymm1_shell);
  zheavy2 = (cz_symm*zsymm + cz_asymm2_shell*zheavy2_shell)/(cz_symm + cz_asymm2_shell);
  //            /* position of valley due to influence of liquid-drop potential */
  nheavy1_eff = (zheavy1 + (a/n * cpol1))*(n/z);
  nheavy2_eff = (zheavy2 + (a/n * cpol2))*(n/z);
  nlight1_eff = n - nheavy1_eff;
  nlight2_eff = n - nheavy2_eff;
  //  /* proton number of light fragments (centre) */
  zlight1 = z - zheavy1;
  //  /* proton number of light fragments (centre) */
  zlight2 = z - zheavy2;
  aheavy1 = nheavy1_eff + zheavy1;
  aheavy2 = nheavy2_eff + zheavy2;
  aheavy1_mean = aheavy1;
  aheavy2_mean = aheavy2;
  alight1 = nlight1_eff + zlight1;
  alight2 = nlight2_eff + zlight2;

  a_levdens = a / xlevdens;
  a_levdens_heavy1 = aheavy1 / xlevdens;
  a_levdens_heavy2 = aheavy2 / xlevdens;
  a_levdens_light1 = alight1 / xlevdens;
  a_levdens_light2 = alight2 / xlevdens;
  gamma = a_levdens / (0.4 * (std::pow(a,1.3333)) );
  gamma_heavy1 = ( a_levdens_heavy1 / (0.4 * (std::pow(aheavy1,1.3333)) ) ) * fgamma1;
  gamma_heavy2 = a_levdens_heavy2 / (0.4 * (std::pow(aheavy2,1.3333)) );

  cz_asymm1_saddle = cz_asymm1_shell + cz_symm;
  cz_asymm2_saddle = cz_asymm2_shell + cz_symm;
	
  // Up to here: Ok! Checked CS 10/10/05      	   

  cn = umass(zsymm,(nsymm+1.),0.0) + umass(zsymm,(nsymm-1.),0.0)
    + 1.44 * (std::pow(zsymm,2))/
    ( (std::pow(r_null,2)) * 
      ( std::pow((asymm+1.0),0.33333) + std::pow((asymm-1.0),0.33333) ) *
      ( std::pow((asymm+1.0),0.33333) + std::pow((asymm-1.0),0.33333) ) )
    - 2.0 * umass(zsymm,nsymm,0.0)
    - 1.44 * (std::pow(zsymm,2))/
    ( ( 2.0 * r_null * (std::pow(asymm,0.33333)) ) * 
      ( 2.0 * r_null * (std::pow(asymm,0.33333)) ) );
	
  // /* shell effect in valley of mode 1 */
  delta_u1 = delta_u1_shell + (std::pow((zheavy1_shell-zheavy1),2))*cz_asymm1_shell;
  // /* shell effect in valley of mode 2 */
  delta_u2 = delta_u2_shell + (std::pow((zheavy2_shell-zheavy2),2))*cz_asymm2_shell;

  //     /* liquid drop energies
  //        at the centres of the different shell effects
  //        with respect to liquid drop at symmetry: */
  epot0_mode1_saddle = (std::pow((zheavy1-zsymm),2)) * cz_symm;
  epot0_mode2_saddle = (std::pow((zheavy2-zsymm),2)) * cz_symm;
  epot0_symm_saddle = 0.0;
      
  if (itest == 1) {
    G4cout << "check zheavy1 = " << zheavy1  << G4endl;
    G4cout << "check zheavy2 = " << zheavy2  << G4endl;
    G4cout << "check zsymm = " << zsymm  << G4endl;
    G4cout << "check czsymm = " << cz_symm  << G4endl;
    G4cout << "check epot0_mode1_saddle = " << epot0_mode1_saddle  << G4endl;
    G4cout << "check epot0_mode2_saddle = " << epot0_mode2_saddle  << G4endl;
    G4cout << "check epot0_symm_saddle = " << epot0_symm_saddle  << G4endl;
    G4cout << "delta_u1 = " << delta_u1 << G4endl;
    G4cout << "delta_u2 = " << delta_u2 << G4endl;
  }
      
  //     /* energies including shell effects
  //        at the centres of the different shell effects
  //        with respect to liquid drop at symmetry: */
  epot_mode1_saddle = epot0_mode1_saddle + delta_u1;
  epot_mode2_saddle = epot0_mode2_saddle + delta_u2;
  epot_symm_saddle = epot0_symm_saddle;
  if (itest == 1) {
    G4cout << "check epot_mode1_saddle = " << epot_mode1_saddle  << G4endl;
    G4cout << "check epot_mode2_saddle = " << epot_mode2_saddle  << G4endl;
    G4cout << "check epot_symm_saddle = " << epot_symm_saddle  << G4endl;
  }

  //     /* Minimum of potential with respect to ld potential at symmetry */
  dueff = min(epot_mode1_saddle,epot_mode2_saddle);
  dueff = min(dueff,epot_symm_saddle);
  dueff = dueff - epot_symm_saddle;

  eld = e + dueff + e_zero_point;
      
  if (itest == 1) {
    G4cout << "check dueff = " << dueff  << G4endl;
    G4cout << "check e = " << e  << G4endl;
    G4cout << "check e_zero_point = " << e_zero_point  << G4endl;
    G4cout << "check eld = " << eld  << G4endl;
  }
  // Up to here: Ok! Checked CS 10/10/05 
     
  //          /* E = energy above lowest effective barrier */
  //          /* Eld = energy above liquid-drop barrier */

  //     /* Due to this treatment the energy E on input means the excitation  */
  //     /* energy above the lowest saddle.                                   */

  //  /* These energies are not used */
  eheavy1 = e * aheavy1 / a;
  eheavy2 = e * aheavy2 / a;
  elight1 = e * alight1 / a;
  elight2 = e * alight2 / a;

  epsilon0_1_saddle = eld - e_zero_point - epot0_mode1_saddle;
  //            /* excitation energy at saddle mode 1 without shell effect */
  epsilon0_2_saddle = eld - e_zero_point - epot0_mode2_saddle;
  //            /* excitation energy at saddle mode 2 without shell effect */

  epsilon_1_saddle = eld - e_zero_point - epot_mode1_saddle;
  //            /* excitation energy at saddle mode 1 with shell effect */
  epsilon_2_saddle = eld - e_zero_point - epot_mode2_saddle;
  //            /* excitation energy at saddle mode 2 with shell effect */
  epsilon_symm_saddle = eld - e_zero_point - epot_symm_saddle;

  //  /* global parameters */
  eexc1_saddle = epsilon_1_saddle;
  eexc2_saddle = epsilon_2_saddle;
  eexc_max = max(eexc1_saddle,eexc2_saddle);
  eexc_max = max(eexc_max,eld);
       
  //         /* EEXC_MAX is energy above the lowest saddle */


  epsilon0_1_scission = eld + e_saddle_scission - epot0_mode1_saddle;
  //                    /* excitation energy without shell effect */
  epsilon0_2_scission = eld + e_saddle_scission - epot0_mode2_saddle;
  //                    /* excitation energy without shell effect */

  epsilon_1_scission = eld + e_saddle_scission - epot_mode1_saddle;
  //                    /* excitation energy at scission */
  epsilon_2_scission = eld+ e_saddle_scission - epot_mode2_saddle;
  //                    /* excitation energy at scission */
  epsilon_symm_scission = eld + e_saddle_scission - epot_symm_saddle;
  //           /* excitation energy of symmetric fragment at scission */

  //     /* Calculate widhts at the saddle:                */

  e_eff1_saddle = epsilon0_1_saddle - delta_u1 * (std::exp((-epsilon_1_saddle*gamma)));
   
  if (e_eff1_saddle > 0.0) {
    wzasymm1_saddle = std::sqrt( (0.5 * 
			     (std::sqrt(1.0/a_levdens*e_eff1_saddle)) /
			     (cz_asymm1_shell * std::exp((-epsilon_1_saddle*gamma)) + cz_symm) ) );
  } 
  else {
    wzasymm1_saddle = 1.0;
  }

  e_eff2_saddle = epsilon0_2_saddle - delta_u2 * (std::exp((-epsilon_2_saddle*gamma)));
  if (e_eff2_saddle > 0.0) {
    wzasymm2_saddle = std::sqrt( (0.5 * 
			     (std::sqrt(1.0/a_levdens*e_eff2_saddle)) /
			     (cz_asymm2_shell * std::exp((-epsilon_2_saddle*gamma)) + cz_symm) ) );
  } 
  else {
    wzasymm2_saddle = 1.0;
  }

  if (eld > e_zero_point) {
    if ( (eld + epsilon_symm_saddle) < 0.0)  {
      G4cout << "<e> eld + epsilon_symm_saddle < 0" << G4endl;
    }
    wzsymm_saddle = std::sqrt( (0.5 * 
			   (std::sqrt(1.0/a_levdens*(eld+epsilon_symm_saddle))) / cz_symm ) );
  } else {
    wzsymm_saddle = 1.0;
  }

  if (itest == 1) {
    G4cout << "wz1(saddle) = " << wzasymm1_saddle << G4endl;
    G4cout << "wz2(saddle) = " << wzasymm2_saddle << G4endl;
    G4cout << "wzsymm(saddle) = " << wzsymm_saddle << G4endl;
  }
     
  //     /* Calculate widhts at the scission point: */
  //     /* fits of ref. Beizin 1991 (Plots brought to GSI by Sergei Zhdanov) */

  wzsymm_scission = wzsymm_saddle;

  if (e_saddle_scission == 0.0) {

    wzasymm1_scission = wzasymm1_saddle;
    wzasymm2_scission = wzasymm2_saddle;

  } 
  else {

    if (nheavy1_eff > 75.0) {
      wzasymm1_scission = (std::sqrt(21.0)) * z/a;
      wzasymm2_scission = (std::sqrt (max( (70.0-28.0)/3.0*(z*z/a-35.0)+28.,0.0 )) ) * z/a;
    } 
    else {
      wzasymm1_scission = wzasymm1_saddle;
      wzasymm2_scission = wzasymm2_saddle;
    }

  }

  wzasymm1_scission = max(wzasymm1_scission,wzasymm1_saddle);
  wzasymm2_scission = max(wzasymm2_scission,wzasymm2_saddle);

  wzasymm1 = wzasymm1_scission * fwidth_asymm1;
  wzasymm2 = wzasymm2_scission * fwidth_asymm2;
  wzsymm = wzsymm_scission;

  /*      if (ITEST == 1) {
	  G4cout << "WZ1(scission) = " << WZasymm1_scission << G4endl;
	  G4cout << "WZ2(scission) = " << WZasymm2_scission << G4endl;
	  G4cout << "WZsymm(scission) = " << WZsymm_scission << G4endl;
	  }
	  if (ITEST == 1) {
	  G4cout << "WZ1(scission) final= " << WZasymm1 << G4endl;
	  G4cout << "WZ2(scission) final= " << WZasymm2 << G4endl;
	  G4cout << "WZsymm(scission) final= " << WZsymm << G4endl;
	  } */
      
  wasymm = wzsymm * a/z;
  waheavy1 = wzasymm1 * a/z;
  waheavy2 = wzasymm2 * a/z;

  wasymm_saddle = wzsymm_saddle * a/z;
  waheavy1_saddle = wzasymm1_saddle * a/z;
  waheavy2_saddle = wzasymm2_saddle * a/z;

  if (itest == 1) {
    G4cout << "wasymm = " << wzsymm << G4endl;
    G4cout << "waheavy1 = " << waheavy1 << G4endl;
    G4cout << "waheavy2 = " << waheavy2 << G4endl;
  }
  // Up to here: Ok! Checked CS 11/10/05
            
  if ( (epsilon0_1_saddle - delta_u1*std::exp((-epsilon_1_saddle*gamma_heavy1))) < 0.0) {
    sasymm1 = -10.0;
  } 
  else {
    sasymm1 = 2.0 * std::sqrt( a_levdens * (epsilon0_1_saddle - 
				       delta_u1*(std::exp((-epsilon_1_saddle*gamma_heavy1))) ) );
  }

  if ( (epsilon0_2_saddle - delta_u2*std::exp((-epsilon_2_saddle*gamma_heavy2))) < 0.0) {
    sasymm2 = -10.0;
  } 
  else {
    sasymm2 = 2.0 * std::sqrt( a_levdens * (epsilon0_2_saddle - 
				       delta_u2*(std::exp((-epsilon_2_saddle*gamma_heavy2))) ) );
  }
              
  if (epsilon_symm_saddle > 0.0) {
    ssymm = 2.0 * std::sqrt( a_levdens*(epsilon_symm_saddle) );
  } 
  else {
    ssymm = -10.0;
  }
      
  if (ssymm > -10.0) {
    ysymm = 1.0;

    if (epsilon0_1_saddle < 0.0) {
      //  /* low energy */
      yasymm1 = std::exp((sasymm1-ssymm)) * wzasymm1_saddle / wzsymm_saddle * 2.0;
      //           /* factor of 2 for symmetry classes */
    } 
    else {
      //        /* high energy */
      ssymm_mode1 = 2.0 * std::sqrt( a_levdens*(epsilon0_1_saddle) );
      yasymm1 = ( std::exp((sasymm1-ssymm)) - std::exp((ssymm_mode1 - ssymm)) )  
	* wzasymm1_saddle / wzsymm_saddle * 2.0;
    }

    if (epsilon0_2_saddle < 0.0) {
      //  /* low energy */
      yasymm2 = std::exp((sasymm2-ssymm)) * wzasymm2_saddle / wzsymm_saddle * 2.0;
      //           /* factor of 2 for symmetry classes */
    } 
    else {
      //        /* high energy */
      ssymm_mode2 = 2.0 * std::sqrt( a_levdens*(epsilon0_2_saddle) );
      yasymm2 = ( std::exp((sasymm2-ssymm)) - std::exp((ssymm_mode2 - ssymm)) )  
	* wzasymm2_saddle / wzsymm_saddle * 2.0;
    }       
    //                            /* difference in the exponent in order */
    //                            /* to avoid numerical overflow         */

  } 
  else {
    if ( (sasymm1 > -10.0) && (sasymm2 > -10.0) ) {
      ysymm = 0.0;
      yasymm1 = std::exp(sasymm1) * wzasymm1_saddle * 2.0;
      yasymm2 = std::exp(sasymm2) * wzasymm2_saddle * 2.0;
    }
  }
        
  //  /* normalize */
  ysum = ysymm + yasymm1 + yasymm2;
  if (ysum > 0.0) {
    ysymm = ysymm / ysum;
    yasymm1 = yasymm1 / ysum;
    yasymm2 = yasymm2 / ysum;
    yasymm = yasymm1 + yasymm2;
  } 
  else {
    ysymm = 0.0;
    yasymm1 = 0.0;
    yasymm2 = 0.0;
    //        /* search minimum threshold and attribute all events to this mode */
    if ( (epsilon_symm_saddle < epsilon_1_saddle) && (epsilon_symm_saddle < epsilon_2_saddle) ) {
      ysymm = 1.0;
    } 
    else {
      if (epsilon_1_saddle < epsilon_2_saddle) {
	yasymm1 = 1.0;
      } 
      else {
	yasymm2 = 1.0;
      }
    }
  }

  if (itest == 1) {
    G4cout << "ysymm normalized= " << ysymm  << G4endl;
    G4cout << "yasymm1 normalized= " << yasymm1  << G4endl;
    G4cout << "yasymm2 normalized= " << yasymm2  << G4endl;
  }
  // Up to here: Ok! Ckecked CS 11/10/05      
      
  //      /* even-odd effect */
  //      /* simple parametrization KHS, Nov. 2000. From Rejmund et al. */
  if ((int)(z) % 2 == 0) {
    r_e_o_exp = -0.017 * (e_saddle_scission + eld) * (e_saddle_scission + eld);
    if ( r_e_o_exp < -307.0) {
      r_e_o_exp = -307.0;
      r_e_o = std::pow(10.0,r_e_o_exp);
    }
    else {
      r_e_o = std::pow(10.0,r_e_o_exp);
    }
  } 
  else {
    r_e_o = 0.0;
  }
    
  //      $LOOP;    /* event loop */
  //     I_COUNT = I_COUNT + 1;

  //     /* random decision: symmetric or asymmetric */
  //     /* IMODE = 1 means asymmetric fission, mode 1,
  //        IMODE = 2 means asymmetric fission, mode 2,
  //        IMODE = 3 means symmetric  */
  //      RMODE = dble(HAZ(k));
  //      rmode = rnd.rndm();  

  // Safety check added to make sure we always select well defined
  // fission mode.
  do {
    rmode = haz(k);
    // Cast for test CS 11/10/05
    //      RMODE = 0.54;    
    //  rmode = 0.54;
    if (rmode < yasymm1) {
      imode = 1;
    } else if ( (rmode > yasymm1) && (rmode < (yasymm1+yasymm2)) ) {
      imode = 2;
    } else if ( (rmode > yasymm1) && (rmode > (yasymm1+yasymm2)) ) {
      imode = 3;
    }
  } while(imode == 0);

  //     /* determine parameters of the Z distribution */
  // force imode (for testing, PK)
  // imode = 3;
  if (imode == 1) {
    z1mean = zheavy1;
    z1width = wzasymm1;
  }
  if (imode == 2) {
    z1mean = zheavy2;
    z1width = wzasymm2;
  }
  if (imode == 3) {
    z1mean = zsymm;
    z1width = wzsymm;
  }

  if (itest == 1) {
    G4cout << "nbre aleatoire tire " << rmode << G4endl;
    G4cout << "fission mode " << imode << G4endl;
    G4cout << "z1mean= " << z1mean << G4endl;
    G4cout << "z1width= " << z1width << G4endl;
  }
		      
  //     /* random decision: Z1 and Z2 at scission: */
  z1 = 1.0;
  z2 = 1.0;
  while  ( (z1<5.0) || (z2<5.0) ) {
    //         Z1 = dble(GAUSSHAZ(K,sngl(Z1mean),sngl(Z1width)));
    //	 z1 = rnd.gaus(z1mean,z1width);
    z1 = gausshaz(k, z1mean, z1width);
    z2 = z - z1;
  }
  if (itest == 1) {
    G4cout << "ff charge sample " << G4endl;
    G4cout << "z1 =  " << z1 << G4endl;
    G4cout << "z2 = " << z2 << G4endl;
  }

  //     CALL EVEN_ODD(Z1,R_E_O,I_HELP);
  //         /* Integer proton number with even-odd effect */
  //     Z1 = REAL(I_HELP)
  //      /* Z1 = INT(Z1+0.5E0); */
  z2 = z - z1;

  //     /* average N of both fragments: */
  if (imode == 1) {
    n1mean = (z1 + cpol1 * a/n) * n/z;
  }
  if (imode == 2) { 
    n1mean = (z1 + cpol2 * a/n) * n/z;
  }
  /*       CASE(99)   ! only for testing;
	   N1UCD = Z1 * N/Z;
	   N2UCD = Z2 * N/Z;
	   re1 = UMASS(Z1,N1UCD,0.6) +;
	   &         UMASS(Z2,N2UCD,0.6) +;
	   &         ECOUL(Z1,N1UCD,0.6,Z2,N2UCD,0.6,d);
	   re2 = UMASS(Z1,N1UCD+1.,0.6) +;
	   &         UMASS(Z2,N2UCD-1.,0.6) +;
	   &         ECOUL(Z1,N1UCD+1.,0.6,Z2,N2UCD-1.,0.6,d);
	   re3 = UMASS(Z1,N1UCD+2.,0.6) +;
	   &         UMASS(Z2,N2UCD-2.,0.6) +;
	   &         ECOUL(Z1,N1UCD+2.,0.6,Z2,N2UCD-2.,0.6,d);
	   eps2 = (re1-2.0*re2+re3) / 2.0;
	   eps1 = re2 - re1 - eps2;
	   DN1_POL = - eps1 / (2.0 * eps2);
	   N1mean = N1UCD + DN1_POL; */
  if (imode == 3) {
    n1ucd = z1 * n/z;
    n2ucd = z2 * n/z;
    re1 = umass(z1,n1ucd,0.6) + umass(z2,n2ucd,0.6) + ecoul(z1,n1ucd,0.6,z2,n2ucd,0.6,d);
    re2 = umass(z1,n1ucd+1.,0.6) + umass(z2,n2ucd-1.,0.6) + ecoul(z1,n1ucd+1.,0.6,z2,n2ucd-1.,0.6,d);
    re3 = umass(z1,n1ucd+2.,0.6) + umass(z2,n2ucd-2.,0.6) + ecoul(z1,n1ucd+2.,0.6,z2,n2ucd-2.,0.6,d);
    eps2 = (re1-2.0*re2+re3) / 2.0;
    eps1 = re2 - re1 - eps2;
    dn1_pol = - eps1 / (2.0 * eps2);
    n1mean = n1ucd + dn1_pol;
  }
  // all fission modes features have been checked CS 11/10/05
  n2mean = n - n1mean;
  z2mean = z - z1mean;
      
  //   /* Excitation energies */
  //     /* formulated in energies in close consistency with the fission model */

  //     /* E_defo = UMASS(Z*0.5E0,N*0.5E0,0.6E0) -
  //                 UMASS(Z*0.5E0,N*0.5E0,0);   */
  //       /* calculates the deformation energy of the liquid drop for
  //          deformation beta = 0.6 which is most probable at scission */

  //    /* N1R and N2R provisionaly taken without fluctuations in
  //       polarisation:                                              */
  n1r = n1mean;
  n2r = n2mean;
  a1r = n1r + z1;
  a2r = n2r + z2;

  if (imode == 1) { /* N = 82 */;
    //! /* Eexc at scission */
    e_scission = max(epsilon_1_scission,1.0);
    if (n1mean > (n * 0.5) ) {
      //! /* 1. fragment is spherical */
      beta1 = 0.0;
      beta2 = 0.6;
      e1exc = epsilon_1_scission * a1r / a;
      e_defo = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
      e2exc = epsilon_1_scission * a2r / a + e_defo;
    }
    else {
      //!                       /* 2. fragment is spherical */
      beta1 = 0.6;
      beta2 = 0.0;
      e_defo = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
      e1exc = epsilon_1_scission * a1r / a + e_defo;
      e2exc = epsilon_1_scission * a2r / a;
    }
  }
	   
  if (imode == 2) { 
    //! /* N appr. 86 */
    e_scission = max(epsilon_2_scission,1.0);
    if (n1mean >  (n * 0.5) ) { 
      //! /* 2. fragment is spherical */
      beta1 = (n1r - nheavy2) * 0.034 + 0.3;
      e_defo = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
      e1exc = epsilon_2_scission * a1r / a + e_defo;
      beta2 = 0.6 - beta1;
      e_defo = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
      e2exc = epsilon_2_scission * a2r / a + e_defo;
    }
    else {
      //!                      /* 1. fragment is spherical */
      beta2 = (n2r - nheavy2) * 0.034 + 0.3;
      e_defo = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
      e2exc = epsilon_2_scission * a2r / a + e_defo;
      beta1 = 0.6 - beta2;
      e_defo = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
      e1exc = epsilon_2_scission * a1r / a + e_defo;
    }
  }
         
  if (imode == 3) { 
    // ! /* Symmetric fission channel */

    //             /* the fit function for beta is the deformation for
    //                optimum energy at the scission point, d = 2 */
    //             /* beta  : deformation of symmetric fragments */
    //             /* beta1 : deformation of first fragment */
    //             /* beta2 : deformation of second fragment */
    beta =  0.177963 + 0.0153241 * zsymm - 0.000162037 * zsymm*zsymm;
    beta1 = 0.177963 + 0.0153241 * z1 - 0.000162037 * z1*z1;
    //            beta1 = 0.6
    e_defo1 = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
    beta2 = 0.177963 + 0.0153241 * z2 - 0.000162037 * z2*z2;
    //            beta2 = 0.6
    e_defo2 = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
    e_asym = umass(z1 , n1r, beta1) + umass(z2, n2r ,beta2)
      + ecoul(z1,n1r,beta1,z2,n2r,beta2,2.0)
      - 2.0 * umass(zsymm,nsymm,beta)
      - ecoul(zsymm,nsymm,beta,zsymm,nsymm,beta,2.0);
    //            E_asym = CZ_symm * (Z1 - Zsymm)**2
    e_scission = max((epsilon_symm_scission - e_asym),1.0);
    //         /*  $LIST(Z1,N1R,Z2,N2R,E_asym,E_scission); */
    e1exc = e_scission * a1r / a + e_defo1;
    e2exc = e_scission * a2r / a + e_defo2;
  }
  // Energies checked for all the modes CS 11/10/05
	  
  //    /* random decision: N1R and N2R at scission, before evaporation: */
  //    /* CN =       UMASS(Zsymm , Nsymm + 1.E0,0) +
  //                UMASS(Zsymm, Nsymm - 1.E0,0)
  //                 + 1.44E0 * (Zsymm)**2 /
  //                 (r_null**2 * ((Asymm+1)**1/3 + (Asymm-1)**1/3)**2 )
  //                 - 2.E0 * UMASS(Zsymm,Nsymm,0)
  //                 - 1.44E0 * (Zsymm)**2 / (r_null * 2.E0 * (Asymm)**1/3)**2; */


  //    /* N1width = std::sqrt(0.5E0 * std::sqrt(1.E0/A_levdens*(Eld+E_saddle_scission)) / CN); */
  //    /* 8. 9. 1998: KHS (see also consideration in the first comment block)
  //       sigma_N(Z=const) = A/Z  * sigma_Z(A=const)
  //       sigma_Z(A=const) = 0.4 to 0.5  (from Lang paper Nucl Phys. A345 (1980) 34)
  //       sigma_N(Z=const) = 0.45 * A/Z  (= 1.16 for 238U)
  //       therefore: SIGZMIN = 1.16                                              */

  if ( (imode == 1) || (imode == 2) ) {
    cn=(umass(z1,n1mean+1.,beta1) + umass(z1,n1mean-1.,beta1)
	+ umass(z2,n2mean+1.,beta2) + umass(z2,n2mean-1.,beta2)
	+ ecoul(z1,n1mean+1.,beta1,z2,n2mean-1.,beta2,2.0)
	+ ecoul(z1,n1mean-1.,beta1,z2,n2mean+1.,beta2,2.0)
	- 2.0 * ecoul(z1,n1mean,beta1,z2,n2mean,beta2,2.0)
	- 2.0 * umass(z1, n1mean, beta1)
	- 2.0 * umass(z2, n2mean, beta2) ) * 0.5;
    //          /* Coulomb energy neglected for the moment! */
    //           IF (E_scission.lt.0.) Then
    //             write(6,*)'<E> E_scission < 0, MODE 1,2'
    //           ENDIF
    //           IF (CN.lt.0.) Then
    //             write(6,*)'CN < 0, MODE 1,2'
    //           ENDIF
    n1width=std::sqrt( (0.5 * (std::sqrt(1.0/a_levdens*(e_scission)))/cn) );
    n1width=max(n1width, sigzmin);

    //          /* random decision: N1R and N2R at scission, before evaporation: */
    n1r = 1.0;
    n2r = 1.0;
    while  ( (n1r<5.0) || (n2r<5.0) ) {
      //             n1r = dble(gausshaz(k,sngl(n1mean),sngl(n1width)));
      //	   	n1r = rnd.gaus(n1mean,n1width);
      n1r = gausshaz(k, n1mean, n1width);
      n2r = n - n1r;
    }	   
    //           N1R = GAUSSHAZ(K,N1mean,N1width)
    if (itest == 1) {
      G4cout << "after neutron sample " << n1r << G4endl;
    }	   
    n1r = (float)( (int)((n1r+0.5)) );
    n2r = n - n1r;
 
    even_odd(z1,r_e_o,i_help);
    //         /*  proton number with even-odd effect */
    z1 = (float)(i_help);
    z2 = z - z1;

    a1r = z1 + n1r;
    a2r = z2 + n2r;
  }
	
  if (imode == 3) {
    //!  /* When(3) */
    if (nzpol > 0.0) {
      //           /* We treat a simultaneous split in Z and N to determine polarisation */
      cz = ( umass(z1-1., n1mean+1.,beta1)
	     +  umass(z2+1., n2mean-1.,beta1) 
	     +  umass(z1+1., n1mean-1.,beta2)
	     +  umass(z2 - 1., n2mean + 1.,beta2)
	     +  ecoul(z1-1.,n1mean+1.,beta1,z2+1.,n2mean-1.,beta2,2.0)
	     +  ecoul(z1+1.,n1mean-1.,beta1,z2-1.,n2mean+1.,beta2,2.0)
	     -  2.0 * ecoul(z1,n1mean,beta1,z2,n2mean,beta2,2.0)
	     -  2.0 * umass(z1, n1mean,beta1)
	     -  2.0 * umass(z2, n2mean,beta2) ) * 0.5;
      //           IF (E_scission.lt.0.) Then
      //             write(6,*) '<E> E_scission < 0, MODE 1,2'
      //           ENDIF
      //           IF (CZ.lt.0.) Then
      //             write(6,*) 'CZ < 0, MODE 1,2'
      //           ENDIF
      za1width=std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*(e_scission)) / cz) );
      za1width=std::sqrt( (max((za1width*za1width-(1.0/12.0)),0.1)) );
      //                        /* Check the value of 0.1 ! */
      //                        /* Shephard correction */
      a1r = z1 + n1mean;
      a1r = (float)((int)((a1r+0.5)));
      a2r = a - a1r;
      //           /* A1R and A2R are integer numbers now */
      //        /* $LIST(A1R,A2R,ZA1WIDTH); */

      n1ucd = n/a * a1r;
      n2ucd = n/a * a2r;
      z1ucd = z/a * a1r;
      z2ucd = z/a * a2r;

      re1 = umass(z1ucd-1.,n1ucd+1.,beta1) + umass(z2ucd+1.,n2ucd-1.,beta2)
	+ ecoul(z1ucd-1.,n1ucd+1.,beta1,z2ucd+1.,n2ucd-1.,beta2,d);
      re2 = umass(z1ucd,n1ucd,beta1) + umass(z2ucd,n2ucd,beta2)
	+ ecoul(z1ucd,n1ucd,beta1,z2ucd,n2ucd,beta2,d);
      re3 = umass(z1ucd+1.,n1ucd-1.,beta1) + umass(z2ucd-1.,n2ucd+1.,beta2) +
	+ ecoul(z1ucd+1.,n1ucd-1.,beta1,z2ucd-1.,n2ucd+1.,beta2,d);
		 
      eps2 = (re1-2.0*re2+re3) / 2.0;
      eps1 = (re3 - re1)/2.0;
      dn1_pol = - eps1 / (2.0 * eps2);
      z1 = z1ucd + dn1_pol;
      if (itest == 1) {
	G4cout << "before proton sample " << z1 << G4endl;
      }  
      //           Z1 = dble(GAUSSHAZ(k,sngl(Z1),sngl(ZA1width)));
      //	   z1 = rnd.gaus(z1,za1width);
      z1 = gausshaz(k, z1, za1width);
      if (itest == 1) {
	G4cout << "after proton sample " << z1 << G4endl;
      }	   
      even_odd(z1,r_e_o,i_help);
      //         /* proton number with even-odd effect */
      z1 = (float)(i_help);
      z2 = (float)((int)( (z - z1 + 0.5)) );

      n1r = a1r - z1;
      n2r = n - n1r;
    } 
    else {
      //           /* First division of protons, then adjustment of neutrons */
      cn = ( umass(z1, n1mean+1.,beta1) + umass(z1, n1mean-1., beta1)
	     + umass(z2, n2mean+1.,beta2) + umass(z2, n2mean-1., beta2)
	     + ecoul(z1,n1mean+1.,beta1,z2,n2mean-1.,beta2,2.0)
	     + ecoul(z1,n1mean-1.,beta1,z2,n2mean+1.,beta2,2.0)
	     - 2.0 * ecoul(z1,n1mean,beta1,z2,n2mean,beta2,2.0)
	     - 2.0 * umass(z1, n1mean, 0.6)
	     - 2.0 * umass(z2, n2mean, 0.6) ) * 0.5;
      //          /* Coulomb energy neglected for the moment! */
      //           IF (E_scission.lt.0.) Then
      //             write(6,*) '<E> E_scission < 0, MODE 1,2'
      //           Endif
      //           IF (CN.lt.0.) Then
      //             write(6,*) 'CN < 0, MODE 1,2'
      //           Endif
      n1width=std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*(e_scission)) / cn) );
      n1width=max(n1width, sigzmin);

      //          /* random decision: N1R and N2R at scission, before evaporation: */
      //           N1R = dble(GAUSSHAZ(k,sngl(N1mean),sngl(N1width)));
      //	   n1r = rnd.gaus(n1mean,n1width);
      n1r = gausshaz(k, n1mean, n1width);
      n1r = (float)( (int)((n1r+0.5)) );
      n2r = n - n1r;

      even_odd(z1,r_e_o,i_help);
      //         /* Integer proton number with even-odd effect */
      z1 = (float)(i_help);
      z2 = z - z1;

      a1r = z1 + n1r;
      a2r = z2 + n2r;
          
    }
  }

  if (itest == 1) {
    G4cout << "remid imode = " << imode << G4endl;
    G4cout << "n1width =  " << n1width << G4endl;
    G4cout << "n1r = " << n1r << G4endl;
    G4cout << "a1r = " << a1r << G4endl;
    G4cout << "n2r = " << n2r << G4endl;
    G4cout << "a2r = " << a2r << G4endl;
  }
  // Up to here: checked CS 11/10/05	
	
  //      /* Extracted from Lang et al. Nucl. Phys. A 345 (1980) 34 */
  e1exc_sigma = 5.5;
  e2exc_sigma = 5.5;

 neufcentquatrevingtsept:
  //       E1final = dble(Gausshaz(k,sngl(E1exc),sngl(E1exc_sigma)));
  //       E2final = dble(Gausshaz(k,sngl(E2exc),sngl(E2exc_sigma)));
  //        e1final = rnd.gaus(e1exc,e1exc_sigma);
  //        e2final = rnd.gaus(e2exc,e2exc_sigma);
  e1final = gausshaz(k, e1exc, e1exc_sigma);
  e2final = gausshaz(k, e2exc, e2exc_sigma);
  if ( (e1final < 0.0) || (e2final < 0.0) ) goto neufcentquatrevingtsept;
  if (itest == 1) {
    G4cout << "sampled exc 1 " << e1final << G4endl;
    G4cout << "sampled exc 2 " << e2final << G4endl;
  }

  //      /* OUTPUT QUANTITIES OF THE EVENT GENERATOR:        */

  //      /* Quantities before neutron evaporation            */

  //      /* Neutron number of prefragments: N1R and N2R      */
  //      /* Atomic number of fragments: Z1 and Z2            */
  //      /* Kinetic energy of fragments: EkinR1, EkinR2      *7

  //      /* Quantities after neutron evaporation:            */

  //      /* Neutron number of fragments: N1 and N2           */
  //      /* Mass number of fragments: A1 and A2              */
  //      /* Atomic number of fragments: Z1 and Z2            */
  //      /* Number of evaporated neutrons: N1R-N1 and N2R-N2 */
  //      /* Kinetic energy of fragments: EkinR1*A1/A1R and
  //                                      EkinR2*A2/A2R       */

  n1 = n1r;
  n2 = n2r;
  a1 = n1 + z1;
  a2 = n2 + z2;
  e1 = e1final;
  e2 = e2final;

  //       /* Pre-neutron-emission total kinetic energy: */
  tker = (z1 * z2 * 1.44) /
    ( r0 * std::pow(a1,0.33333) * (1.0 + 2.0/3.0 * beta1) +
      r0 * std::pow(a2,0.33333) * (1.0 + 2.0/3.0 * beta2) + 2.0 );
  //       /* Pre-neutron-emission kinetic energy of 1. fragment: */
  ekinr1 = tker * a2 / a;
  //       /* Pre-neutron-emission kinetic energy of 2. fragment: */
  ekinr2 = tker * a1 / a;

  v1 = std::sqrt( (ekinr1/a1) ) * 1.3887;
  v2 = std::sqrt( (ekinr2/a2) ) * 1.3887;

  if (itest == 1) {
    G4cout << "ekinr1 " << ekinr1 << G4endl;
    G4cout << "ekinr2 " << ekinr2 << G4endl;
  }

 milledeux:       
  //**************************
  //*** only symmetric fission
  //**************************
  // Symmetric fission: Ok! Checked CS 10/10/05
  if ( (icz == -1) || (a1 < 0.0) || (a2 < 0.0) ) {
    //           IF (z.eq.92) THEN
    //              write(6,*)'symmetric fission'
    //              write(6,*)'Z,A,E,A1,A2,icz,Atot',Z,A,E,A1,A2,icz,Atot
    //           END IF

    if (itest == 1) {
      G4cout << "milledeux: liquid-drop option "  << G4endl;
    }

    n = a-z;
    //  proton number in symmetric fission (centre) *
    zsymm  = z / 2.0;
    nsymm  = n / 2.0;
    asymm = nsymm + zsymm;

    a_levdens = a / xlevdens;

    masscurv = 2.0;
    cz_symm = 8.0 / std::pow(z,2) * masscurv;

    wzsymm = std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*e) / cz_symm) ) ;

    if (itest == 1) {
      G4cout << " symmetric high energy fission " << G4endl;
      G4cout << "wzsymm " << wzsymm << G4endl;
    }

    z1mean = zsymm;
    z1width = wzsymm;

    // random decision: Z1 and Z2 at scission: */
    z1 = 1.0;
    z2 = 1.0;
    while  ( (z1 < 5.0) || (z2 < 5.0) ) {
      //           z1 = dble(gausshaz(kkk,sngl(z1mean),sngl(z1width)));
      //	   z1 = rnd.gaus(z1mean,z1width);
      z1 = gausshaz(kkk, z1mean, z1width);
      z2 = z - z1;
    }

    if (itest == 1) {
      G4cout << " z1 " << z1 << G4endl;
      G4cout << " z2 " << z2 << G4endl;
    }
    if (itest == 1) {
      G4cout << " zsymm " << zsymm << G4endl;
      G4cout << " nsymm " << nsymm << G4endl;
      G4cout << " asymm " << asymm << G4endl;
    }
    //    	CN =  UMASS(Zsymm , Nsymm + 1.E0) + UMASS(Zsymm, Nsymm - 1.E0)
    //    #            + 1.44E0 * (Zsymm)**2 /
    //    #            (r_null**2 * ((Asymm+1)**(1./3.) +
    //    #            (Asymm-1)**(1./3.))**2 )
    //    #            - 2.E0 * UMASS(Zsymm,Nsymm)
    //    #            - 1.44E0 * (Zsymm)**2 /
    //    #            (r_null * 2.E0 * (Asymm)**(1./3.))**2

    n1ucd = z1 * n/z;
    n2ucd = z2 * n/z;
    re1 = umass(z1,n1ucd,0.6) + umass(z2,n2ucd,0.6) +
      ecoul(z1,n1ucd,0.6,z2,n2ucd,0.6,2.0);
    re2 = umass(z1,n1ucd+1.,0.6) + umass(z2,n2ucd-1.,0.6) +
      ecoul(z1,n1ucd+1.,0.6,z2,n2ucd-1.,0.6,2.0);
    re3 = umass(z1,n1ucd+2.,0.6) + umass(z2,n2ucd-2.,0.6) +
      ecoul(z1,n1ucd+2.,0.6,z2,n2ucd-2.,0.6,2.0);
    reps2 = (re1-2.0*re2+re3)/2.0;
    reps1 = re2 - re1 -reps2;
    rn1_pol = -reps1/(2.0*reps2);
    n1mean = n1ucd + rn1_pol;
    n2mean = n - n1mean;

    if (itest == 1) {
      G4cout << " n1mean " << n1mean << G4endl;
      G4cout << " n2mean " << n2mean << G4endl;
    }

    cn = (umass(z1,n1mean+1.,0.0) + umass(z1,n1mean-1.,0.0) +
	  + umass(z2,n2mean+1.,0.0) + umass(z2,n2mean-1.,0.0)
	  - 2.0 * umass(z1,n1mean,0.0) +
	  - 2.0 * umass(z2,n2mean,0.0) ) * 0.5;
    //      This is an approximation! Coulomb energy is neglected.

    n1width = std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*e) / cn) );

    if (itest == 1) {
      G4cout << " cn " << cn << G4endl;
      G4cout << " n1width " << n1width << G4endl;
    }
	
    // random decision: N1R and N2R at scission, before evaporation: */
    //       N1R = dfloat(NINT(GAUSSHAZ(KKK,sngl(N1mean),sngl(N1width))));
    //      n1r = (float)( (int)(rnd.gaus(n1mean,n1width)) );
    n1r = (float)( (int)(gausshaz(k, n1mean,n1width)) );
    n2r = n - n1r;
    // Mass of first and second fragment */
    a1 = z1 + n1r;
    a2 = z2 + n2r;

    e1 = e*a1/(a1+a2);
    e2 = e - e*a1/(a1+a2);
    if (itest == 1) {
      G4cout << " n1r " << n1r << G4endl;
      G4cout << " n2r " << n2r << G4endl;
    }

  }

  if (itest == 1) {
    G4cout << " a1 " << a1 << G4endl;
    G4cout << " z1 " << z1 << G4endl;
    G4cout << " a2 " << a2 << G4endl;
    G4cout << " z2 " << z2 << G4endl;
    G4cout << " e1 " << e1 << G4endl;
    G4cout << " e2 " << e << G4endl;
  }

  //       /* Pre-neutron-emission total kinetic energy: */
  tker = (z1 * z2 * 1.44) /
    ( r0 * std::pow(a1,0.33333) * (1.0 + 2.0/3.0 * beta1) +
      r0 * std::pow(a2,0.33333) * (1.0 + 2.0/3.0 * beta2) + 2.0 );
  //       /* Pre-neutron-emission kinetic energy of 1. fragment: */
  ekin1 = tker * a2 / a;
  //       /* Pre-neutron-emission kinetic energy of 2. fragment: */
  ekin2 = tker * a1 / a;

  v1 = std::sqrt( (ekin1/a1) ) * 1.3887;
  v2 = std::sqrt( (ekin2/a2) ) * 1.3887;

  if (itest == 1) {
    G4cout << " kinetic energies " << G4endl;
    G4cout << " ekin1 " << ekin1 << G4endl;
    G4cout << " ekin2 " << ekin2 << G4endl;
  }
}

//       SUBROUTINE TRANSLAB(GAMREM,ETREM,CSREM,NOPART,NDEC)
void G4Abla::translab(G4double gamrem, G4double etrem, G4double csrem[4], G4int nopart, G4int ndec)
{
  // c Ce subroutine transforme dans un repere 1 les impulsions pcv des 
  // c particules acv, zcv et de cosinus directeurs xcv, ycv, zcv calculees 
  // c dans un repere 2.    
  // c La transformation de lorentz est definie par GAMREM (gamma) et
  // c ETREM (eta). La direction  du repere 2 dans 1 est donnees par les 
  // c cosinus directeurs ALREM,BEREM,GAREM (axe oz du repere 2).
  // c L'axe oy(2) est fixe par le produit vectoriel oz(1)*oz(2).
  // c Le calcul est fait pour les particules de NDEC a iv du common volant.
  // C Resultats dans le NTUPLE (common VAR_NTP) decale de NOPART (cascade).
    
  //       REAL*8  GAMREM,ETREM,ER,PLABI(3),PLABF(3),R(3,3)
  //       real*8  MASSE,PTRAV2,CSREM(3),UMA,MELEC,EL
  //       real*4 acv,zpcv,pcv,xcv,ycv,zcv
  //       common/volant/acv(200),zpcv(200),pcv(200),xcv(200),
  //      s              ycv(200),zcv(200),iv
      
  // 	parameter (max=250)                                                                       
  // 	real*4 EXINI,ENERJ,BIMPACT,PLAB,TETLAB,PHILAB,ESTFIS
  // 	integer AVV,ZVV,JREMN,KFIS,IZFIS,IAFIS
  //         common/VAR_NTP/MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP,
  //      +MULNTOT,BIMPACT,JREMN,KFIS,ESTFIS,IZFIS,IAFIS,NTRACK,
  //      +ITYPCASC(max),AVV(max),ZVV(max),ENERJ(max),PLAB(max),
  //      +TETLAB(max),PHILAB(max)
      
  //       DATA UMA,MELEC/931.4942,0.511/

  // C Matrice de rotation dans le labo:
  G4double sitet = std::sqrt(std::pow(csrem[1],2)+std::pow(csrem[2],2));
  G4double cstet = 0.0, siphi = 0.0, csphi = 0.0;
  G4double R[4][4];
  for(G4int init_i = 0; init_i < 4; init_i++) {
    for(G4int init_j = 0; init_j < 4; init_j++) {
      R[init_i][init_j] = 0.0;
    }
  }

  if(sitet > 1.0e-6) { //then
    cstet = csrem[3];
    siphi = csrem[2]/sitet;
    csphi = csrem[1]/sitet;	

    R[1][1] = cstet*csphi;
    R[1][2] = -siphi;
    R[1][3] = sitet*csphi;
    R[2][1] = cstet*siphi;
    R[2][2] = csphi;
    R[2][3] = sitet*siphi;
    R[3][1] = -sitet;
    R[3][2] = 0.0;
    R[3][3] = cstet;
  }
  else {
    R[1][1] = 1.0;
    R[1][2] = 0.0;
    R[1][3] = 0.0;
    R[2][1] = 0.0;
    R[2][2] = 1.0;
    R[2][3] = 0.0;
    R[3][1] = 0.0;
    R[3][2] = 0.0;
    R[3][3] = 1.0;
  } //endif

  G4int intp = 0;
  G4double el = 0.0;
  G4double masse = 0.0;
  G4double er = 0.0;
  G4double plabi[4];
  G4double ptrav2 = 0.0;
  G4double plabf[4];
  G4double bidon = 0.0;
  for(G4int init_i = 0; init_i < 4; init_i++) {
    plabi[init_i] = 0.0;
    plabf[init_i] = 0.0;
  }

  for(G4int i = ndec; i <= volant->iv; i++) { //do i=ndec,iv
    intp = i + nopart;
    varntp->ntrack = varntp->ntrack + 1;
    if(nint(volant->acv[i]) == 0 && nint(volant->zpcv[i]) == 0) {
      if(verboseLevel > 2) {
	G4cout <<"Error: Particles with A = 0 Z = 0 detected! " << G4endl;
      }
      continue;
    }
    if(varntp->ntrack >= VARNTPSIZE) {
      if(verboseLevel > 2) {
	G4cout <<"Error! Output data structure not big enough!" << G4endl;
      }
    }
    varntp->avv[intp] = nint(volant->acv[i]);
    varntp->zvv[intp] = nint(volant->zpcv[i]);
    varntp->itypcasc[intp] = 0;    
    // transformation de lorentz remnan --> labo:
    if (varntp->avv[intp] == -1) { //then
      masse = 138.00;	// cugnon
      // c		if (avv(intp).eq.1)  masse=938.2796	!cugnon
      // c		if (avv(intp).eq.4)  masse=3727.42	!ok
    }
    else {
      mglms(double(volant->acv[i]),double(volant->zpcv[i]),0, &el);
      // assert(isnan(el) == false);
      masse = volant->zpcv[i]*938.27 + (volant->acv[i] - volant->zpcv[i])*939.56 + el;
    } //end if
	
    er = std::sqrt(std::pow(volant->pcv[i],2) + std::pow(masse,2));
    // assert(isnan(er) == false);
    plabi[1] = volant->pcv[i]*(volant->xcv[i]);
    plabi[2] = volant->pcv[i]*(volant->ycv[i]);
    plabi[3] = er*etrem + gamrem*(volant->pcv[i])*(volant->zcv[i]);
	
    ptrav2 = std::pow(plabi[1],2) + std::pow(plabi[2],2) + std::pow(plabi[3],2);
    // assert(isnan(ptrav2) == false);
    varntp->plab[intp] = std::sqrt(ptrav2); 
    varntp->enerj[intp] = std::sqrt(ptrav2 + std::pow(masse,2)) - masse;

    // Rotation dans le labo:
    for(G4int j = 1; j <= 3; j++) { //do j=1,3
      plabf[j] = 0.0;
      for(G4int k = 1; k <= 3; k++) { //do k=1,3
	plabf[j] = plabf[j] + R[k][j]*plabi[k]; // :::Fixme::: (indices?)
      } // end do
    }  // end do
    // C impulsions dans le nouveau systeme copiees dans /volant/
    volant->pcv[i] = varntp->plab[intp];
    ptrav2 = std::sqrt(std::pow(plabf[1],2) + std::pow(plabf[2],2) + std::pow(plabf[3],2));
    if(ptrav2 >= 1.0e-6) { //then
      volant->xcv[i] = plabf[1]/ptrav2;
      volant->ycv[i] = plabf[2]/ptrav2;
      volant->zcv[i] = plabf[3]/ptrav2;
    }
    else {
      volant->xcv[i] = 1.0;
      volant->ycv[i] = 0.0;
      volant->zcv[i] = 0.0;
    } //endif
    // impulsions dans le nouveau systeme copiees dans /VAR_NTP/	
    if(varntp->plab[intp] >= 1.0e-6) { //then
      bidon = plabf[3]/(varntp->plab[intp]);
      // assert(isnan(bidon) == false);
      if(bidon > 1.0) { 
	bidon = 1.0;
      }
      if(bidon < -1.0) {
	bidon = -1.0;
      }
      varntp->tetlab[intp] = std::acos(bidon);
      sitet = std::sin(varntp->tetlab[intp]);
      varntp->philab[intp] = std::atan2(plabf[2],plabf[1]);        
      varntp->tetlab[intp] = varntp->tetlab[intp]*57.2957795;
      varntp->philab[intp] = varntp->philab[intp]*57.2957795;
    }
    else {
      varntp->tetlab[intp] = 90.0;
      varntp->philab[intp] = 0.0;
    } // endif
  } // end do
}
// C-------------------------------------------------------------------------

//       SUBROUTINE TRANSLABPF(MASSE1,T1,P1,CTET1,PHI1,GAMREM,ETREM,R,
//      s   PLAB1,GAM1,ETA1,CSDIR)
void G4Abla::translabpf(G4double masse1, G4double t1, G4double p1, G4double ctet1,
			G4double phi1, G4double gamrem, G4double etrem, G4double R[][4],
			G4double *plab1, G4double *gam1, G4double *eta1, G4double csdir[])
{
  // C Calcul de l'impulsion du PF (PLAB1, cos directeurs CSDIR(3)) dans le
  // C systeme remnant et des coefs de Lorentz GAM1,ETA1 de passage  
  // c du systeme PF --> systeme remnant.
  // c 
  // C Input: MASSE1, T1 (energie cinetique), CTET1,PHI1 (cosTHETA et PHI)
  // C                    (le PF dans le systeme du Noyau de Fission (NF)).
  // C	 GAMREM,ETREM les coefs de Lorentz systeme NF --> syst remnant, 
  // C        R(3,3) la matrice de rotation systeme NF--> systeme remnant.
  // C
  // C      
  //      	REAL*8 MASSE1,T1,P1,CTET1,PHI1,GAMREM,ETREM,R(3,3),
  //      s   PLAB1,GAM1,ETA1,CSDIR(3),ER,SITET,PLABI(3),PLABF(3)
     
  G4double er = t1 + masse1;
	
  G4double sitet = std::sqrt(1.0 - std::pow(ctet1,2));

  G4double plabi[4];
  G4double plabf[4];
  for(G4int init_i = 0; init_i < 4; init_i++) {
    plabi[init_i] = 0.0;
    plabf[init_i] = 0.0;
  }

  // C ----Transformation de Lorentz Noyau fissionnant --> Remnant:	
  plabi[1] = p1*sitet*std::cos(phi1);
  plabi[2] = p1*sitet*std::sin(phi1);
  plabi[3] = er*etrem + gamrem*p1*ctet1;
  
  // C ----Rotation du syst Noyaut Fissionant vers syst remnant:
  for(G4int j = 1; j <= 3; j++) { // do j=1,3
    plabf[j] = 0.0;
    for(G4int k = 1; k <= 3; k++) { //do k=1,3
      //      plabf[j] = plabf[j] + R[j][k]*plabi[k];
      plabf[j] = plabf[j] + R[k][j]*plabi[k];
    } //end do
  } //end do
  // C ----Cosinus directeurs et coefs de la transf de Lorentz dans le
  // c     nouveau systeme:	
  (*plab1) = std::pow(plabf[1],2) + std::pow(plabf[2],2) + std::pow(plabf[3],2);
  (*gam1) = std::sqrt(std::pow(masse1,2) + (*plab1))/masse1;
  (*plab1) = std::sqrt((*plab1));
  (*eta1) = (*plab1)/masse1;

  if((*plab1) <= 1.0e-6) { //then
    csdir[1] = 0.0;
    csdir[2] = 0.0;
    csdir[3] = 1.0;
  }
  else {   
    for(G4int i = 1; i <= 3; i++) { //do i=1,3
      csdir[i] = plabf[i]/(*plab1);
    } // end do
  } //endif
}

//       SUBROUTINE LOR_AB(GAM,ETA,Ein,Pin,Eout,Pout)
void G4Abla::lorab(G4double gam, G4double eta, G4double ein, G4double pin[],
		   G4double *eout, G4double pout[])
{
  // C  Transformation de lorentz brute pour vérifs.
  // C	P(3) = P_longitudinal (transformé)
  // C	P(1) et P(2) = P_transvers (non transformés)
  //       DIMENSION Pin(3),Pout(3)
  //       REAL*8 GAM,ETA,Ein

  pout[1] = pin[1];
  pout[2] = pin[2]; 
  (*eout) = gam*ein + eta*pin[3];
  pout[3] = eta*ein + gam*pin[3];
}

//       SUBROUTINE ROT_AB(R,Pin,Pout)
void G4Abla::rotab(G4double R[4][4], G4double pin[4], G4double pout[4])
{
  // C  Rotation d'un vecteur
  //       DIMENSION Pin(3),Pout(3)
  //       REAL*8 R(3,3)
      
  for(G4int i = 1; i <= 3; i++) { // do i=1,3
    pout[i] = 0.0;
    for(G4int j = 1; j <= 3; j++) { //do j=1,3
      //      pout[i] = pout[i] + R[i][j]*pin[j];
      pout[i] = pout[i] + R[j][i]*pin[j];
    } // enddo
  } //enddo
}

// Methods related to the internal ABLA random number generator. In
// the future the random number generation must be factored into its
// own class

void G4Abla::standardRandom(G4double *rndm, G4long *seed)
{
  (*seed) = (*seed); // Avoid warning during compilation.
  // Use Geant4 G4UniformRand
  (*rndm) = G4UniformRand();
}

G4double G4Abla::haz(G4int k)
{
  const G4int pSize = 110;
  static G4double p[pSize];
  static G4long ix = 0, i = 0;
  static G4double x = 0.0, y = 0.0, a = 0.0, haz = 0.0;
  //  k =< -1 on initialise                                        
  //  k = -1 c'est reproductible                                   
  //  k < -1 || k > -1 ce n'est pas reproductible

  // Zero is invalid random seed. Set proper value from our random seed collection:
  if(ix == 0) {
    ix = hazard->ial;
  }

  if (k <= -1) { //then                                             
    if(k == -1) { //then                                            
      ix = 0;
    }
    else {
      x = 0.0;
      y = secnds(int(x));
      ix = int(y * 100 + 43543000);
      if(mod(ix,2) == 0) {
	ix = ix + 1;
      }
    }

    // Here we are using random number generator copied from INCL code
    // instead of the CERNLIB one! This causes difficulties for
    // automatic testing since the random number generators, and thus
    // the behavior of the routines in C++ and FORTRAN versions is no
    // longer exactly the same!
    standardRandom(&x, &ix);
    for(G4int i = 0; i < pSize; i++) { //do i=1,110                                                 
      standardRandom(&(p[i]), &ix);
    }
    standardRandom(&a, &ix);
    k = 0;
  }

  i = nint(100*a)+1;
  haz = p[i];
  standardRandom(&a, &ix);
  p[i] = a;

  hazard->ial = ix;
  
  return haz;
}


G4double G4Abla::gausshaz(int k, double xmoy, double sig)
{
  // Gaussian random numbers:

  //   1005       C*** TIRAGE ALEATOIRE DANS UNE GAUSSIENNE DE LARGEUR SIG ET MOYENNE XMOY
  static G4int  iset = 0;
  static G4double v1,v2,r,fac,gset,gausshaz;

  if(iset == 0) { //then                                              
    do {
      v1 = 2.0*haz(k) - 1.0;
      v2 = 2.0*haz(k) - 1.0;
      r = std::pow(v1,2) + std::pow(v2,2);
    } while(r >= 1);

    fac = std::sqrt(-2.*std::log(r)/r);
    // assert(isnan(fac) == false);
    gset = v1*fac;
    gausshaz = v2*fac*sig+xmoy;
    iset = 1;
  }
  else {
    gausshaz=gset*sig+xmoy;
    iset=0;
  }

  return gausshaz;                                                         
}


// Utilities

G4double G4Abla::min(G4double a, G4double b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Abla::min(G4int a, G4int b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4double G4Abla::max(G4double a, G4double b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Abla::max(G4int a, G4int b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Abla::nint(G4double number)
{
  G4double intpart = 0.0;
  G4double fractpart = 0.0;
  fractpart = std::modf(number, &intpart);
  if(number == 0) {
    return 0;
  }
  if(number > 0) {
    if(fractpart < 0.5) {
      return int(std::floor(number));
    }
    else {
      return int(std::ceil(number));
    }
  }
  if(number < 0) {
    if(fractpart < -0.5) {
      return int(std::floor(number));
    }
    else {
      return int(std::ceil(number));
    }
  }

  return int(std::floor(number));
}

G4int G4Abla::secnds(G4int x)
{
  time_t mytime;
  tm *mylocaltime;

  time(&mytime);
  mylocaltime = localtime(&mytime);

  if(x == 0) {
    return(mylocaltime->tm_hour*60*60 + mylocaltime->tm_min*60 + mylocaltime->tm_sec);
  }
  else {
    return(mytime - x);
  }
}

G4int G4Abla::mod(G4int a, G4int b)
{
  if(b != 0) {
    return (a - (a/b)*b);
  }
  else {
    return 0;
  } 
}

G4double G4Abla::dmod(G4double a, G4double b)
{
  if(b != 0) {
    return (a - (a/b)*b);
  }
  else {
    return 0.0;
  } 
}

G4double G4Abla::dint(G4double a)
{
  G4double value = 0.0;

  if(a < 0.0) {
    value = double(std::ceil(a));
  }
  else {
    value = double(std::floor(a));
  }

  return value;
}

G4int G4Abla::idint(G4double a)
{
  G4int value = 0;

  if(a < 0) {
    value = int(std::ceil(a));
  }
  else {
    value = int(std::floor(a));
  }

  return value;
}

G4int G4Abla::idnint(G4double value)
{
  G4double valueCeil = int(std::ceil(value));
  G4double valueFloor = int(std::floor(value));

  if(std::fabs(value - valueCeil) < std::fabs(value - valueFloor)) {
    return int(valueCeil);
  }
  else {
    return int(valueFloor);
  }
}

G4double G4Abla::dmin1(G4double a, G4double b, G4double c)
{
  if(a < b && a < c) {
    return a;
  }
  if(b < a && b < c) {
    return b;
  }
  if(c < a && c < b) {
    return c;
  }
  return a;
}

G4double G4Abla::utilabs(G4double a)
{
  if(a > 0) {
    return a;
  }
  if(a < 0) {
    return (-1*a);
  }
  if(a == 0) {
    return a;
  }

  return a;
}
