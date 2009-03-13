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
// $Id: G4FissionBarrier.cc,v 1.6 2009-03-13 18:57:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// J.M Quesada 08/Mar/2009 New fission barrier has been added 
// (A.J.Sierk, PRC33, 1985,2039) 


#include "G4FissionBarrier.hh"

G4FissionBarrier::G4FissionBarrier(const G4FissionBarrier & ) : G4VFissionBarrier()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4FissionBarrier::copy_constructor meant to not be accessable.");
}


const G4FissionBarrier & G4FissionBarrier::operator=(const G4FissionBarrier & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4FissionBarrier::operator= meant to not be accessable.");
    return *this;
}

G4bool G4FissionBarrier::operator==(const G4FissionBarrier & ) const 
{
    return false;
}

G4bool G4FissionBarrier::operator!=(const G4FissionBarrier & ) const 
{
    return true;
}



G4double G4FissionBarrier::FissionBarrier(const G4int A, const G4int Z, const G4double U)
// JMQ 08/03/2009 New fission barrier
//Includes Sierk's parameterization as option (now default)   
{
  G4bool Sierk=true;
  if(Sierk){
    G4int  Lav=0;
    G4double B=SierkFissionBarrier(A,Z,Lav);   
    G4double BU=B/(1.0 + std::sqrt(U/(2.0*static_cast<G4double>(A))));
    //	G4cout<<"A = "<<A<<"  Z = "<<Z<<"  Barrier_Sierk = "<<B<<" .. B_Sierk(E*="<<U<<"MeV)="<<BU<<G4endl;
    //	G4cout<<"A = "<<A<<"  Z = "<<Z<<"  Barrier_Sierk(L="<<Lav<<")= "<<B<<G4endl;
    //  G4cout<<"   Sierk's fission barrier  (with L=0 & Ex dependence) is used"<<G4endl;
    //return B;
    return BU;
  }
  else{

    // Compute fission barrier according with Barashenkov's prescription for A >= 65 
    if (A >= 65){
      G4double B=BarashenkovFissionBarrier(A,Z);
      G4double BU=B/(1.0 + std::sqrt(U/(2.0*static_cast<G4double>(A))));
      //      G4cout<<"A = "<<A<<"  Z = "<<Z<<"  Barrier_Bara = "<<B<<" .. B_Bara(E*="<<U<<"MeV)="<<BU<<G4endl;
      //    G4cout<<"   Barashevkov's fission barrier is used"<<G4endl;
      return BU;
    }
    else return 100.0*GeV;
  }
}


G4double G4FissionBarrier::BarashenkovFissionBarrier(const G4int A, const G4int Z)
  // Calculates Fission Barrier heights 
{
  // Liquid drop model parameters for
  // surface energy of a spherical nucleus
  const G4double aSurf = 17.9439*MeV;
  // and coulomb energy
  const G4double aCoul = 0.7053*MeV;
  const G4int N = A - Z;
  const G4double k = 1.7826;

  // fissibility parameter
  G4double x = (aCoul/(2.0*aSurf))*(static_cast<G4double>(Z)*static_cast<G4double>(Z))/static_cast<G4double>(A);
  x /= (1.0 - k*(static_cast<G4double>(N-Z)/static_cast<G4double>(A))*
	(static_cast<G4double>(N-Z)/static_cast<G4double>(A)));
  
  // Liquid drop model part of Fission Barrier
  G4double BF0 = aSurf*std::pow(static_cast<G4double>(A),2./3.);
  if (x <= 2./3.) BF0 *= 0.38*(3./4.-x);
  else BF0 *= 0.83*(1. - x)*(1. - x)*(1. - x);


  // 
  G4double D = 1.248*MeV;
  D *= (static_cast<G4double>(N)-2.0*(N/2)) + (static_cast<G4double>(Z)-2.0*(Z/2));

  //    G4cout<<"x="<<x<<" BF0="<<BF0<<" D="<<D<<" deltaZ+deltaN ="<<SellPlusPairingCorrection(Z,N)<<G4endl;
  return BF0 + D - SellPlusPairingCorrection(Z,N);
}

G4double G4FissionBarrier::SierkFissionBarrier(const G4int Ia, const G4int Iz, G4int Il  )
  // JMQ 23/02/2009 Calculates fission barrier heights according to Sierk parameterization
{

  G4double a,a1,a2, aa,aj,ak, z, zz, amax, amin, amax2, amin2, bfis0, Bfis,Segs,Selmax;
  G4int i, j, m,k,l;
  G4double emncof[5][4],elmcof[5][4],emxcof[7][5],elzcof[7][7], pa[7], pz[7],pl[10];
  // G4double emncofi[4][5],elmcofi[4][5],emxcofi[5][7],elzcofi[7][7];
  G4double el,el20, el80, ell, elmax;
  G4double q, qa, qb, sel20, sel80;
  G4double x,y;
  G4double egs, egs1[5][7], egs2[5][7], egs3[5][7], egs4[5][7], egs5[5][7], egscof[5][7][5];
  // G4double egs1i[7][5], egs2i[7][5], egs3i[7][5], egs4i[7][5], egs5i[7][5];

  static const G4double 
    emncofi[4][5] = {{-9.01100e+2, -1.40818e+3, 2.77000e+3, -7.06695e+2, 8.89867e+2},
		     {1.35355e+4, -2.03847e+4, 1.09384e+4, -4.86297e+3,-6.18603e+2},
		     {-3.26367e+3, 1.62447e+3, 1.36856e+3, 1.31731e+3, 1.53372e+2},
		     {7.48863e+3, -1.21581e+4, 5.50281e+3, -1.33630e+3, 5.05367e-02}};

  static const G4double 
    elmcofi[4][5]={{1.84542e+3, -5.64002e+3, 5.66730e+3, -3.15150e+3, 9.54160e+2},
		   {-2.24577e+3, 8.56133e+3, -9.67348e+3, 5.81744e+3, -1.86997e+3},
		   {2.79772e+3, -8.73073e+3, 9.19706e+3, -4.91900e+3, 1.37283e+3}, 
		   {-3.01866e+1, 1.41161e+3, -2.85919e+3, 2.13016e+3, -6.49072e+2}};

  static const G4double 
    emxcofi[5][7] = {{-4.10652732e6, 1.00064947e7, -1.09533751e7, 
		      7.84797252e6, -3.78574926e6, 1.12237945e6, -1.77561170e5}, 
		     {1.08763330e7, -2.63758245e7, 2.85472400e7, 
		      -2.01107467e7, 9.48373641e6, -2.73438528e6, 4.13247256e5}, 
		     {-8.76530903e6, 2.14250513e7, -2.35799595e7, 
		      1.70161347e7, -8.23738190e6, 2.42447957e6, -3.65427239e5},
		     {6.30258954e6, -1.52999004e7, 1.65640200e7, 
		      -1.16695776e7, 5.47369153e6, -1.54986342e6, 2.15409246e5},
		     {-1.45539891e6, 3.64961835e6, -4.21267423e6, 
		      3.24312555e6, -1.67927904e6, 5.23795062e5, -7.66576599e4}};
  
  static const G4double 
    elzcofi[7][7] ={{5.11819909e+5, -1.30303186e+6, 1.90119870e+6, 
		     -1.20628242e+6, 5.68208488e+5, 5.48346483e+4, -2.45883052e+4}, 
		    {-1.13269453e+6, 2.97764590e+6, -4.54326326e+6, 
		     3.00464870e+6, -1.44989274e+6, -1.02026610e+5, 6.27959815e+4}, 
		    {1.37543304e+6, -3.65808988e+6, 5.47798999e+6, 
		     -3.78109283e+6, 1.84131765e+6, 1.53669695e+4, -6.96817834e+4}, 
		    {-8.56559835e+5, 2.48872266e+6, -4.07349128e+6, 
		     3.12835899e+6, -1.62394090e+6, 1.19797378e+5, 4.25737058e+4}, 
		    {3.28723311e+5, -1.09892175e+6, 2.03997269e+6, 
		     -1.77185718e+6, 9.96051545e+5, -1.53305699e+5, -1.12982954e+4}, 
		    {4.15850238e+4, 7.29653408e+4, -4.93776346e+5, 
		     6.01254680e+5, -4.01308292e+5, 9.65968391e+4, -3.49596027e+3}, 
		    {-1.82751044e+5, 3.91386300e+5, -3.03639248e+5, 
		     1.15782417e+5, -4.24399280e+3, -6.11477247e+3, 3.66982647e+2}};


  static const G4double 
    egs1i[7][5]={{- 1.781665232e6, -2.849020290e6, 9.546305856e5, 2.453904278e5, 3.656148926e5}, 
		 {4.358113622e6, 6.960182192e6, -2.381941132e6, -6.262569370e5, -9.026606463e5}, 
		 { -4.804291019e6, -7.666333374e6, 2.699742775e6, 7.415602390e5,  1.006008724e6}, 
		 {3.505397297e6, 5.586825123e6, -2.024820713e6, -5.818008462e5, -7.353683218e5}, 
		 {-1.740990985e6, -2.759325148e6, 1.036253535e6, 3.035749715e5, 3.606919356e5}, 
		 {5.492532874e5, 8.598827288e5, -3.399809581e5, -9.852362945e4, -1.108872347e5}, 
		 {-9.229576432e4, -1.431344258e5,5.896521547e4, 1.772385043e4, 1.845424227e4}};

  static const G4double 
    egs2i[7][5]={{4.679351387e6, 7.707630513e6, -2.718115276e6, -9.845252314e5, -1.107173456e6}, 
		 {-1.137635233e7, -1.870617878e7, 6.669154225e6, 2.413451470e6, 2.691480439e6}, 
		 {1.237627138e7, 2.030222826e7, -7.334289876e6, -2.656357635e6, -2.912593917e6}, 
		 {-8.854155353e6, -1.446966194e7, 5.295832834e6, 1.909275233e6, 2.048899787e6}, 
		 {4.290642787e6, 6.951223648e6, -2.601557110e6, -9.129731614e5, -9.627344865e5}, 
		 {-1.314924218e6, -2.095971932e6, 8.193066795e5, 2.716279969e5, 2.823297853e5}, 
		 {2.131536582e5, 3.342907992e5, -1.365390745e5, -4.417841315e4, -4.427025540e4}};

  static const G4double 
    egs3i[7][5]={{- 3.600471364e6, -5.805932202e6, 1.773029253e6, 4.064280430e5, 7.419581557e5}, 
		 {8.829126250e6, 1.422377198e7,-4.473342834e6, -1.073350611e6, -1.845960521e6}, 
		 {-9.781712604e6, -1.575666314e7, 5.161226883e6, 1.341287330e6, 2.083994843e6}, 
		 {7.182555931e6, 1.156915972e7, -3.941330542e6, -1.108259560e6, -1.543982755e6}, 
		 {-3.579820035e6, -5.740079339e6, 2.041827680e6, 5.981648181e5, 7.629263278e5}, 
		 {1.122573403e6, 1.777161418e6, -6.714631146e5, -1.952833263e5, -2.328129775e5}, 
		 {-1.839672155e5, -2.871137706e5, 1.153532734e5, 3.423868607e4, 3.738902942e4}};

  static const G4double 
    egs4i[7][5]={{2.421750735e6, 4.107929841e6, -1.302310290e6, -5.267906237e5, -6.197966854e5}, 
		 {-5.883394376e6,-9.964568970e6, 3.198405768e6, 1.293156541e6, 1.506909314e6}, 
		 {6.387411818e6, 1.079547152e7, -3.517981421e6, -1.424705631e6, -1.629099740e6}, 
		 {-4.550695232e6, -7.665548805e6, 2.530844204e6, 1.021187317e6, 1.141553709e6}, 
		 {2.182540324e6, 3.646532772e6, -1.228378318e6, -4.813626449e5, -5.299974544e5}, 
		 {-6.518758807e5, -1.070414288e6, 3.772592079e5, 1.372024952e5, 1.505359294e5}, 
		 {9.952777968e4, 1.594230613e5, -6.029082719e4, -2.023689807e4, -2.176008230e4}};

  static const G4double 
    egs5i[7][5]={{- 4.902668827e5, -8.089034293e5, 1.282510910e5, -1.704435174e4, 8.876109934e4}, 
		 {1.231673941e6, 2.035989814e6, -3.727491110e5, 4.071377327e3, -2.375344759e5}, 
		 {-1.429330809e6, -2.376692769e6, 5.216954243e5, 7.268703575e4, 3.008350125e5}, 
		 {1.114306796e6, 1.868800148e6, -4.718718351e5, -1.215904582e5, -2.510379590e5}, 
		 {-5.873353309e5, -9.903614817e5, 2.742543392e5, 9.055579135e4, 1.364869036e5}, 
		 {1.895325584e5, 3.184776808e5, -9.500485442e4, -3.406036086e4, -4.380685984e4}, 
		 {-2.969272274e4, -4.916872669e4,1.596305804e4, 5.741228836e3, 6.669912421e3}};
 
  for(i=0;i<4;i++) {
    for(j=0;j<5;j++) {
      emncof[j][i]=emncofi[i][j];
      elmcof[j][i]=elmcofi[i][j];
    }
  }    

  for(i=0;i<7;i++) {
    for(j=0;j<7;j++) {
      elzcof[j][i]=elzcofi[i][j];
    }
  }
  
  for(i=0;i<7;i++) {
    for(j=0;j<5;j++) {
      egs1[j][i]=egs1i[i][j];
      egs2[j][i]=egs2i[i][j];
      egs3[j][i]=egs3i[i][j];
      egs4[j][i]=egs4i[i][j];
      egs5[j][i]=egs5i[i][j];
    }
  }

  for(i=0;i<5;i++) {
    for(j=0;j<7;j++) {
      emxcof[j][i]=emxcofi[i][j];
      egscof[i][j][0]=egs1[i][j];
      egscof[i][j][1]=egs2[i][j];
      egscof[i][j][2]=egs3[i][j];
      egscof[i][j][3]=egs4[i][j];
      egscof[i][j][4]=egs5[i][j];
    }
  }
 
  if(Iz>=19)
    {
      if(Iz>111){
	// G4cout<<" *  *  *  *  BARFIT CALLED WITH  Z  LESS THAN 19 OR  GREATER THAN 111.  BFIS IS SET TO 0.0.  *  *  *  "<<G4endl;
      }
      else if(Iz>102 && Il>0){
	//       G4cout<<"*  *  *  *  BARFIT CALLED WITH  Z  GREATER THAN 102 AND  L  NOT EQUAL TO ZERO.  BFIS IS SET TO 0.0.  *  *  *  *"<<G4endl;
      }
      else {
	z = (G4double)Iz;
	a = (G4double)Ia;
	el= (G4double)Il;
	amin = 1.2*z + 0.01*z*z;
	amax = 5.8*z - 0.024*z*z;
	if (a<amin || a>amax){ 
	  //G4cout<<"*  *  *  *  BARFIT CALLED WITH  A = "<<Ia<<"OUTSIDE THE ALLOWED VALUES FOR Z " <<G4endl;
	  std::ostringstream errOs;
	  errOs << "BARFIT CALLED WITH  A = "<<Ia<<"OUTSIDE THE ALLOWED VALUES FOR Z"  <<G4endl;
	  throw G4HadronicException(__FILE__, __LINE__, errOs.str());
	  return 0.;
	}
	else  {
	  aa = a/400.;
	  zz = z/100.;
	  bfis0 = 0.;
	  LPOLY(zz, 7, pz);
	  //       cout << "pz(7) = " << pz[6] << endl;
	  LPOLY(aa, 7, pa);
	  for  (i=0;i< 7;i++)
	   for (j = 0;j< 7;j++){
	     bfis0 = bfis0 + elzcof[j][i]*pz[j]*pa[i];
	     //		cout<<" elzcof("<<j+1<<","<<i+1<<")="<< elzcof[j][i]<<"pz("<<j+1<<")="<<pz[j]<<" pa("<<i+1<<")="<<pa[i]<<" bfis0= "<<bfis0<<endl;
	   }
	  egs=0.;
	  Segs=egs;
	  Bfis = bfis0;
	  amin2 = 1.4*z + 0.009*z*z;
	  amax2 = 20. + 3.0*z;

	  if((a<amin2 - 5.e0 || a>amax2 + 10.e0) && Il>0){
	    std::ostringstream errOs;
	    errOs << "BARFIT CALLED WITH  A = "<<Ia<<"OUTSIDE THE ALLOWED VALUES FOR Z for NON ZERO Il="<<Il<<G4endl;
	    throw G4HadronicException(__FILE__, __LINE__, errOs.str());
	    return 0.;
	  }else 
	    {               
	      el80 = 0.;
	      el20 = 0.;
	      elmax = 0.;
	      for (i = 0; i<4;i++)		  
		for (j = 0;j<5;j++){
		  el80 = el80 + (elmcof[j][i])*pz[j]*pa[i];
		  el20 = el20 + (emncof[j][i])*pz[j]*pa[i];
		}
		                      
	      sel80 = el80;
	      sel20 = el20;
              
	      for (i = 0;i< 5;i++)
		for (j = 0;j< 7;j++)
		  elmax = elmax + emxcof[j][i]*pz[j]*pa[i];
                     
	      Selmax = elmax;
	      if (Il<1) return (Bfis);
	      x = sel20/Selmax;
	      y = sel80/Selmax;
	      if(el<=sel20){
		q = 0.2/(sel20*sel20*sel80*sel80*(sel20 - sel80));
		qa = q*(4.*sel80*sel80*sel80 - sel20*sel20*sel20);
		qb = -q*(4.*sel80*sel80 - sel20*sel20);
		Bfis = Bfis*(1. + qa*el*el + qb*el*el*el);
	      } else {
		aj = (( - 20.*x*x*x*x*x) + 25.*x*x*x*x - 4.)*(y - 1.)*(y-1)*y*y;
                ak = (( - 20.*y*y*y*y*y) + 25.*y*y*y*y - 1.)*(x - 1.)*(x-1)*x*x;
		q = 0.2/((y - x)*((1.-x)*(1.-y)*x*y)*((1.-x)*(1.-y)*x*y));
		qa = q*(aj*y - ak*x);
		qb = -q*(aj*(2.*y + 1.) - ak*(2.*x + 1.));
		z = el/Selmax;
		a1 = 4.*z*z*z*z*z - 5.*z*z*z*z + 1.;
		a2 = qa*(2.*z + 1.);
		Bfis = Bfis*(a1 + (z - 1.)*(a2 + qb*z)*z*z*(z - 1.));
	      }
	      if(Bfis<=0.0e0 ||  el>Selmax)Bfis = 0.0;
	      if(el>Selmax && Il!=1000) return Bfis;
	      //
	      //                 NOW CALCULATE ROTATING GROUND-STATE ENERGY
	      //
	      ell = el/elmax;
	      if(Il==1000)ell = 1.e0;
	      LPOLY(ell, 9, pl);
	      for (k = 0;k< 5;k++) {
		for (l = 0;l< 7;l++) {
		  for (m = 0; m<5;m++) {
		    egs = egs + egscof[m][l][k]*pz[l]*pa[k]*pl[2*m];
		  }
		}
	      }
	      Segs = egs;
	      if(Segs<0.0)Segs = 0.0;
	      //	 G4cout<<"A = "<<Ia<<"  Z = "<<Iz<<"  Bar_fis = "<<Bfis<<G4endl;
	      return Bfis;
	    }
	}
      }
 
      Bfis = 0.0;
      Segs = 0.0;
      Selmax = 0.0;
      //  G4cout<<"A = "<<Ia<<"  Z = "<<Iz<<"  Bar_fis = "<<Bfis<<G4endl;
      return Bfis;
    }
  // G4cout<<"*  *  *  *  BARFIT CALLED WITH  Z  LESS THAN 19 BFIS IS SET TO 100.0.  *  *  *  * "<<G4endl;
  Bfis = 100.;
  Segs =0.0;
  Selmax=0.0;

  //G4cout<<"A = "<<Ia<<"  Z = "<<Iz<<"  Bar_fis = "<<Bfis<<G4endl;
  return Bfis;
}
    
//JMQ 20/02/2009 new method for legendre polynomials (Sierk's barrier) 
 void G4FissionBarrier::LPOLY( G4double X, G4int  N, G4double* Pl)
{
  Pl[0] = 1.0;
  Pl[1] = X;
 
  for( G4int i=2; i< N;i++)
    Pl[i] = ((2*(i+1) - 3)*X*Pl[i-1] - ((i+1)-2)*Pl[i-2])/i;
}
