//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4InitXscPAI.cc,v 1.2 2004-04-15 09:13:16 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4InitXscPAI.cc -- class implementation file
//
// GEANT 4 class implementation file
//
// For information related to this code, please, contact
// the Geant4 Collaboration.
//
// R&D: Vladimir.Grichine@cern.ch
//
// History:
//



#include "G4InitXscPAI.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4Poisson.hh"
#include "G4Integrator.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4SandiaTable.hh"



// Local class constants

const G4double G4InitXscPAI::fDelta = 0.005 ; // energy shift from interval border

//////////////////////////////////////////////////////////////////
//
// Constructor
//

G4InitXscPAI::G4InitXscPAI(G4MaterialCutsCouple* matCC)
{
  G4int i, j, matIndex;
 
  fDensity         = matCC->GetMaterial()->GetDensity();
  fElectronDensity = matCC->GetMaterial()->GetElectronDensity();
  matIndex         = matCC->GetMaterial()->GetIndex();

  fSandia          = new G4SandiaTable(matIndex);
  fIntervalNumber  = fSandia->GetMaxInterval()-1;

  fMatSandiaMatrix = new G4OrderedTable();
 
  for (i = 0; i < fIntervalNumber; i++)
  {
    fMatSandiaMatrix->push_back(new G4DataVector(5,0.));
  }	         	
  for (G4int i = 0; i < fIntervalNumber; i++)
  {
    (*(*fMatSandiaMatrix)[i])[0] = fSandia->GetSandiaMatTable(i,0);

    for(j = 1; j < 5 ; j++)
    {
      (*(*fMatSandiaMatrix)[i])[j] = fSandia->GetSandiaMatTable(i,j)*fDensity;
    }     
  }
  KillCloseIntervals();
  Normalisation();

}




////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4InitXscPAI::~G4InitXscPAI()
{
  
}

////////////////////////////////////////////////////////////////////////
//
// Kill close intervals, recalculate fIntervalNumber

void G4InitXscPAI::KillCloseIntervals()
{
  G4int i, j, k;
  G4double energy1, energy2; 
	         	
  for( i = 0 ; i < fIntervalNumber - 1 ; i++ )
  {
    energy1 = (*(*fMatSandiaMatrix)[i])[0];
    energy2 = (*(*fMatSandiaMatrix)[i+1])[0];

    if( energy2 - energy1 > 1.5*fDelta*(energy1 + energy2) )  continue ;
    else
    {
      for(j = i; j < fIntervalNumber-1; j++)
      {
        for( k = 0; k < 5; k++ )
        {
          (*(*fMatSandiaMatrix)[j])[k] = (*(*fMatSandiaMatrix)[j+1])[k];
	}
      }
      fIntervalNumber-- ;
      i-- ;
    }
  }
  
}

////////////////////////////////////////////////////////////////////////
//
// Kill close intervals, recalculate fIntervalNumber

void G4InitXscPAI::Normalisation()
{
  G4int i, j;
  G4double energy1, energy2, delta, cof; // , shift;

  energy1 = (*(*fMatSandiaMatrix)[fIntervalNumber-1])[0];
  energy2 = 2.*(*(*fMatSandiaMatrix)[fIntervalNumber-1])[0];

 
  cof = RutherfordIntegral(fIntervalNumber-1,energy1,energy2);
	         	
  for( i = fIntervalNumber-2; i >= 0; i-- )
  {
    energy1 = (*(*fMatSandiaMatrix)[i])[0];
    energy2 = (*(*fMatSandiaMatrix)[i+1])[0];

    cof += RutherfordIntegral(i,energy1,energy2);
    // G4cout<<"norm. cof = "<<cof<<G4endl;
  }
  fNormalizationCof  = 2*pi*pi*hbarc*hbarc*fine_structure_const/electron_mass_c2 ;
  fNormalizationCof *= fElectronDensity;
  delta = fNormalizationCof - cof;
  fNormalizationCof /= cof;
  //  G4cout<<"G4InitXscPAI::fNormalizationCof/cof = "<<fNormalizationCof
  //    <<";  at delta ="<<delta<<G4endl ;

  for (G4int i = 0; i < fIntervalNumber; i++) // renormalisation on QM sum rule
  {
    for(j = 1; j < 5 ; j++)
    {
      (*(*fMatSandiaMatrix)[i])[j] *= fNormalizationCof;
    }     
  }
  /* 
  if(delta > 0) // shift the first energy interval
  {
    for(i=1;i<100;i++)
    {
      energy1 = (1.-i/100.)*(*(*fMatSandiaMatrix)[0])[0];
      energy2 = (*(*fMatSandiaMatrix)[0])[0];
      shift   = RutherfordIntegral(0,energy1,energy2);
      G4cout<<shift<<"\t";
      if(shift >= delta) break;
    }
    (*(*fMatSandiaMatrix)[0])[0] = energy1;
    cof += shift;
  }
  else if(delta < 0)
  {
    for(i=1;i<100;i++)
    {
      energy1 = (*(*fMatSandiaMatrix)[0])[0];
      energy2 = (*(*fMatSandiaMatrix)[0])[0] + 
              ( (*(*fMatSandiaMatrix)[0])[0] - (*(*fMatSandiaMatrix)[0])[0] )*i/100.;
      shift   = RutherfordIntegral(0,energy1,energy2);
      if( shift >= abs(delta) ) break;
    }
    (*(*fMatSandiaMatrix)[0])[0] = energy2;
    cof -= shift;
  }
  G4cout<<G4cout<<"G4InitXscPAI::fNormalizationCof/cof = "<<fNormalizationCof/cof
        <<";  at delta ="<<delta<<"  and i = "<<i<<G4endl ;
 */ 
}





////////////////////////////////////////////////////////////////////
//
// Integration over electrons that could be considered
// quasi-free at energy transfer of interest

G4double G4InitXscPAI::RutherfordIntegral( G4int k,
				            G4double x1,
			  	            G4double x2   )
{
   G4double  c1, c2, c3, a1, a2, a3, a4 ;

   a1 = (*(*fMatSandiaMatrix)[k])[1]; 
   a2 = (*(*fMatSandiaMatrix)[k])[2]; 
   a3 = (*(*fMatSandiaMatrix)[k])[3]; 
   a4 = (*(*fMatSandiaMatrix)[k])[4]; 
   // G4cout<<"RI: x1 = "<<x1<<"; "<<"x2 = "<<x2<<G4endl;   
   c1 = (x2 - x1)/x1/x2 ;
   c2 = (x2 - x1)*(x2 + x1)/x1/x1/x2/x2 ;
   c3 = (x2 - x1)*(x1*x1 + x1*x2 + x2*x2)/x1/x1/x1/x2/x2/x2 ;
   // G4cout<<" RI: c1 = "<<c1<<"; "<<"c2 = "<<c2<<"; "<<"c3 = "<<c3<<G4endl;   
   
   return  a1*log(x2/x1) + a2*c1 + a3*c2/2 + a4*c3/3 ;

}   // end of RutherfordIntegral 

///////////////////////////////////////////////////////////////
//
//  Integrate photo-absorption cross-section from I1 up to omega

G4double G4InitXscPAI::IntegralTerm(G4double omega)
{
  G4int i;
  G4double energy1, energy2, result = 0.; 
	         	
  for( i = 0; i <= fIntervalTmax; i++ )
  {
    if(i == fIntervalTmax) 
    {
      energy1 = (*(*fMatSandiaMatrix)[i])[0];
      result += RutherfordIntegral(i,energy1,omega);
    }
    else 
    {
      if( omega <= (*(*fMatSandiaMatrix)[i+1])[0])
      {
        energy1 = (*(*fMatSandiaMatrix)[i])[0];
        result += RutherfordIntegral(i,energy1,omega);
        break;
      }
      else
      {
        energy1 = (*(*fMatSandiaMatrix)[i])[0];
        energy2 = (*(*fMatSandiaMatrix)[i+1])[0];
        result += RutherfordIntegral(i,energy1,energy2);
      }
    }
    // G4cout<<"IntegralTerm<<"("<<omega<<")"<<" = "<<result<<G4endl;
  }
  return result;
}


////////////////////////////////////////////////////////////////
//
// Imaginary part of dielectric constant
// (G4int k - interval number, G4double en1 - energy point)

G4double G4InitXscPAI::ImPartDielectricConst( G4int    k ,
			                       G4double energy1 )
{
   G4double energy2,energy3,energy4,a1,a2,a3,a4,result;

   a1 = (*(*fMatSandiaMatrix)[k])[1]; 
   a2 = (*(*fMatSandiaMatrix)[k])[2]; 
   a3 = (*(*fMatSandiaMatrix)[k])[3]; 
   a4 = (*(*fMatSandiaMatrix)[k])[4]; 

   energy2 = energy1*energy1;
   energy3 = energy2*energy1;
   energy4 = energy3*energy1;
   
   result = a1/energy1+a2/energy2+a3/energy3+a4/energy4 ;  
   result *=hbarc/energy1 ;
   
   return result ;

}  // end of ImPartDielectricConst 


//////////////////////////////////////////////////////////////////////////////
//
// Real part of dielectric constant minus unit: epsilon_1 - 1
// (G4double enb - energy point)
//

G4double G4InitXscPAI::RePartDielectricConst(G4double enb)
{
  G4int i;       
   G4double x0, x02, x03, x04, x05, x1, x2, a1,a2,a3,a4,xx1 ,xx2 , xx12,
            c1, c2, c3, cof1, cof2, xln1, xln2, xln3, result ;

   x0 = enb ;
   result = 0 ;
   
   for( i = 0; i < fIntervalNumber-1; i++)
   {
      x1 = (*(*fMatSandiaMatrix)[i])[0];
      x2 = (*(*fMatSandiaMatrix)[i+1])[0] ;

      a1 = (*(*fMatSandiaMatrix)[i])[1]; 
      a2 = (*(*fMatSandiaMatrix)[i])[2]; 
      a3 = (*(*fMatSandiaMatrix)[i])[3]; 
      a4 = (*(*fMatSandiaMatrix)[i])[4];
 
      if( abs(x0-x1) < 0.5*(x0+x1)*fDelta ) 
      {
        if(x0 >= x1) x0 = x1*(1+fDelta);
        else         x0 = x1*(1-fDelta);
      } 
      if( abs(x0-x2) < 0.5*(x0+x2)*fDelta ) 
      {
        if(x0 >= x2) x0 = x2*(1+fDelta);
        else         x0 = x2*(1-fDelta);
      }
      xx1 = x1 - x0 ;
      xx2 = x2 - x0 ;
      xx12 = xx2/xx1 ;
      
      if( xx12 < 0 ) xx12 = -xx12;
      
      xln1 = log(x2/x1) ;
      xln2 = log(xx12) ;
      xln3 = log((x2 + x0)/(x1 + x0)) ;

      x02 = x0*x0 ;
      x03 = x02*x0 ;
      x04 = x03*x0 ;
      x05 = x04*x0;

      c1  = (x2 - x1)/x1/x2 ;
      c2  = (x2 - x1)*(x2 +x1)/x1/x1/x2/x2 ;
      c3  = (x2 -x1)*(x1*x1 + x1*x2 + x2*x2)/x1/x1/x1/x2/x2/x2 ;

      result -= (a1/x02 + a3/x04)*xln1 ;
      result -= (a2/x02 + a4/x04)*c1 ;
      result -= a3*c2/2/x02 ;
      result -= a4*c3/3/x02 ;

      cof1 = a1/x02 + a3/x04 ;
      cof2 = a2/x03 + a4/x05 ;

      result += 0.5*(cof1 +cof2)*xln2 ;
      result += 0.5*(cof1 - cof2)*xln3 ;
   } 
   result *= 2*hbarc/pi ;
   
   return result ;

}   // end of RePartDielectricConst 

//////////////////////////////////////////////////////////////////////
//
// PAI differential cross-section in terms of
// simplified Allison's equation
//

G4double G4InitXscPAI::DifPAIxSection( G4double omega )
{
  G4int i = fCurrentInterval;
  G4double  betaGammaSq = fBetaGammaSq;       
  G4double integralTerm = IntegralTerm(omega);
  G4double be2,cof,x1,x2,x3,x4,x5,x6,x7,x8,result ;
  G4double epsilonRe = RePartDielectricConst(omega);
  G4double epsilonIm = ImPartDielectricConst(i,omega);
  G4double be4 ;
  G4double betaBohr2 = fine_structure_const*fine_structure_const ;
  G4double betaBohr4 = betaBohr2*betaBohr2*4.0 ;
  be2 = betaGammaSq/(1 + betaGammaSq) ;
  be4 = be2*be2 ;
 
   cof = 1 ;
   x1 = log(2*electron_mass_c2/omega) ;

   if( betaGammaSq < 0.01 ) x2 = log(be2) ;
   else
   {
     x2 = -log( (1/betaGammaSq - epsilonRe)*
	        (1/betaGammaSq - epsilonRe) + 
	        epsilonIm*epsilonIm )/2 ;
   }
   if( epsilonIm == 0.0 || betaGammaSq < 0.01 )
   {
     x6=0 ;
   }
   else
   {
     x3 = -epsilonRe + 1/betaGammaSq ;
     x5 = -1 - epsilonRe + be2*((1 +epsilonRe)*(1 + epsilonRe) +
	   epsilonIm*epsilonIm) ;

     x7 = atan2(epsilonIm,x3) ;
     x6 = x5 * x7 ;
   }
    // if(fImPartDielectricConst[i] == 0) x6 = 0 ;
   
   x4 = ((x1 + x2)*epsilonIm + x6)/hbarc ;
   //   if( x4 < 0.0 ) x4 = 0.0 ;
   x8 = (1 + epsilonRe)*(1 + epsilonRe) + 
        epsilonIm*epsilonIm;

   result = (x4 + cof*integralTerm/omega/omega) ;
   if(result < 1.0e-8) result = 1.0e-8 ;
   result *= fine_structure_const/be2/pi ;
   //   result *= (1-exp(-beta/betaBohr))*(1-exp(-beta/betaBohr)) ;
   //  result *= (1-exp(-be2/betaBohr2)) ;
   result *= (1-exp(-be4/betaBohr4)) ;
   if(fDensity >= 0.1)
   { 
      result /= x8 ;
   }
   return result ;

} // end of DifPAIxSection 

//////////////////////////////////////////////////////////////////////
//
// Differential PAI dEdx(omega)=omega*dNdx(omega)
//

G4double G4InitXscPAI::DifPAIdEdx( G4double omega )
{
  G4double dEdx = omega*DifPAIxSection(omega);
  return dEdx;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculation od dN/dx of collisions with creation of Cerenkov pseudo-photons

G4double G4InitXscPAI::PAIdNdxCerenkov( G4double omega  )
{        
  G4int i = fCurrentInterval;
  G4double  betaGammaSq = fBetaGammaSq;       
  G4double epsilonRe = RePartDielectricConst(omega);
  G4double epsilonIm = ImPartDielectricConst(i,omega);

  G4double cof, logarithm, x3, x5, argument, modul2, dNdxC ; 
  G4double be2, be4, betaBohr2,betaBohr4,cofBetaBohr ;

   cof         = 1.0 ;
   cofBetaBohr = 4.0 ;
   betaBohr2   = fine_structure_const*fine_structure_const ;
   betaBohr4   = betaBohr2*betaBohr2*cofBetaBohr ;

   be2 = betaGammaSq/(1 + betaGammaSq) ;
   be4 = be2*be2 ;

   if( betaGammaSq < 0.01 ) logarithm = log(1.0+betaGammaSq) ; // 0.0 ;
   else
   {
     logarithm  = -log( (1/betaGammaSq - epsilonRe)*
	                (1/betaGammaSq - epsilonRe) + 
	                epsilonIm*epsilonIm )*0.5 ;
     logarithm += log(1+1.0/betaGammaSq) ;
   }

   if( epsilonIm == 0.0 || betaGammaSq < 0.01 )
   {
     argument = 0.0 ;
   }
   else
   {
     x3 = -epsilonRe + 1.0/betaGammaSq ;
     x5 = -1.0 - epsilonRe +
          be2*((1.0 +epsilonRe)*(1.0 + epsilonRe) +
	  epsilonIm*epsilonIm) ;
     if( x3 == 0.0 ) argument = 0.5*pi;
     else            argument = atan2(epsilonIm,x3) ;
     argument *= x5  ;
   }   
   dNdxC = ( logarithm*epsilonIm + argument )/hbarc ;
  
   if(dNdxC < 1.0e-8) dNdxC = 1.0e-8 ;

   dNdxC *= fine_structure_const/be2/pi ;

   dNdxC *= (1-exp(-be4/betaBohr4)) ;

   if(fDensity >= 0.1)
   { 
      modul2 = (1.0 + epsilonRe)*(1.0 + epsilonRe) + 
                    epsilonIm*epsilonIm;
      dNdxC /= modul2 ;
   }
   return dNdxC ;

} // end of PAIdNdxCerenkov 

//////////////////////////////////////////////////////////////////////////
//
// Calculation od dN/dx of collisions with creation of longitudinal EM
// excitations (plasmons, delta-electrons)

G4double G4InitXscPAI::PAIdNdxPlasmon( G4double omega )
{        
  G4int i = fCurrentInterval;
  G4double  betaGammaSq = fBetaGammaSq;       
  G4double integralTerm = IntegralTerm(omega);
  G4double epsilonRe = RePartDielectricConst(omega);
  G4double epsilonIm = ImPartDielectricConst(i,omega);

   G4double cof, resonance, modul2, dNdxP ;
   G4double be2, be4, betaBohr2, betaBohr4, cofBetaBohr ;

   cof = 1 ;
   cofBetaBohr = 4.0 ;
   betaBohr2   = fine_structure_const*fine_structure_const ;
   betaBohr4   = betaBohr2*betaBohr2*cofBetaBohr ;

   be2 = betaGammaSq/(1 + betaGammaSq) ;
   be4 = be2*be2 ;
 
   resonance  = log(2*electron_mass_c2*be2/omega) ;  
   resonance *= epsilonIm/hbarc ;


   dNdxP = ( resonance + cof*integralTerm/omega/omega ) ;

   if( dNdxP < 1.0e-8 ) dNdxP = 1.0e-8 ;

   dNdxP *= fine_structure_const/be2/pi ;
   dNdxP *= (1-exp(-be4/betaBohr4)) ;

   if( fDensity >= 0.1 )
   { 
     modul2 = (1 + epsilonRe)*(1 + epsilonRe) + 
               epsilonIm*epsilonIm;
     dNdxP /= modul2 ;
   }
   return dNdxP ;

} // end of PAIdNdxPlasmon 

////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI integral cross-section
// = specific primary ionisation, 1/cm
// 

G4double G4InitXscPAI::IntegralPAIxSection(G4double omega, G4double bg2, G4double Tmax)
{
  G4int i;
  G4double energy1, energy2, result=0.;

  fBetaGammaSq = bg2;
  fTmax = Tmax;
	         	
  for( i = fIntervalNumber-1; i >= 0; i-- )
  {
    if( Tmax >= (*(*fMatSandiaMatrix)[i])[0] ) break;
  }
  if (i < 0) i = 0; // Tmax should be more than 
                    // first ionisation potential
  fIntervalTmax = i;

  G4Integrator<G4InitXscPAI,G4double(G4InitXscPAI::*)(G4double)> integral;

  for( i = fIntervalTmax; i >= 0; i-- )
  {
    fCurrentInterval = i;

    if(omega >= (*(*fMatSandiaMatrix)[i])[0])
    {
      if(i == fIntervalTmax) 
      {
        return integral.Legendre10(this,&G4InitXscPAI::DifPAIxSection,omega,fTmax);
      }
      else 
      {
        result += integral.Legendre10(this,&G4InitXscPAI::DifPAIxSection,omega,
                     (*(*fMatSandiaMatrix)[i+1])[0]);
        break;
      }
    }
    else
    {
      if(i == fIntervalTmax) 
      {
        result += integral.Legendre10(this,&G4InitXscPAI::DifPAIxSection,
                (*(*fMatSandiaMatrix)[i])[0],fTmax);
      }
      else
      {
        energy1 = (*(*fMatSandiaMatrix)[i])[0];
        energy2 = (*(*fMatSandiaMatrix)[i+1])[0];
        result += integral.Legendre10(this,&G4InitXscPAI::DifPAIxSection,
                           energy1,energy2);
      }
    }
    // G4cout<<"IntegralPAIdNdx<<"("<<omega<<")"<<" = "<<result<<G4endl;
  }
  return result;
}


////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI integral dEdx
// = mean energy loss per unit length, keV/cm
// 

G4double G4InitXscPAI::IntegralPAIdEdx(G4double omega, G4double bg2, G4double Tmax)
{
  G4int i;
  G4double energy1, energy2, result=0.;

  fBetaGammaSq = bg2;
  fTmax = Tmax;
	         	
  for( i = fIntervalNumber-2; i >= 0; i-- )
  {
    if(Tmax >= (*(*fMatSandiaMatrix)[i])[0]) break;
  }
  if (i < 0) i = 0; // Tmax should be more than 
                    // first ionisation potential
  fIntervalTmax = i;

  G4Integrator<G4InitXscPAI,G4double(G4InitXscPAI::*)(G4double)> integral;

  for( i = fIntervalTmax; i >= 0; i-- )
  {
    fCurrentInterval = i;

    if(omega >= (*(*fMatSandiaMatrix)[i])[0])
    {
      if(i == fIntervalTmax) 
      {
        return integral.Legendre10(this,&G4InitXscPAI::DifPAIdEdx,omega,fTmax);
      }
      else 
      {
        result += integral.Legendre10(this,&G4InitXscPAI::DifPAIdEdx,omega,
                     (*(*fMatSandiaMatrix)[i+1])[0]);
        break;
      }
    }
    else
    {
      if(i == fIntervalTmax) 
      {
        result += integral.Legendre10(this,&G4InitXscPAI::DifPAIdEdx,
                (*(*fMatSandiaMatrix)[i])[0],fTmax);
      }
      else
      {
        energy1 = (*(*fMatSandiaMatrix)[i])[0];
        energy2 = (*(*fMatSandiaMatrix)[i+1])[0];
        result += integral.Legendre10(this,&G4InitXscPAI::DifPAIdEdx,
                           energy1,energy2);
      }
    }
    // G4cout<<"IntegralPAIdEdx<<"("<<omega<<")"<<" = "<<result<<G4endl;
  }
  return result;
}

////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI Cerenkov integral cross-section
// fIntegralCrenkov[1] = specific Crenkov ionisation, 1/cm
// and fIntegralCerenkov[0] = mean Cerenkov loss per cm  in keV/cm

void G4InitXscPAI::IntegralCerenkov()
{

}   // end of IntegralCerenkov 

////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI Plasmon integral cross-section
// fIntegralPlasmon[1] = splasmon primary ionisation, 1/cm
// and fIntegralPlasmon[0] = mean plasmon loss per cm  in keV/cm

void G4InitXscPAI::IntegralPlasmon()
{

}   // end of IntegralPlasmon


/////////////////////////////////////////////////////////////////////////
//
//

G4double G4InitXscPAI::GetStepEnergyLoss( G4double step )
{  
  G4double loss = 0.0 ;
  loss *= step;

  return loss ;
}

/////////////////////////////////////////////////////////////////////////
//
//

G4double G4InitXscPAI::GetStepCerenkovLoss( G4double step )
{  
  G4double loss = 0.0 ;
  loss *= step;

  return loss ;
}

/////////////////////////////////////////////////////////////////////////
//
//

G4double G4InitXscPAI::GetStepPlasmonLoss( G4double step )
{  


  G4double loss = 0.0 ;
  loss *= step;
  return loss ;
}

   
//   
// end of G4InitXscPAI implementation file 
//
////////////////////////////////////////////////////////////////////////////

