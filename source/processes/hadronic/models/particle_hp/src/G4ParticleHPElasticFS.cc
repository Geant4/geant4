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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 25-08-06 New Final State type (refFlag==3 , Legendre (Low Energy) + Probability (High Energy) ) 
//          is added by T. KOI
// 080904 Add Protection for negative energy results in very low energy ( 1E-6 eV ) scattering by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPElasticFS.hh"
#include "G4ParticleHPManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4IonTable.hh"
#include "G4ParticleHPDataUsed.hh"
#include "G4Pow.hh"
#include "zlib.h"

void G4ParticleHPElasticFS::Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String &, G4ParticleDefinition* )
  {
    G4String tString = "/FS";
    G4bool dbool;
    G4ParticleHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), M, dirName, tString, dbool);
    G4String filename = aFile.GetName();
    SetAZMs( A, Z, M, aFile ); 
    //theBaseA = aFile.GetA();
    //theBaseZ = aFile.GetZ();
    if(!dbool)
    {
      hasAnyData = false;
      hasFSData = false; 
      hasXsec = false;
      return;
    }
   //130205 For compressed data files 
   std::istringstream theData(std::ios::in);
   G4ParticleHPManager::GetInstance()->GetDataStream(filename,theData);
   //130205 END
    theData >> repFlag >> targetMass >> frameFlag;
    if(repFlag==1)
    {
      G4int nEnergy;
      theData >> nEnergy; 
      theCoefficients = new G4ParticleHPLegendreStore(nEnergy);
      theCoefficients->InitInterpolation(theData);
      G4double temp, energy;
      G4int tempdep, nLegendre;
      G4int i, ii;
      for (i=0; i<nEnergy; i++)
      {
        theData >> temp >> energy >> tempdep >> nLegendre;
        energy *=eV;
        theCoefficients->Init(i, energy, nLegendre);
        theCoefficients->SetTemperature(i, temp);
        G4double coeff=0;
        for(ii=0; ii<nLegendre; ii++)
        {
          // load legendre coefficients.
          theData >> coeff;
          theCoefficients->SetCoeff(i, ii+1, coeff); // @@@HPW@@@
        }
      }
    }
    else if (repFlag==2)
    {
      G4int nEnergy;
      theData >> nEnergy;
      theProbArray = new G4ParticleHPPartial(nEnergy, nEnergy);
      theProbArray->InitInterpolation(theData);
      G4double temp, energy;
      G4int tempdep, nPoints;
      for(G4int i=0; i<nEnergy; i++)
      {
        theData >> temp >> energy >> tempdep >> nPoints;
        energy *= eV;
        theProbArray->InitInterpolation(i, theData);
        theProbArray->SetT(i, temp);
        theProbArray->SetX(i, energy);
        G4double prob, costh;
        for(G4int ii=0; ii<nPoints; ii++)
        {
          // fill probability arrays.
          theData >> costh >> prob;
          theProbArray->SetX(i, ii, costh);
          theProbArray->SetY(i, ii, prob);
        }
        theProbArray->DoneSetXY( i );
      }
    }
    else if ( repFlag==3 )
    {
       G4int nEnergy_Legendre;
       theData >> nEnergy_Legendre; 
       if ( nEnergy_Legendre <= 0 ) {
          std::stringstream iss;
          iss << "G4ParticleHPElasticFS::Init Data Error repFlag is 3 but nEnergy_Legendre <= 0";
          iss << "Z, A and M of problematic file is " << theNDLDataZ << ", " << theNDLDataA << " and " << theNDLDataM << " respectively.";
          throw G4HadronicException(__FILE__, __LINE__, iss.str() );
       }
       theCoefficients = new G4ParticleHPLegendreStore( nEnergy_Legendre );
       theCoefficients->InitInterpolation( theData );
       G4double temp, energy;
       G4int tempdep, nLegendre;
       //G4int i, ii;
       for ( G4int i = 0 ; i < nEnergy_Legendre ; i++ )
       {
          theData >> temp >> energy >> tempdep >> nLegendre;
          energy *=eV;
          theCoefficients->Init( i , energy , nLegendre );
          theCoefficients->SetTemperature( i , temp );
          G4double coeff = 0;
          for (G4int ii = 0 ; ii < nLegendre ; ii++ )
          {
             // load legendre coefficients.
             theData >> coeff;
             theCoefficients->SetCoeff(i, ii+1, coeff); // @@@HPW@@@
          }
       } 

       tE_of_repFlag3 = energy; 

       G4int nEnergy_Prob;
       theData >> nEnergy_Prob;
       theProbArray = new G4ParticleHPPartial( nEnergy_Prob , nEnergy_Prob );
       theProbArray->InitInterpolation( theData );
       G4int nPoints;
       for ( G4int i=0 ; i < nEnergy_Prob ; i++ )
       {
          theData >> temp >> energy >> tempdep >> nPoints;

          energy *= eV;

//        consistency check
          if ( i == 0 )
             //if ( energy != tE_of_repFlag3 ) //110620TK This is too tight for 32bit machines 
             if ( std::abs( energy - tE_of_repFlag3 ) / tE_of_repFlag3 > 1.0e-15 )
                G4cout << "Warning Transition Energy of repFlag3 is not consistent." << G4endl; 

          theProbArray->InitInterpolation( i , theData );
          theProbArray->SetT( i , temp );
          theProbArray->SetX( i , energy );
          G4double prob, costh;
          for( G4int ii = 0 ; ii < nPoints ; ii++ )
          {
             // fill probability arrays.
             theData >> costh >> prob;
             theProbArray->SetX( i , ii , costh );
             theProbArray->SetY( i , ii , prob );
          }
          theProbArray->DoneSetXY( i );
       }
    }
    else if (repFlag==0)
    {
      theData >> frameFlag;
    }
    else
    {
      G4cout << "unusable number for repFlag: repFlag="<<repFlag<<G4endl;
      throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPElasticFS::Init -- unusable number for repFlag");
    }
   //130205 For compressed data files(theData changed from ifstream to istringstream)
   //theData.close();
  }
  G4HadFinalState * G4ParticleHPElasticFS::ApplyYourself(const G4HadProjectile & theTrack)
  {  
//    G4cout << "G4ParticleHPElasticFS::ApplyYourself+"<<G4endl;
   if ( theResult.Get() == NULL ) theResult.Put( new G4HadFinalState );
   theResult.Get()->Clear();
    G4double eKinetic = theTrack.GetKineticEnergy();
    const G4HadProjectile *incidentParticle = &theTrack;
    G4ReactionProduct theNeutron( const_cast<G4ParticleDefinition *>(incidentParticle->GetDefinition() ));
    theNeutron.SetMomentum( incidentParticle->Get4Momentum().vect() );
    theNeutron.SetKineticEnergy( eKinetic );
//    G4cout << "G4ParticleHPElasticFS::ApplyYourself++"<<eKinetic<<" "<<G4endl;
//    G4cout << "CMSVALUES 0 "<<theNeutron.GetTotalMomentum()<<G4endl;
    
    G4ReactionProduct theTarget; 
    G4Nucleus aNucleus;
    G4ThreeVector neuVelo = (1./incidentParticle->GetDefinition()->GetPDGMass())*theNeutron.GetMomentum();
    theTarget = aNucleus.GetBiasedThermalNucleus( targetMass, neuVelo, theTrack.GetMaterial()->GetTemperature());
    //t    theTarget.SetDefinition( G4IonTable::GetIonTable()->GetIon( G4int(theBaseZ), G4int(theBaseA) , 0.0 ) );  //TESTPHP
//     G4cout << "Nucleus-test"<<" "<<targetMass<<" ";
//     G4cout << theTarget.GetMomentum().x()<<" ";
//     G4cout << theTarget.GetMomentum().y()<<" ";
//     G4cout << theTarget.GetMomentum().z()<<G4endl;
    
    // neutron and target defined as reaction products.

// prepare lorentz-transformation to Lab.

    G4ThreeVector the3Neutron = theNeutron.GetMomentum();
    G4double nEnergy = theNeutron.GetTotalEnergy();
    G4ThreeVector the3Target = theTarget.GetMomentum();
//    cout << "@@@" << the3Target<<G4endl;
    G4double tEnergy = theTarget.GetTotalEnergy();
    G4ReactionProduct theCMS;
    G4double totE = nEnergy+tEnergy;
    G4ThreeVector the3CMS = the3Target+the3Neutron;
    theCMS.SetMomentum(the3CMS);
    G4double cmsMom = std::sqrt(the3CMS*the3CMS);
    G4double sqrts = std::sqrt((totE-cmsMom)*(totE+cmsMom));
    theCMS.SetMass(sqrts);
    theCMS.SetTotalEnergy(totE);
    
    // data come as fcn of n-energy in nuclear rest frame
    G4ReactionProduct boosted;
    boosted.Lorentz(theNeutron, theTarget);
    eKinetic = boosted.GetKineticEnergy(); // get kinetic energy for scattering
    G4double cosTh = -2;
    if(repFlag == 1)
    {
      cosTh = theCoefficients->SampleElastic(eKinetic);
    }
    
    else if (repFlag==2)
    {
      cosTh = theProbArray->Sample(eKinetic);
    }
    else if (repFlag==3)
    {
       if ( eKinetic <= tE_of_repFlag3 )
       {
          cosTh = theCoefficients->SampleElastic(eKinetic);
       }
       else
       {
          cosTh = theProbArray->Sample(eKinetic);
       }
    }
    else if (repFlag==0)
    {
      cosTh = 2.*G4UniformRand()-1.;
    }
    else
    {
      G4cout << "unusable number for repFlag: repFlag="<<repFlag<<G4endl;
      throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPElasticFS::Init -- unusable number for repFlag");
    }
    if(cosTh<-1.1) { return 0; }
    G4double phi = twopi*G4UniformRand();
    G4double theta = std::acos(cosTh);
    G4double sinth = std::sin(theta);
    if (frameFlag == 1) // final state data given in target rest frame.
    {
      // we have the scattering angle, now we need the energy, then do the
      // boosting.
      // relativistic elastic scattering energy angular correlation:
      theNeutron.Lorentz(theNeutron, theTarget);
      G4double e0 = theNeutron.GetTotalEnergy();
      G4double p0 = theNeutron.GetTotalMomentum();
      G4double mN = theNeutron.GetMass();
      G4double mT = theTarget.GetMass();
      G4double eE = e0+mT;
      G4double ap = (mT+eE)*(mT-eE) + (p0+mN)*(p0-mN);
      G4double a = 4*(eE+p0*cosTh)*(eE-p0*cosTh);
      G4double b = 4*ap*p0*cosTh;
      G4double c = (2.*eE*mN-ap)*(2.*eE*mN+ap);
      G4double en = (-b+std::sqrt(b*b - 4*a*c) )/(2*a);
      G4ThreeVector tempVector(en*sinth*std::cos(phi), en*sinth*std::sin(phi), en*std::cos(theta) );
      theNeutron.SetMomentum(tempVector);
      theNeutron.SetTotalEnergy(std::sqrt(en*en+theNeutron.GetMass()*theNeutron.GetMass()));
      // first to lab     
      theNeutron.Lorentz(theNeutron, -1.*theTarget);
      // now to CMS
      theNeutron.Lorentz(theNeutron, theCMS);
      theTarget.SetMomentum(-theNeutron.GetMomentum());
      theTarget.SetTotalEnergy(theNeutron.GetTotalEnergy());
      // and back to lab
      theNeutron.Lorentz(theNeutron, -1.*theCMS);
      theTarget.Lorentz(theTarget, -1.*theCMS);      
//111005 Protection for not producing 0 kinetic energy target
      if ( theNeutron.GetKineticEnergy() <= 0 ) theNeutron.SetTotalEnergy ( theNeutron.GetMass() * ( 1 + G4Pow::GetInstance()->powA( 10 , -15.65 ) ) );
      if ( theTarget.GetKineticEnergy() <= 0 ) theTarget.SetTotalEnergy ( theTarget.GetMass() * ( 1 + G4Pow::GetInstance()->powA( 10 , -15.65 ) ) );
    }
    else if (frameFlag == 2) // CMS
    {
      theNeutron.Lorentz(theNeutron, theCMS);
      theTarget.Lorentz(theTarget, theCMS);
      G4double en = theNeutron.GetTotalMomentum(); // already in CMS.
      G4ThreeVector cmsMom_tmp=theNeutron.GetMomentum(); // for neutron direction in CMS
      G4double cms_theta=cmsMom_tmp.theta();
      G4double cms_phi=cmsMom_tmp.phi();
      G4ThreeVector tempVector;
      tempVector.setX(std::cos(theta)*std::sin(cms_theta)*std::cos(cms_phi)
                      +std::sin(theta)*std::cos(phi)*std::cos(cms_theta)*std::cos(cms_phi)
                      -std::sin(theta)*std::sin(phi)*std::sin(cms_phi)  );
      tempVector.setY(std::cos(theta)*std::sin(cms_theta)*std::sin(cms_phi)
                      +std::sin(theta)*std::cos(phi)*std::cos(cms_theta)*std::sin(cms_phi)
                      +std::sin(theta)*std::sin(phi)*std::cos(cms_phi)  );
      tempVector.setZ(std::cos(theta)*std::cos(cms_theta)
                      -std::sin(theta)*std::cos(phi)*std::sin(cms_theta)  );
      tempVector *= en;
      theNeutron.SetMomentum(tempVector);
      theTarget.SetMomentum(-tempVector);
      G4double tP = theTarget.GetTotalMomentum();
      G4double tM = theTarget.GetMass();
      theTarget.SetTotalEnergy(std::sqrt((tP+tM)*(tP+tM)-2.*tP*tM));

/*
For debug purpose. 
Same transformation G4ReactionProduct.Lorentz() by 4vectors
{
G4LorentzVector n4p = G4LorentzVector ( theNeutron.GetMomentum() , theNeutron.GetKineticEnergy() + theNeutron.GetMass() );    
G4cout << "before " << ( n4p.e() - n4p.m() ) / eV<< G4endl;
G4LorentzVector cm4p = G4LorentzVector ( theCMS.GetMomentum() , theCMS.GetKineticEnergy() + theCMS.GetMass() );    
n4p.boost( cm4p.boostVector() );
G4cout << cm4p/eV << G4endl;
G4cout << "after " <<  ( n4p.e() - n4p.m() ) / eV<< G4endl;
}
*/

      theNeutron.Lorentz(theNeutron, -1.*theCMS);
//080904 Add Protection for very low energy (1e-6eV) scattering 
      if ( theNeutron.GetKineticEnergy() <= 0 )
      {
         //theNeutron.SetMomentum( G4ThreeVector(0) ); 
         //theNeutron.SetTotalEnergy ( theNeutron.GetMass() );
//110822 Protection for not producing 0 kinetic energy neutron
         theNeutron.SetTotalEnergy ( theNeutron.GetMass() * ( 1 + G4Pow::GetInstance()->powA( 10 , -15.65 ) ) );
      }

      theTarget.Lorentz(theTarget, -1.*theCMS);
//080904 Add Protection for very low energy (1e-6eV) scattering 
      if ( theTarget.GetKineticEnergy() < 0 )
      {
         //theTarget.SetMomentum( G4ThreeVector(0) ); 
         //theTarget.SetTotalEnergy ( theTarget.GetMass()  );
//110822 Protection for not producing 0 kinetic energy target
         theTarget.SetTotalEnergy ( theTarget.GetMass() * ( 1 + G4Pow::GetInstance()->powA( 10 , -15.65 ) ) );
      }
    }
    else
    {
      G4cout <<"Value of frameFlag (1=LAB, 2=CMS): "<<frameFlag;
      throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPElasticFS::ApplyYourSelf frameflag incorrect");
    }
    // now all in Lab
// nun den recoil generieren...und energy change, momentum change angeben.
    theResult.Get()->SetEnergyChange(theNeutron.GetKineticEnergy());
    theResult.Get()->SetMomentumChange(theNeutron.GetMomentum().unit());
    G4DynamicParticle* theRecoil = new G4DynamicParticle;
    theRecoil->SetDefinition( G4IonTable::GetIonTable()->GetIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA), 0 ) );
    theRecoil->SetMomentum(theTarget.GetMomentum());
    theResult.Get()->AddSecondary(theRecoil);
//    G4cout << "G4ParticleHPElasticFS::ApplyYourself 10+"<<G4endl;
    // postpone the tracking of the primary neutron
     theResult.Get()->SetStatusChange(suspend);
    return theResult.Get();
  }
