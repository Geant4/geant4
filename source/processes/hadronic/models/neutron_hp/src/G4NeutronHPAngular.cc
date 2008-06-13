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
// 070523 bug fix for G4FPE_DEBUG on by A. Howard ( and T. Koi)
// 080612 bug fix contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #5
//
#include "G4NeutronHPAngular.hh"

void G4NeutronHPAngular::Init(std::ifstream & aDataFile)
{
//  G4cout << "here we are entering the Angular Init"<<G4endl;
  aDataFile >> theAngularDistributionType >> targetMass;
  aDataFile >> frameFlag;
  if(theAngularDistributionType == 0)
  {
    theIsoFlag = true;
  }
  else if(theAngularDistributionType==1)
  {
    G4int nEnergy;
    aDataFile >> nEnergy;  
    theCoefficients = new G4NeutronHPLegendreStore(nEnergy);
    theCoefficients->InitInterpolation(aDataFile);
    G4double temp, energy;
    G4int tempdep, nLegendre;
    G4int i, ii;
    for (i=0; i<nEnergy; i++)
    {
      aDataFile >> temp >> energy >> tempdep >> nLegendre;
      energy *=eV;
      theCoefficients->Init(i, energy, nLegendre);
      theCoefficients->SetTemperature(i, temp);
      G4double coeff=0;
      for(ii=0; ii<nLegendre; ii++)
      {
        aDataFile >> coeff;
        theCoefficients->SetCoeff(i, ii+1, coeff);
      }
    }
  }
  else if (theAngularDistributionType==2)
  {
    G4int nEnergy;
    aDataFile >> nEnergy;
    theProbArray = new G4NeutronHPPartial(nEnergy, nEnergy);
    theProbArray->InitInterpolation(aDataFile);
    G4double temp, energy;
    G4int tempdep;
    for(G4int i=0; i<nEnergy; i++)
    {
      aDataFile >> temp >> energy >> tempdep;
      energy *= eV;
      theProbArray->SetT(i, temp);
      theProbArray->SetX(i, energy);
      theProbArray->InitData(i, aDataFile);
    }
  }
  else
  {
    theIsoFlag = false;
    G4cout << "unknown distribution found for Angular"<<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, "unknown distribution needs implementation!!!");
  }    
}

void G4NeutronHPAngular::SampleAndUpdate(G4ReactionProduct & aHadron)
{
  if(theIsoFlag)
  {
//  G4cout << "Angular result "<<aHadron.GetTotalMomentum()<<" ";
// @@@ add code for isotropic emission in CMS.
      G4double costheta = 2.*G4UniformRand()-1;
      G4double theta = std::acos(costheta);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = std::sin(theta);
      G4double en = aHadron.GetTotalMomentum();
      G4ThreeVector temp(en*sinth*std::cos(phi), en*sinth*std::sin(phi), en*std::cos(theta) );
      aHadron.SetMomentum( temp );
      aHadron.Lorentz(aHadron, -1.*theTarget);
  }
  else
  {
    if(frameFlag == 1) // LAB
    {
      G4double en = aHadron.GetTotalMomentum();
      G4ReactionProduct boosted;
      boosted.Lorentz(theNeutron, theTarget);
      G4double kineticEnergy = boosted.GetKineticEnergy();
      G4double cosTh = 0.0; 
      if(theAngularDistributionType == 1) cosTh = theCoefficients->SampleMax(kineticEnergy); 
      if(theAngularDistributionType == 2) cosTh = theProbArray->Sample(kineticEnergy); 
      G4double theta = std::acos(cosTh);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = std::sin(theta);
      G4ThreeVector temp(en*sinth*std::cos(phi), en*sinth*std::sin(phi), en*std::cos(theta) );
      aHadron.SetMomentum( temp );
    }
    else if(frameFlag == 2) // costh in CMS
    {
      G4ReactionProduct boostedN;
      boostedN.Lorentz(theNeutron, theTarget);
      G4double kineticEnergy = boostedN.GetKineticEnergy();

      G4double cosTh = 0.0; 
      if(theAngularDistributionType == 1) cosTh = theCoefficients->SampleMax(kineticEnergy); 
      if(theAngularDistributionType == 2) cosTh = theProbArray->Sample(kineticEnergy); 

//080612TK bug fix contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) 
/*
    if(theAngularDistributionType == 1) // LAB
    {
      G4double en = aHadron.GetTotalMomentum();
      G4ReactionProduct boosted;
      boosted.Lorentz(theNeutron, theTarget);
      G4double kineticEnergy = boosted.GetKineticEnergy();
      G4double cosTh = theCoefficients->SampleMax(kineticEnergy);
      G4double theta = std::acos(cosTh);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = std::sin(theta);
      G4ThreeVector temp(en*sinth*std::cos(phi), en*sinth*std::sin(phi), en*std::cos(theta) );
      aHadron.SetMomentum( temp );
    }
    else if(theAngularDistributionType == 2) // costh in CMS
    {
*/

//      G4ReactionProduct boostedN;
//      boostedN.Lorentz(theNeutron, theTarget);
//      G4double kineticEnergy = boostedN.GetKineticEnergy();
//      G4double cosTh = theProbArray->Sample(kineticEnergy); 

      G4double theta = std::acos(cosTh);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = std::sin(theta);      
      
      G4ThreeVector temp(sinth*std::cos(phi), sinth*std::sin(phi), std::cos(theta) ); //CMS

//080612TK bug fix contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #5
/*
      G4double en = aHadron.GetTotalEnergy(); // Target rest
      
      // get trafo from Target rest frame to CMS
      G4ReactionProduct boostedT;
      boostedT.Lorentz(theTarget, theTarget);
      
      G4ThreeVector the3Neutron = boostedN.GetMomentum();
      G4double nEnergy = boostedN.GetTotalEnergy();
      G4ThreeVector the3Target = boostedT.GetMomentum();
      G4double tEnergy = boostedT.GetTotalEnergy();
      G4double totE = nEnergy+tEnergy;
      G4ThreeVector the3trafo = -the3Target-the3Neutron;
      G4ReactionProduct trafo; // for transformation from CMS to target rest frame
      trafo.SetMomentum(the3trafo);
      G4double cmsMom = std::sqrt(the3trafo*the3trafo);
      G4double sqrts = std::sqrt((totE-cmsMom)*(totE+cmsMom));
      trafo.SetMass(sqrts);
      trafo.SetTotalEnergy(totE);
      
      G4double gamma = trafo.GetTotalEnergy()/trafo.GetMass();
      G4double cosalpha = temp*trafo.GetMomentum()/trafo.GetTotalMomentum()/temp.mag();
      G4double fac = cosalpha*trafo.GetTotalMomentum()/trafo.GetMass();
      fac*=gamma;
      
      G4double mom;
//    For G4FPE_DEBUG ON
      G4double mom2 = ( en*fac*en*fac - 
                   (fac*fac - gamma*gamma)*
                   (en*en - gamma*gamma*aHadron.GetMass()*aHadron.GetMass())
                );
      if ( mom2 > 0.0 ) 
        mom = std::sqrt( mom2 );
      else
        mom = 0.0; 

      mom = -en*fac - mom;
      mom /= (fac*fac-gamma*gamma);
      temp = mom*temp;
      
      aHadron.SetMomentum( temp ); // now all in CMS
      aHadron.SetTotalEnergy( std::sqrt( mom*mom + aHadron.GetMass()*aHadron.GetMass() ) );
      aHadron.Lorentz(aHadron, trafo); // now in target rest frame
*/
      // Determination of the hadron kinetic energy in CMS
      // aHadron.GetKineticEnergy() is actually the residual kinetic energy in CMS (or target frame)
      // kineticEnergy is incident neutron kinetic energy  in CMS (or target frame)  
      G4double QValue = aHadron.GetKineticEnergy() - kineticEnergy;
      G4double A1     =   theTarget.GetMass()/boostedN.GetMass(); 
      G4double A1prim =   aHadron.GetMass()/ boostedN.GetMass();
      G4double kinE   = (A1+1-A1prim)/(A1+1)/(A1+1)*(A1*kineticEnergy+(1+A1)*QValue);
      G4double totalE = kinE + aHadron.GetMass();
      G4double mom2   = totalE*totalE - aHadron.GetMass()*aHadron.GetMass();
      G4double mom;
      if ( mom2 > 0.0 ) mom = std::sqrt( mom2 );
      else mom = 0.0;     

      aHadron.SetMomentum( mom*temp ); // Set momentum in CMS
      aHadron.SetKineticEnergy(kinE);  // Set kinetic energy in CMS

      // get trafo from Target rest frame to CMS
      G4ReactionProduct boostedT;
      boostedT.Lorentz(theTarget, theTarget);
      
      G4ThreeVector the3Neutron = boostedN.GetMomentum();
      G4double nEnergy = boostedN.GetTotalEnergy();
      G4ThreeVector the3Target = boostedT.GetMomentum();
      G4double tEnergy = boostedT.GetTotalEnergy();
      G4double totE = nEnergy+tEnergy;
      G4ThreeVector the3trafo = -the3Target-the3Neutron;
      G4ReactionProduct trafo; // for transformation from CMS to target rest frame
      trafo.SetMomentum(the3trafo);
      G4double cmsMom = std::sqrt(the3trafo*the3trafo);
      G4double sqrts = std::sqrt((totE-cmsMom)*(totE+cmsMom));
      trafo.SetMass(sqrts);
      trafo.SetTotalEnergy(totE);

      aHadron.Lorentz(aHadron, trafo);

    }
    else
    {
      throw G4HadronicException(__FILE__, __LINE__, "Tried to sample non isotropic neutron angular");
    }
  }
  aHadron.Lorentz(aHadron, -1.*theTarget); 
//  G4cout << aHadron.GetMomentum()<<" ";
//  G4cout << aHadron.GetTotalMomentum()<<G4endl;
}
