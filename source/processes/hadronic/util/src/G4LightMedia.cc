// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LightMedia.cc,v 1.2 1999-12-15 14:53:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Light Media Charge and/or Strangeness Exchange
 // J.L. Chuma, TRIUMF, 21-Feb-1997
 // Last modified: 13-Mar-1997

#include "G4LightMedia.hh"
#include "Randomize.hh"

  G4DynamicParticle *
   G4LightMedia::PionPlusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus & targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aPiZero = G4PionZero::PionZero();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    if( targetParticle->GetDefinition() == aNeutron ) {
      
      // for pi+ n reactions, change some of the elastic cross section to pi0 p
      
      const G4double cech[] = {0.33,0.27,0.29,0.31,0.27,0.18,0.13,0.10,0.09,0.07};
      G4int iplab = G4int(G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*5.0 ));
      if( G4UniformRand() > cech[iplab]/pow(atomicNumber,0.42) ) {
        G4DynamicParticle* resultant = new G4DynamicParticle;
        resultant->SetDefinition( aPiZero );
        // targetParticle->SetDefinition( aProton );
        delete targetParticle;
        return resultant;
      }
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::PionMinusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::KaonPlusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aKaonZS = G4KaonZeroShort::KaonZeroShort();
    G4ParticleDefinition* aKaonZL = G4KaonZeroLong::KaonZeroLong();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    if( targetParticle->GetDefinition() == aNeutron ) {

      // for k+ n reactions, change some of the elastic cross section to k0 p
      
      const G4double cech[] = {0.33,0.27,0.29,0.31,0.27,0.18,0.13,0.10,0.09,0.07};
      G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*5.0 ) );
      if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
        G4DynamicParticle* resultant = new G4DynamicParticle;
        if( G4UniformRand() < 0.5 )
          resultant->SetDefinition( aKaonZS );
        else
          resultant->SetDefinition( aKaonZL );
        // targetParticle->SetDefinition( aProton );
        delete targetParticle;
        return resultant;
      }
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::KaonZeroShortExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aKaonPlus = G4KaonPlus::KaonPlus();
    G4ParticleDefinition* aKaonZL = G4KaonZeroLong::KaonZeroLong();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    if( targetParticle->GetDefinition() == aProton ) {
      
      // for k0 p reactions, change some of the elastic cross section to k+ n
      
      const G4double cech[] = {0.33,0.27,0.29,0.31,0.27,0.18,0.13,0.10,0.09,0.07};
      G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*5.0 ) );
      if( G4UniformRand() > cech[iplab]/pow(atomicNumber,0.42) ) {
        G4DynamicParticle* resultant = new G4DynamicParticle;
        resultant->SetDefinition( aKaonPlus );
        // targetParticle->SetDefinition( aNeutron );
        delete targetParticle;
        return resultant;
      }
    } else {
      if( G4UniformRand() >= 0.5 ) {
        G4DynamicParticle* resultant = new G4DynamicParticle;
        resultant->SetDefinition( aKaonZL );
        delete targetParticle;
        return resultant;
      }
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::KaonZeroLongExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aKaonZS = G4KaonZeroShort::KaonZeroShort();
    
    if( G4UniformRand() >= 0.5 ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      resultant->SetDefinition( aKaonZS );
      return resultant;
    }
    return (G4DynamicParticle*)NULL;    
  }
 
  G4DynamicParticle *
   G4LightMedia::KaonMinusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::ProtonExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    if( targetParticle->GetDefinition() == aNeutron ) {
      const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.};
      G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
      if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
        G4DynamicParticle* resultant = new G4DynamicParticle;
        resultant->SetDefinition( aNeutron );
        // targetParticle->SetDefinition( aProton );
        delete targetParticle;
        return resultant;
      }
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::AntiProtonExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* anAntiNeutron = G4AntiNeutron::AntiNeutron();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    if( targetParticle->GetDefinition() == aProton ) {
      const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.};
      G4int iplab = G4int( incidentParticle->GetTotalMomentum()/GeV*10.0 );
      if( iplab > 9 )iplab = G4int( incidentParticle->GetTotalMomentum()/GeV ) + 9;
      if( iplab > 19 )iplab = 19;
      if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.75) ) {
        G4DynamicParticle* resultant = new G4DynamicParticle;
        resultant->SetDefinition( anAntiNeutron );
        // targetParticle->SetDefinition( aNeutron );
        delete targetParticle;
        return resultant;
      }
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::NeutronExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    if( targetParticle->GetDefinition() == aProton ) {
      const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.};
      G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
      if( G4UniformRand() > cech[iplab]/pow(atomicNumber,0.42) ) {
        G4DynamicParticle* resultant = new G4DynamicParticle;
        resultant->SetDefinition( aProton );
        // targetParticle->SetDefinition( aNeutron );
        delete targetParticle;
        return resultant;
      }
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::AntiNeutronExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* anAntiProton = G4AntiProton::AntiProton();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    if( targetParticle->GetDefinition() == aNeutron ) {
      const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
      G4int iplab = G4std::min( 9, G4int( incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
      if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.75) ) {
        G4DynamicParticle* resultant = new G4DynamicParticle;
        resultant->SetDefinition( anAntiProton );
        // targetParticle->SetDefinition( aProton );
        delete targetParticle;
        return resultant;
      }
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::LambdaExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aSigmaPlus = G4SigmaPlus::SigmaPlus();
    G4ParticleDefinition* aSigmaMinus = G4SigmaMinus::SigmaMinus();
    G4ParticleDefinition* aSigmaZero = G4SigmaZero::SigmaZero();
    G4ParticleDefinition* aLambda = G4Lambda::Lambda();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      G4int irn = G4int( G4UniformRand()/0.2 );
      if( targetParticle->GetDefinition() == aNeutron ) {
        
        // LN --> S0 N, LN --> S- P, LN --> N L, LN --> N S0, LN --> P S-
        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aSigmaZero );
           break;
         case 1:
           resultant->SetDefinition( aSigmaMinus );
           // targetParticle->SetDefinition( aProton );
           break;
         case 2:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 3:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
        }
      } else {  // target particle is a proton
        
        // LP --> S+ N, LP --> S0 P, LP --> P L, LP --> P S0, LP --> N S+
        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aSigmaPlus );
           // targetParticle->SetDefinition( aNeutron );
           break;
         case 1:
           resultant->SetDefinition( aSigmaZero );
           break;
         case 2:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 3:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         default:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
 G4DynamicParticle *
  G4LightMedia::AntiLambdaExchange(
   const G4DynamicParticle* incidentParticle,
   const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* anAntiSigmaPlus = G4AntiSigmaPlus::AntiSigmaPlus();
    G4ParticleDefinition* anAntiSigmaMinus = G4AntiSigmaMinus::AntiSigmaMinus();
    G4ParticleDefinition* anAntiSigmaZero = G4AntiSigmaZero::AntiSigmaZero();
    G4ParticleDefinition* anAntiLambda = G4AntiLambda::AntiLambda();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      G4int irn = G4int( G4UniformRand()/0.2 );
      if( targetParticle->GetDefinition() == aNeutron ) {
        
        // LB N --> S+B P, LB N --> S0B N, LB N --> N LB,
        // LB N --> N S0B, LB N --> P S+B
        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( anAntiSigmaPlus );
           // targetParticle->SetDefinition( aProton );
           break;
         case 1:
           resultant->SetDefinition( anAntiSigmaZero );
           break;
         case 2:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( anAntiLambda );
           break;
         case 3:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( anAntiSigmaZero );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( anAntiSigmaPlus );
           break;
        }
      } else {  // target particle is a proton
        
        // LB P --> S0B P, LB P --> S-B N, LB P --> P LB,
        // LB P --> P S0B, LB P --> N S-B
        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( anAntiSigmaZero );
           break;
         case 1:
           resultant->SetDefinition( anAntiSigmaMinus );
           // targetParticle->SetDefinition( aNeutron );
           break;
         case 2:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( anAntiLambda );
           break;
         case 3:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( anAntiSigmaZero );
           break;
         default:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( anAntiSigmaMinus );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::SigmaPlusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aLambda = G4Lambda::Lambda();
    G4ParticleDefinition* aSigmaZero = G4SigmaZero::SigmaZero();
    G4ParticleDefinition* aSigmaPlus = G4SigmaPlus::SigmaPlus();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      
      // introduce charge and strangeness exchange reactions
      
      G4int irn = G4int( G4UniformRand()/0.2 );
      if( targetParticle->GetDefinition() == aNeutron ) {
        
        //  S+ N --> S0 P, S+ N --> L P, S+ N --> N S+, S+ N --> P S0, S+ N --> P L
        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aProton );
           break;
         case 1:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aProton );
           break;
         case 2:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
         case 3:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aLambda );
           break;
        }
      } else {  // target particle is a proton
        
        // S+ P --> P S+
        
        resultant->SetDefinition( aProton );
        // targetParticle->SetDefinition( aSigmaPlus );
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::SigmaMinusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aLambda = G4Lambda::Lambda();
    G4ParticleDefinition* aSigmaZero = G4SigmaZero::SigmaZero();
    G4ParticleDefinition* aSigmaMinus = G4SigmaMinus::SigmaMinus();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      
      // introduce charge and strangeness exchange reactions
      
      G4int irn = G4int( G4UniformRand()/0.2 );
      if( targetParticle->GetDefinition() == aNeutron ) {
        
        // S- N --> N S-
        
        resultant->SetDefinition( aNeutron );
        // targetParticle->SetDefinition( aSigmaMinus );
      } else {  // target particle is a proton
        
        //  S+ N --> S0 P, S+ N --> L P, S+ N --> N S+, S+ N --> P S0, S+ N --> P L
        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aNeutron );
           break;
         case 1:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aNeutron );
           break;
         case 2:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
         case 3:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         default:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aLambda );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::AntiSigmaPlusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* anAntiLambda = G4AntiLambda::AntiLambda();
    G4ParticleDefinition* anAntiSigmaZero = G4AntiSigmaZero::AntiSigmaZero();
    G4ParticleDefinition* anAntiSigmaPlus = G4AntiSigmaPlus::AntiSigmaPlus();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      G4int irn = G4int( G4UniformRand()/0.2 );
      if( targetParticle->GetDefinition() == aNeutron ) {
        
        // S+B N --> N S+B
        
        resultant->SetDefinition( aNeutron );
        // targetParticle->SetDefinition( anAntiSigmaPlus );
      } else {  // target particle is a proton
        
        // S+ N --> S0 P, S+ N --> L P, S+ N --> N S+, S+ N --> P S0, S+ N --> P L
        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( anAntiLambda );
           // targetParticle->SetDefinition( aNeutron );
           break;
         case 1:
           resultant->SetDefinition( anAntiSigmaZero );
           // targetParticle->SetDefinition( aNeutron );
           break;
         case 2:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( anAntiLambda );
           break;
         case 3:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( anAntiSigmaZero );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( anAntiLambda );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::AntiSigmaMinusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* anAntiLambda = G4AntiLambda::AntiLambda();
    G4ParticleDefinition* anAntiSigmaZero = G4AntiSigmaZero::AntiSigmaZero();
    G4ParticleDefinition* anAntiSigmaMinus = G4AntiSigmaMinus::AntiSigmaMinus();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      G4int irn = G4int( G4UniformRand()/0.2 );
      if( targetParticle->GetDefinition() == aNeutron ) {
        
        // S-B N --> LB P, S-B N --> S0B P, S-B N --> N S-B,
        // S-B N --> P LB, S-B N --> P S0B        
        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( anAntiLambda );
           // targetParticle->SetDefinition( aProton );
           break;
         case 1:
           resultant->SetDefinition( anAntiSigmaZero );
           // targetParticle->SetDefinition( aProton );
           break;
         case 2:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( anAntiSigmaMinus );
           break;
         case 3:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( anAntiLambda );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( anAntiSigmaZero );
           break;
        }
      } else {  // target particle is a proton
        
        // S-B P --> P S-B
        
        resultant->SetDefinition( aProton );
        // targetParticle->SetDefinition( anAntiSigmaMinus );
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::XiZeroExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aLambda = G4Lambda::Lambda();
    G4ParticleDefinition* aSigmaZero = G4SigmaZero::SigmaZero();
    G4ParticleDefinition* aSigmaMinus = G4SigmaMinus::SigmaMinus();
    G4ParticleDefinition* aSigmaPlus = G4SigmaPlus::SigmaPlus();
    G4ParticleDefinition* aXiMinus = G4XiMinus::XiMinus();
    G4ParticleDefinition* aXiZero = G4XiZero::XiZero();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      if( targetParticle->GetDefinition() == aNeutron ) {
        G4int irn = G4int( G4UniformRand()*7.0 );
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 1:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 2:
           resultant->SetDefinition( aXiMinus );
           // targetParticle->SetDefinition( aProton );
           break;
         case 3:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         case 4:
           resultant->SetDefinition( aSigmaPlus );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
         case 5:
           resultant->SetDefinition( aSigmaMinus );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
         default:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aXiZero );
           break;
        }
      } else {  // target particle is a proton
        G4int irn = G4int( G4UniformRand()*5.0 );
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aSigmaPlus );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 1:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
         case 2:
           resultant->SetDefinition( aSigmaPlus );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aXiZero );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::XiMinusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aLambda = G4Lambda::Lambda();
    G4ParticleDefinition* aSigmaZero = G4SigmaZero::SigmaZero();
    G4ParticleDefinition* aSigmaMinus = G4SigmaMinus::SigmaMinus();
    G4ParticleDefinition* aXiMinus = G4XiMinus::XiMinus();
    G4ParticleDefinition* aXiZero = G4XiZero::XiZero();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      if( targetParticle->GetDefinition() == aNeutron ) {
        G4int irn = G4int( G4UniformRand()*5.0 );
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         case 1:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
         case 2:
           resultant->SetDefinition( aSigmaMinus );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
         default:
           resultant->SetDefinition( aSigmaMinus );
           // targetParticle->SetDefinition( aLambda );
           break;
        }
      } else {  // target particle is a proton
        G4int irn = G4int( G4UniformRand()*7.0 );        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aXiZero );
           // targetParticle->SetDefinition( aNeutron );
           break;
         case 1:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aXiZero );
           break;
         case 2:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 4:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 5:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aXiMinus );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::AntiXiZeroExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    // NOTE:  The FORTRAN version of the cascade, CASAXO, simply called the
    //        routine for the XiZero particle.  Hence, the Exchange function
    //        below is just a copy of the Exchange from the XiZero particle
 
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aLambda = G4Lambda::Lambda();
    G4ParticleDefinition* aSigmaZero = G4SigmaZero::SigmaZero();
    G4ParticleDefinition* aSigmaMinus = G4SigmaMinus::SigmaMinus();
    G4ParticleDefinition* aSigmaPlus = G4SigmaPlus::SigmaPlus();
    G4ParticleDefinition* aXiMinus = G4XiMinus::XiMinus();
    G4ParticleDefinition* aXiZero = G4XiZero::XiZero();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      if( targetParticle->GetDefinition() == aNeutron ) {
        G4int irn = G4int( G4UniformRand()*7.0 );
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 1:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 2:
           resultant->SetDefinition( aXiMinus );
           // targetParticle->SetDefinition( aProton );
           break;
         case 3:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         case 4:
           resultant->SetDefinition( aSigmaPlus );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
         case 5:
           resultant->SetDefinition( aSigmaMinus );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
         default:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aXiZero );
           break;
        }
      } else {  // target particle is a proton
        G4int irn = G4int( G4UniformRand()*5.0 );
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aSigmaPlus );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 1:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
         case 2:
           resultant->SetDefinition( aSigmaPlus );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aXiZero );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::AntiXiMinusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    // NOTE:  The FORTRAN version of the cascade, CASAXM, simply called the
    //        routine for the XiMinus particle.  Hence, the Exchange function
    //        below is just a copy of the Exchange from the XiMinus particle
 
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aLambda = G4Lambda::Lambda();
    G4ParticleDefinition* aSigmaZero = G4SigmaZero::SigmaZero();
    G4ParticleDefinition* aSigmaMinus = G4SigmaMinus::SigmaMinus();
    G4ParticleDefinition* aXiMinus = G4XiMinus::XiMinus();
    G4ParticleDefinition* aXiZero = G4XiZero::XiZero();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      if( targetParticle->GetDefinition() == aNeutron ) {
        G4int irn = G4int( G4UniformRand()*5.0 );
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         case 1:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
         case 2:
           resultant->SetDefinition( aSigmaMinus );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
         default:
           resultant->SetDefinition( aSigmaMinus );
           // targetParticle->SetDefinition( aLambda );
           break;
        }
      } else {  // target particle is a proton
        G4int irn = G4int( G4UniformRand()*7.0 );        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aXiZero );
           // targetParticle->SetDefinition( aNeutron );
           break;
         case 1:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( aXiZero );
           break;
         case 2:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 4:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 5:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( aXiMinus );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::OmegaMinusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aLambda = G4Lambda::Lambda();
    G4ParticleDefinition* aSigmaZero = G4SigmaZero::SigmaZero();
    G4ParticleDefinition* aSigmaMinus = G4SigmaMinus::SigmaMinus();
    G4ParticleDefinition* aSigmaPlus = G4SigmaPlus::SigmaPlus();
    G4ParticleDefinition* aXiMinus = G4XiMinus::XiMinus();
    G4ParticleDefinition* aXiZero = G4XiZero::XiZero();
    G4ParticleDefinition* anOmegaMinus = G4OmegaMinus::OmegaMinus();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      
      // introduce charge and strangeness exchange reactions
      
      if( targetParticle->GetDefinition() == aNeutron ) {
        G4int irn = G4int( G4UniformRand()*7.0 );
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aXiZero );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
         case 1:
           resultant->SetDefinition( aSigmaMinus );
           // targetParticle->SetDefinition( aXiZero );
           break;
         case 2:
           resultant->SetDefinition( aXiMinus );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         case 4:
           resultant->SetDefinition( aXiMinus );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 5:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         default:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( anOmegaMinus );
           break;
        }
      } else {  // target particle is a proton
        G4int irn = G4int( G4UniformRand()*7.0 );        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aXiZero );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 1:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aXiZero );
           break;
         case 2:
           resultant->SetDefinition( aXiZero );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aXiZero );
           break;
         case 4:
           resultant->SetDefinition( aXiMinus );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
         case 5:
           resultant->SetDefinition( aSigmaPlus );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( anOmegaMinus );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }
 
  G4DynamicParticle *
   G4LightMedia::AntiOmegaMinusExchange(
    const G4DynamicParticle* incidentParticle,
    const G4Nucleus& targetNucleus )
  {
    // NOTE:  The FORTRAN version of the cascade, CASAOM, simply called the
    //        routine for the OmegaMinus particle.  Hence, the Exchange function
    //        below is just a copy of the Exchange from the OmegaMinus particle.
    
    G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition* aProton = G4Proton::Proton();
    G4ParticleDefinition* aLambda = G4Lambda::Lambda();
    G4ParticleDefinition* aSigmaZero = G4SigmaZero::SigmaZero();
    G4ParticleDefinition* aSigmaMinus = G4SigmaMinus::SigmaMinus();
    G4ParticleDefinition* aSigmaPlus = G4SigmaPlus::SigmaPlus();
    G4ParticleDefinition* aXiMinus = G4XiMinus::XiMinus();
    G4ParticleDefinition* aXiZero = G4XiZero::XiZero();
    G4ParticleDefinition* anOmegaMinus = G4OmegaMinus::OmegaMinus();
    
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    
    G4DynamicParticle* targetParticle = targetNucleus.ReturnTargetParticle();
    
    const G4double cech[] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
    G4int iplab = G4int( G4std::min( 9.0, incidentParticle->GetTotalMomentum()/GeV*2.5 ) );
    if( G4UniformRand() <= cech[iplab]/pow(atomicNumber,0.42) ) {
      G4DynamicParticle* resultant = new G4DynamicParticle;
      
      // introduce charge and strangeness exchange reactions
      
      if( targetParticle->GetDefinition() == aNeutron ) {
        G4int irn = G4int( G4UniformRand()*7.0 );
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aXiZero );
           // targetParticle->SetDefinition( aSigmaMinus );
           break;
         case 1:
           resultant->SetDefinition( aSigmaMinus );
           // targetParticle->SetDefinition( aXiZero );
           break;
         case 2:
           resultant->SetDefinition( aXiMinus );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         case 4:
           resultant->SetDefinition( aXiMinus );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 5:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         default:
           resultant->SetDefinition( aNeutron );
           // targetParticle->SetDefinition( anOmegaMinus );
           break;
        }
      } else {  // target particle is a proton
        G4int irn = G4int( G4UniformRand()*7.0 );        
        switch( irn ) {
         case 0:
           resultant->SetDefinition( aXiZero );
           // targetParticle->SetDefinition( aSigmaZero );
           break;
         case 1:
           resultant->SetDefinition( aSigmaZero );
           // targetParticle->SetDefinition( aXiZero );
           break;
         case 2:
           resultant->SetDefinition( aXiZero );
           // targetParticle->SetDefinition( aLambda );
           break;
         case 3:
           resultant->SetDefinition( aLambda );
           // targetParticle->SetDefinition( aXiZero );
           break;
         case 4:
           resultant->SetDefinition( aXiMinus );
           // targetParticle->SetDefinition( aSigmaPlus );
           break;
         case 5:
           resultant->SetDefinition( aSigmaPlus );
           // targetParticle->SetDefinition( aXiMinus );
           break;
         default:
           resultant->SetDefinition( aProton );
           // targetParticle->SetDefinition( anOmegaMinus );
           break;
        }
      }
      delete targetParticle;
      return resultant;
    }
    delete targetParticle;
    return (G4DynamicParticle*)NULL;
  }

 /* end of file */
 
