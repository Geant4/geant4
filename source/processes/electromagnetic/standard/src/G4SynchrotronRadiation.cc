// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SynchrotronRadiation.cc,v 1.3 2000-11-01 15:30:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, 
//      21-5-98 V.Grichine
//      
//                    
// 
//
//
// --------------------------------------------------------------

#include "G4SynchrotronRadiation.hh"
 
///////////////////////////////////////////////////////////////////////    
//
//  Constructor
//
 
G4SynchrotronRadiation::G4SynchrotronRadiation(const G4String& processName)

  : G4VDiscreteProcess("SynchrotronRadiation"),      // initialization
    LowestKineticEnergy (10.*keV),
    HighestKineticEnergy (100.*TeV),
    TotBin(200),
    theGamma (G4Gamma::Gamma() ),
    theElectron ( G4Electron::Electron() ),
    thePositron ( G4Positron::Positron() )
{
 ;
}
 
/////////////////////////////////////////////////////////////////////////
//
// Destructor
//
 
G4SynchrotronRadiation::~G4SynchrotronRadiation()
{
     ;
}
 
 
/////////////////////////////// METHODS /////////////////////////////////
//
//
// Production of synchrotron X-ray photon
// GEANT4 internal units.
// 

G4VParticleChange* 
G4SynchrotronRadiation::PostStepDoIt(const G4Track& trackData,
                                     const G4Step&  stepData      )

{
  G4int i ;
  aParticleChange.Initialize(trackData);
  G4Material* aMaterial=trackData.GetMaterial() ;

  const G4DynamicParticle* aDynamicParticle=trackData.GetDynamicParticle();

  G4double gamma = aDynamicParticle->GetTotalEnergy()/
                   (aDynamicParticle->GetMass()              ) ;
  if(gamma <= 1.0e3) 
  {
    return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
  }
  G4TransportationManager* transportMgr;
  G4FieldManager* globalFieldMgr;

  transportMgr = G4TransportationManager::GetTransportationManager() ;
  
  globalFieldMgr = transportMgr->GetFieldManager() ;

  G4bool FieldExists = globalFieldMgr->DoesFieldExist() ;
  G4ThreeVector  FieldValue;
  const G4Field*   pField = 0 ;
  if (FieldExists)
  {  
    pField = globalFieldMgr->GetDetectorField() ;
    G4ThreeVector  globPosition = trackData.GetPosition() ;
    G4double  globPosVec[3], FieldValueVec[3] ;
    globPosVec[0] = globPosition.x() ;
    globPosVec[1] = globPosition.y() ;
    globPosVec[2] = globPosition.z() ;

    pField->GetFieldValue( globPosVec, FieldValueVec ) ;
    FieldValue = G4ThreeVector( FieldValueVec[0], 
                                   FieldValueVec[1], 
                                   FieldValueVec[2]   ) ;
       
    G4ThreeVector unitMomentum = aDynamicParticle->GetMomentumDirection(); 
    G4ThreeVector unitMcrossB = FieldValue.cross(unitMomentum) ;
    G4double perpB = unitMcrossB.mag() ;
    if(perpB > 0.0)
    {
      // M-C of synchrotron photon energy

      G4double random = G4UniformRand() ;
      for(i=0;i<200;i++)
      {
        if(random >= fIntegralProbabilityOfSR[i]) break ;
      }
      G4double energyOfSR = 0.0001*i*i*fEnergyConst*gamma*gamma*perpB ;

      // check against insufficient energy

      if(energyOfSR <= 0.0)
      {
        return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
      }
      G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();
      G4ParticleMomentum 
      particleDirection = aDynamicParticle->GetMomentumDirection();

      // Gamma production cut in this material
      G4double 
      gammaEnergyCut = (G4Gamma::GetCutsInEnergy())[aMaterial->GetIndex()];

      if (energyOfSR <= gammaEnergyCut)  // only energy loss, no SR photon
      {
        aParticleChange.SetMomentumChange( particleDirection );
        aParticleChange.SetEnergyChange( kineticEnergy - energyOfSR );
        aParticleChange.SetLocalEnergyDeposit (0.); 
        aParticleChange.SetNumberOfSecondaries(0);
        return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
      }
      // SR photon has energy more than the current material cut
      // M-C of its direction
      
      G4double Teta = G4UniformRand()/gamma ;    // Very roughly

      G4double Phi  = twopi * G4UniformRand() ;

      G4double dirx = sin(Teta)*cos(Phi) , 
               diry = sin(Teta)*sin(Phi) , 
               dirz = cos(Teta) ;

      G4ThreeVector gammaDirection ( dirx, diry, dirz);
      gammaDirection.rotateUz(particleDirection);   
 
      // create G4DynamicParticle object for the SR photon
 
      G4DynamicParticle* aGamma= new G4DynamicParticle ( G4Gamma::Gamma(),
                                                         gammaDirection, 
                                                         energyOfSR      );

      aParticleChange.SetNumberOfSecondaries(1);
      aParticleChange.AddSecondary(aGamma); 

      // Update the incident particle 
   
      G4double newKinEnergy = kineticEnergy - energyOfSR ;
      
      if (newKinEnergy > 0.)
      {
        aParticleChange.SetMomentumChange( particleDirection );
        aParticleChange.SetEnergyChange( newKinEnergy );
        aParticleChange.SetLocalEnergyDeposit (0.); 
      } 
      else
      { 
        aParticleChange.SetEnergyChange( 0. );
        aParticleChange.SetLocalEnergyDeposit (0.);
        G4double charge = aDynamicParticle->GetDefinition()->GetPDGCharge();   
        if (charge<0.)
        {
           aParticleChange.SetStatusChange(fStopAndKill) ;
	}
        else      
        {
          aParticleChange.SetStatusChange(fStopButAlive) ;
	}
      } 
    }       
    else
    {
      return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
    }
  }     
  return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
}


G4double 
G4SynchrotronRadiation::GetPhotonEnergy( const G4Track& trackData,
                                         const G4Step&  stepData      )

{
  G4int i ;
  G4double energyOfSR = -1.0 ;
  G4Material* aMaterial=trackData.GetMaterial() ;

  const G4DynamicParticle* aDynamicParticle=trackData.GetDynamicParticle();

  G4double gamma = aDynamicParticle->GetTotalEnergy()/
                   (aDynamicParticle->GetMass()              ) ;

  G4TransportationManager* transportMgr;
  G4FieldManager* globalFieldMgr;

  transportMgr = G4TransportationManager::GetTransportationManager() ;
  
  globalFieldMgr = transportMgr->GetFieldManager() ;

  G4bool FieldExists = globalFieldMgr->DoesFieldExist() ;
  G4ThreeVector  FieldValue;
  const G4Field*   pField = 0 ;
  if (FieldExists)
  {  
    pField = globalFieldMgr->GetDetectorField() ;
    G4ThreeVector  globPosition = trackData.GetPosition() ;
    G4double  globPosVec[3], FieldValueVec[3] ;
    globPosVec[0] = globPosition.x() ;
    globPosVec[1] = globPosition.y() ;
    globPosVec[2] = globPosition.z() ;

    pField->GetFieldValue( globPosVec, FieldValueVec ) ;
    FieldValue = G4ThreeVector( FieldValueVec[0], 
                                   FieldValueVec[1], 
                                   FieldValueVec[2]   ) ;
       
    G4ThreeVector unitMomentum = aDynamicParticle->GetMomentumDirection(); 
    G4ThreeVector unitMcrossB = FieldValue.cross(unitMomentum) ;
    G4double perpB = unitMcrossB.mag() ;
    if(perpB > 0.0)
    {
      // M-C of synchrotron photon energy

      G4double random = G4UniformRand() ;
      for(i=0;i<200;i++)
      {
        if(random >= fIntegralProbabilityOfSR[i]) break ;
      }
      energyOfSR = 0.0001*i*i*fEnergyConst*gamma*gamma*perpB ;

      // check against insufficient energy

      if(energyOfSR <= 0.0)
      {
        return -1.0 ;
      }
      G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();
      G4ParticleMomentum 
      particleDirection = aDynamicParticle->GetMomentumDirection();

      // Gamma production cut in this material
      G4double 
      gammaEnergyCut = (G4Gamma::GetCutsInEnergy())[aMaterial->GetIndex()];

      // SR photon has energy more than the current material cut
      // M-C of its direction
      
      G4double Teta = G4UniformRand()/gamma ;    // Very roughly

      G4double Phi  = twopi * G4UniformRand() ;
    }       
    else
    {
      return -1.0 ;
    }
  }     
  return energyOfSR ;
}




////////////////////////////////////////////////////////////////////
//
// Constant for calculation of mean free path
//

const G4double
G4SynchrotronRadiation::fLambdaConst = sqrt(3.0)*electron_mass_c2/
                                       (2.5*fine_structure_const*eplus*c_light) ;

/////////////////////////////////////////////////////////////////////
//
// Constant for calculation of characterictic energy
//

const G4double
G4SynchrotronRadiation::fEnergyConst = 1.5*c_light*c_light*eplus*hbar_Planck/
                                       electron_mass_c2  ;

////////////////////////////////////////////////////////////////////
//
// Array of integral probability of synchrotron photons:
//
// the corresponding energy = 0.0001*i*i*(characteristic energy)
//

const G4double
G4SynchrotronRadiation::fIntegralProbabilityOfSR[200] =
{
  1.000000e+00,	9.428859e-01,	9.094095e-01,	8.813971e-01,	8.565154e-01,
  8.337008e-01,	8.124961e-01,	7.925217e-01,	7.735517e-01,	7.554561e-01,
  7.381233e-01,	7.214521e-01,	7.053634e-01,	6.898006e-01,	6.747219e-01,
  6.600922e-01,	6.458793e-01,	6.320533e-01,	6.185872e-01,	6.054579e-01,
  5.926459e-01,	5.801347e-01,	5.679103e-01,	5.559604e-01,	5.442736e-01,
  5.328395e-01,	5.216482e-01,	5.106904e-01,	4.999575e-01,	4.894415e-01,
  4.791351e-01,	4.690316e-01,	4.591249e-01,	4.494094e-01,	4.398800e-01,
  4.305320e-01,	4.213608e-01,	4.123623e-01,	4.035325e-01,	3.948676e-01,
  3.863639e-01,	3.780179e-01,	3.698262e-01,	3.617858e-01,	3.538933e-01,
  3.461460e-01,	3.385411e-01,	3.310757e-01,	3.237474e-01,	3.165536e-01,
  3.094921e-01,	3.025605e-01,	2.957566e-01,	2.890784e-01,	2.825237e-01,
  2.760907e-01,	2.697773e-01,	2.635817e-01,	2.575020e-01,	2.515365e-01,
  2.456834e-01,	2.399409e-01,	2.343074e-01,	2.287812e-01,	2.233607e-01,
  2.180442e-01,	2.128303e-01,	2.077174e-01,	2.027040e-01,	1.977885e-01,
  1.929696e-01,	1.882457e-01,	1.836155e-01,	1.790775e-01,	1.746305e-01,
  1.702730e-01,	1.660036e-01,	1.618212e-01,	1.577243e-01,	1.537117e-01,
  1.497822e-01,	1.459344e-01,	1.421671e-01,	1.384791e-01,	1.348691e-01,
  1.313360e-01,	1.278785e-01,	1.244956e-01,	1.211859e-01,	1.179483e-01,
  1.147818e-01,	1.116850e-01,	1.086570e-01,	1.056966e-01,	1.028026e-01,
  9.997405e-02,	9.720975e-02,	9.450865e-02,	9.186969e-02,	8.929179e-02,
  8.677391e-02,	8.431501e-02,	8.191406e-02,	7.957003e-02,	7.728192e-02,
  7.504872e-02,	7.286944e-02,	7.074311e-02,	6.866874e-02,	6.664538e-02,
  6.467208e-02,	6.274790e-02,	6.087191e-02,	5.904317e-02,	5.726079e-02,
  5.552387e-02,	5.383150e-02,	5.218282e-02,	5.057695e-02,	4.901302e-02,
  4.749020e-02,	4.600763e-02,	4.456450e-02,	4.315997e-02,	4.179325e-02,
  4.046353e-02,	3.917002e-02,	3.791195e-02,	3.668855e-02,	3.549906e-02,
  3.434274e-02,	3.321884e-02,	3.212665e-02,	3.106544e-02,	3.003452e-02,
  2.903319e-02,	2.806076e-02,	2.711656e-02,	2.619993e-02,	2.531021e-02,
  2.444677e-02,	2.360897e-02,	2.279620e-02,	2.200783e-02,	2.124327e-02,
  2.050194e-02,	1.978324e-02,	1.908662e-02,	1.841151e-02,	1.775735e-02,
  1.712363e-02,	1.650979e-02,	1.591533e-02,	1.533973e-02,	1.478250e-02,
  1.424314e-02,	1.372117e-02,	1.321613e-02,	1.272755e-02,	1.225498e-02,
  1.179798e-02,	1.135611e-02,	1.092896e-02,	1.051609e-02,	1.011712e-02,
  9.731635e-03,	9.359254e-03,	8.999595e-03,	8.652287e-03,	8.316967e-03,
  7.993280e-03,	7.680879e-03,	7.379426e-03,	7.088591e-03,	6.808051e-03,
  6.537491e-03,	6.276605e-03,	6.025092e-03,	5.782661e-03,	5.549027e-03,
  5.323912e-03,	5.107045e-03,	4.898164e-03,	4.697011e-03,	4.503336e-03,
  4.316896e-03,	4.137454e-03,	3.964780e-03,	3.798649e-03,	3.638843e-03,
  3.485150e-03,	3.337364e-03,	3.195284e-03,	3.058715e-03,	2.927469e-03,
  2.801361e-03,	2.680213e-03,	2.563852e-03,	2.452110e-03,	2.344824e-03
} ;



// end of G4SynchrotronRadiation.cc











