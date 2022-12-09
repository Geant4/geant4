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
//---------------------------------------------------------------------------
//
// ClassName:   G4GenericBiasingPhysics
//
// Author:      M. Verderi (Sept.10.2013)
// Modified:
// 07/11/2014, M. Verderi : fix bug of PhysicsBias(...) which was not taking
//             into account the vector of processes passed, but biasing all.
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GenericBiasingPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4BiasingHelper.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4ParallelGeometriesLimiterProcess.hh"


// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4GenericBiasingPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4GenericBiasingPhysics::G4GenericBiasingPhysics(const G4String& name)
:  G4VPhysicsConstructor(name),
  fPhysBiasAllCharged(false),    fNonPhysBiasAllCharged(false),
  fPhysBiasAllChargedISL(false), fNonPhysBiasAllChargedISL(false),
  fPhysBiasAllNeutral(false),    fNonPhysBiasAllNeutral(false),
  fPhysBiasAllNeutralISL(false), fNonPhysBiasAllNeutralISL(false),
  fVerbose(false)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4GenericBiasingPhysics::~G4GenericBiasingPhysics()
{;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::PhysicsBias(const G4String& particleName)
{
  fBiasedParticles.push_back(particleName);
  std::vector< G4String > dummy;
  fBiasedProcesses.push_back(dummy);
  fBiasAllProcesses.push_back(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::PhysicsBias(const G4String& particleName, const std::vector< G4String >& processNames)
{
  fBiasedParticles.push_back(particleName);
  fBiasedProcesses.push_back(processNames);
  fBiasAllProcesses.push_back(false); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::NonPhysicsBias(const G4String& particleName)
{
  fNonPhysBiasedParticles.push_back(particleName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::Bias(const G4String& particleName)
{
  PhysicsBias(particleName);
  NonPhysicsBias(particleName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::Bias(const G4String& particleName, const std::vector< G4String >& processNames)
{
  PhysicsBias(particleName, processNames);
  NonPhysicsBias(particleName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GenericBiasingPhysics::PhysicsBiasAddPDGRange( G4int PDGlow, G4int PDGhigh, G4bool includeAntiParticle )
{
  if ( PDGlow > PDGhigh ) G4cout << " G4GenericBiasingPhysics::PhysicsBiasAddPDGRange(...) :  PDGlow > PDGhigh, call ignored." << G4endl;
  fPhysBiasByPDGRangeLow .push_back( PDGlow  );
  fPhysBiasByPDGRangeHigh.push_back( PDGhigh );
  if ( includeAntiParticle )
    {
      fPhysBiasByPDGRangeLow .push_back( -PDGhigh );
      fPhysBiasByPDGRangeHigh.push_back( -PDGlow  );
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GenericBiasingPhysics::NonPhysicsBiasAddPDGRange( G4int PDGlow, G4int PDGhigh, G4bool includeAntiParticle )
{
  if ( PDGlow > PDGhigh ) G4cout << " G4GenericBiasingPhysics::NonPhysicsBiasAddPDGRange(...) :  PDGlow > PDGhigh, call ignored." << G4endl;
  fNonPhysBiasByPDGRangeLow .push_back( PDGlow  );
  fNonPhysBiasByPDGRangeHigh.push_back( PDGhigh );
  if ( includeAntiParticle )
    {
      fNonPhysBiasByPDGRangeLow .push_back( -PDGhigh );
      fNonPhysBiasByPDGRangeHigh.push_back( -PDGlow  );
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GenericBiasingPhysics::BiasAddPDGRange( G4int PDGlow, G4int PDGhigh, G4bool includeAntiParticle )
{
  if ( PDGlow > PDGhigh ) G4cout << " G4GenericBiasingPhysics::BiasAddPDGRange(...) :  PDGlow > PDGhigh, call ignored." << G4endl;
  PhysicsBiasAddPDGRange   ( PDGlow, PDGhigh, includeAntiParticle );
  NonPhysicsBiasAddPDGRange( PDGlow, PDGhigh, includeAntiParticle );
}

void G4GenericBiasingPhysics::PhysicsBiasAllCharged( G4bool includeShortLived )
{
  fPhysBiasAllCharged    = true; 
  fPhysBiasAllChargedISL = includeShortLived;
} 
void G4GenericBiasingPhysics::NonPhysicsBiasAllCharged( G4bool includeShortLived )
{ 
  fNonPhysBiasAllCharged    = true; 
  fNonPhysBiasAllChargedISL = includeShortLived;
}
void G4GenericBiasingPhysics::BiasAllCharged( G4bool includeShortLived )
{  
  fPhysBiasAllCharged       = true;
  fNonPhysBiasAllCharged    = true;
  fPhysBiasAllChargedISL    = includeShortLived;
  fNonPhysBiasAllChargedISL = includeShortLived;
}
void G4GenericBiasingPhysics::PhysicsBiasAllNeutral( G4bool includeShortLived )
{ 
  fPhysBiasAllNeutral    = true;  
  fPhysBiasAllNeutralISL = includeShortLived;
} 
void G4GenericBiasingPhysics::NonPhysicsBiasAllNeutral( G4bool includeShortLived )
{
  fNonPhysBiasAllNeutral    = true;
  fNonPhysBiasAllNeutralISL = includeShortLived;
}
void G4GenericBiasingPhysics::BiasAllNeutral( G4bool includeShortLived )
{  
  fPhysBiasAllNeutral       = true;
  fNonPhysBiasAllNeutral    = true;
  fPhysBiasAllNeutralISL    = includeShortLived;
  fNonPhysBiasAllNeutralISL = includeShortLived;
}


void G4GenericBiasingPhysics::AddParallelGeometry( const G4String& particleName, const G4String&                parallelGeometryName  )
{
  // -- add particle, caring of possible duplication:
  G4bool isKnown = false;
  for ( G4String knownParticle : fParticlesWithParallelGeometries )
    {
      if ( knownParticle == particleName )
	{
	  isKnown = true;
	  break;
	}
    }

  // -- add the geometry, caring for possible duplication of this geometry, for this particle:
  if ( !isKnown ) fParticlesWithParallelGeometries.push_back( particleName );
  std::vector< G4String >& geometries = fParallelGeometriesForParticle[particleName];

  isKnown = false;
  for ( G4String knownGeometry : geometries )
    {
      if ( knownGeometry == parallelGeometryName )
	{
	  isKnown = true;
	  break;
	}
    }
  if ( !isKnown ) geometries.push_back( parallelGeometryName );
  
}

void G4GenericBiasingPhysics::AddParallelGeometry( const G4String& particleName, const std::vector< G4String >& parallelGeometryNames )
{
  for ( G4String geometry : parallelGeometryNames ) AddParallelGeometry( particleName, geometry );
}

void G4GenericBiasingPhysics::AddParallelGeometry( G4int PDGlow, G4int PDGhigh,  const G4String&                parallelGeometryName , G4bool includeAntiParticle )
{
  if ( PDGlow > PDGhigh )
    {
      G4cout << "G4GenericBiasingPhysics::AddParallelGeometry( G4int PDGlow, G4int PDGhigh, const G4String& parallelGeometryName , G4bool includeAntiParticle = true ), PDGlow > PDGhigh : call ignored" << G4endl;
      return;
    }
  
  fPDGlowParallelGeometries .push_back( PDGlow );
  fPDGhighParallelGeometries.push_back( PDGhigh );
  G4int rangeIndex = G4int(fPDGlowParallelGeometries.size() - 1);
  fPDGrangeParallelGeometries[rangeIndex].push_back( parallelGeometryName );

  if ( includeAntiParticle )
    {
      fPDGlowParallelGeometries .push_back( -PDGhigh );
      fPDGhighParallelGeometries.push_back( -PDGlow  );
      rangeIndex = G4int(fPDGlowParallelGeometries.size() - 1);
      fPDGrangeParallelGeometries[rangeIndex].push_back( parallelGeometryName );
    }
  
}

void G4GenericBiasingPhysics::AddParallelGeometry( G4int PDGlow, G4int PDGhigh,  const std::vector< G4String >& parallelGeometryNames, G4bool includeAntiParticle )
{
  if ( PDGlow > PDGhigh )
    {
      G4cout << "G4GenericBiasingPhysics::AddParallelGeometry( G4int PDGlow, G4int PDGhigh, const std::vector< G4String >& parallelGeometryNames, G4bool includeAntiParticle = true ), PDGlow > PDGhigh : call ignored" << G4endl;
      return;
    }
  
  for ( G4String geometry : parallelGeometryNames ) AddParallelGeometry( PDGlow, PDGhigh, geometry, includeAntiParticle );
}

void G4GenericBiasingPhysics::AddParallelGeometryAllCharged(                     const G4String&                parallelGeometryName , G4bool includeShortLived )
{
  G4bool isKnown = false;
  for ( G4String geometry : fParallelGeometriesForCharged )
    {
      if ( geometry == parallelGeometryName )
	{
	  isKnown = true;
	  break;
	}
    }
  if ( !isKnown )
    {
      fParallelGeometriesForCharged   .push_back( parallelGeometryName );
      fAllChargedParallelGeometriesISL.push_back( includeShortLived );
    }
}

void G4GenericBiasingPhysics::AddParallelGeometryAllCharged(                     const std::vector< G4String >& parallelGeometryNames, G4bool includeShortLived )
{
  for ( G4String geometry : parallelGeometryNames ) AddParallelGeometryAllCharged( geometry, includeShortLived );
}

void G4GenericBiasingPhysics::AddParallelGeometryAllNeutral(                     const G4String&                parallelGeometryName , G4bool includeShortLived )
{
  G4bool isKnown = false;
  for ( G4String geometry : fParallelGeometriesForNeutral )
    {
      if ( geometry == parallelGeometryName )
	{
	  isKnown = true;
	  break;
	}
    }
  if ( !isKnown )
    {
      fParallelGeometriesForNeutral   .push_back( parallelGeometryName );
      fAllNeutralParallelGeometriesISL.push_back( includeShortLived );
    } 
}

void G4GenericBiasingPhysics::AddParallelGeometryAllNeutral(                     const std::vector< G4String >& parallelGeometryNames, G4bool includeShortLived )
{
  for ( G4String geometry : parallelGeometryNames ) AddParallelGeometryAllNeutral( geometry, includeShortLived ); 
}

  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::ConstructParticle()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GenericBiasingPhysics::ConstructProcess()
{
  
  // -- bias setup per individual particle name:
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  
  while( (*particleIterator)() )
    {
      G4ParticleDefinition*     particle = particleIterator->value();
      G4String              particleName = particle->GetParticleName();
      G4ProcessManager*         pmanager = particle->GetProcessManager();
      
      // -- include non physics process interface for biasing:
      if ( std::find(fNonPhysBiasedParticles.begin(),
		     fNonPhysBiasedParticles.end(),
		     particleName             )  != fNonPhysBiasedParticles.end() )
	{
	  G4BiasingHelper::ActivateNonPhysicsBiasing(pmanager);
	}
      
      // -- wrap biased physics processes, all processes or only user selected:
      std::vector< G4String >::const_iterator particleIt =
	std::find(fBiasedParticles.begin(),
		  fBiasedParticles.end(),
		  particleName             );
      if ( particleIt == fBiasedParticles.end() ) continue;
      
      std::vector < G4String >& biasedProcesses = fBiasedProcesses [ particleIt - fBiasedParticles.begin() ];
      G4bool biasAll                            = fBiasAllProcesses[ particleIt - fBiasedParticles.begin() ];
      
      if ( biasAll )
	{
	  G4ProcessVector*  vprocess = pmanager->GetProcessList();
	  for (G4int ip = 0 ; ip < (G4int)vprocess->size() ; ++ip)
	    {
	      G4VProcess* process = (*vprocess)[ip];
	      biasedProcesses.push_back( process->GetProcessName() );
	    }
	}
      
      G4bool restartLoop(true);
      while ( restartLoop )
	{
	  for (std::size_t ip = 0 ; ip < biasedProcesses.size() ; ++ip)
	    {
	      G4bool activ = G4BiasingHelper::ActivatePhysicsBiasing(pmanager, biasedProcesses[ip] );
	      restartLoop = activ;
	      if ( restartLoop ) break;
	    }
	}
      
    }

  
  // -- bias setup per group:
  particleIterator->reset();

  while( (*particleIterator)() )
    {
      G4ParticleDefinition*     particle = particleIterator->value();
      G4String              particleName = particle->GetParticleName();
      G4ProcessManager*         pmanager = particle->GetProcessManager();
      
      // -- exclude particles invidually specified by name:
      if ( std::find( fNonPhysBiasedParticles.begin(),
		      fNonPhysBiasedParticles.end(),
		      particleName                    ) != fNonPhysBiasedParticles.end() ) continue;
      
      if ( std::find( fBiasedParticles.begin(),
		      fBiasedParticles.end(),
		      particleName                    ) != fBiasedParticles.end() ) continue;


      G4bool physBias(false), nonPhysBias(false);
      
      auto PDG = particle->GetPDGEncoding();

      // -- include particle if in right PDG range:
      for ( std::size_t i = 0 ; i < fPhysBiasByPDGRangeLow.size() ; i++ )
	if ( ( PDG >= fPhysBiasByPDGRangeLow[i] ) && ( PDG <= fPhysBiasByPDGRangeHigh[i] ) )
	  {
	    physBias = true;
	    break;
	  }
      for ( std::size_t i = 0 ; i < fNonPhysBiasByPDGRangeLow.size() ; i++ )
	if ( ( PDG >= fNonPhysBiasByPDGRangeLow[i] ) && ( PDG <= fNonPhysBiasByPDGRangeHigh[i] ) )
	  {
	    nonPhysBias = true;
	    break;
	  }
      
      // -- if particle has not yet any biasing, include it on charge criteria:
      if ( ( physBias == false ) && ( nonPhysBias == false ) )
	{
	  if ( std::abs( particle->GetPDGCharge() ) > DBL_MIN )
	    {
	      if ( fPhysBiasAllCharged    ) if (    fPhysBiasAllChargedISL || !particle->IsShortLived() )    physBias = true;
	      if ( fNonPhysBiasAllCharged ) if ( fNonPhysBiasAllChargedISL || !particle->IsShortLived() ) nonPhysBias = true;
	    }
	  else
	    {
	      if ( fPhysBiasAllNeutral    ) if (    fPhysBiasAllNeutralISL || !particle->IsShortLived() )    physBias = true;
	      if ( fNonPhysBiasAllNeutral ) if ( fNonPhysBiasAllNeutralISL || !particle->IsShortLived() ) nonPhysBias = true;
	    }
	}
      
      
      if ( nonPhysBias ) G4BiasingHelper::ActivateNonPhysicsBiasing(pmanager);
      
      if ( physBias )
	{
	  std::vector < G4String > biasedProcesses;
	  G4ProcessVector*  vprocess = pmanager->GetProcessList();
	  for (G4int ip = 0 ; ip < (G4int)vprocess->size() ; ++ip)
	    {
	      G4VProcess* process = (*vprocess)[ip];
	      biasedProcesses.push_back( process->GetProcessName() );
	    }
	  
	  G4bool restartLoop(true);
	  while ( restartLoop )
	    {
	      for (std::size_t ip = 0 ; ip < biasedProcesses.size() ; ++ip)
		{
		  G4bool activ = G4BiasingHelper::ActivatePhysicsBiasing(pmanager, biasedProcesses[ip] );
		  restartLoop = activ;
		  if ( restartLoop ) break;
		}
	    }
	}
      
    }



  // -- Associate parallel geometries:
  AssociateParallelGeometries();
  

  // -- tells what is done:
  if ( fVerbose )
    {
      // -- print:
      particleIterator->reset();
      
      while( (*particleIterator)() )
	{
	  G4ParticleDefinition*     particle = particleIterator->value();
	  G4String              particleName = particle->GetParticleName();
	  G4ProcessManager*         pmanager = particle->GetProcessManager();
	  
	  G4bool isBiased(false);
	  G4String processNames;
	  G4int icount(0);
	  
	  G4ProcessVector*  vprocess = pmanager->GetProcessList();
	  for (G4int ip = 0 ; ip < (G4int)vprocess->size() ; ++ip)
	    {
	      G4VProcess* process = (*vprocess)[ip];
	      G4BiasingProcessInterface* pb = dynamic_cast< G4BiasingProcessInterface* >(process);
	      if ( pb != nullptr )
		{
		  isBiased = true;
		  if ( icount < 3 )
		    {
		      processNames += pb->GetProcessName();
		      processNames += " ";
		    }
		  else
		    {
		      processNames += "\n                     ";
		      processNames += pb->GetProcessName();
		      processNames += " ";
		      icount = 0;
		    }
		  icount++;
		}
	    }
	  if ( isBiased )
	    {
	      if ( particle->IsShortLived() )
		G4cout << std::setw(14) << particleName << " **** : " << processNames << G4endl;
	      else
		G4cout << std::setw(18) << particleName << " : " << processNames << G4endl;
	    }
	}
    }
}



void G4GenericBiasingPhysics::AssociateParallelGeometries()
{
  
  // -- parallel geometries for individual particles:
    auto particleIterator=GetParticleIterator();
    particleIterator->reset();
  
  while( (*particleIterator)() )
    {
      G4ParticleDefinition*     particle = particleIterator->value();
      G4String              particleName = particle->GetParticleName();
      G4ProcessManager*         pmanager = particle->GetProcessManager();

      G4bool requested = false;
      for ( G4String requestedParticles : fParticlesWithParallelGeometries )
	{
	  if ( requestedParticles == particleName )
	    {
	      requested = true;
	      break;
	    }
	}
      if ( requested )
	{
	  // -- insert biasing process for handling parallel geometries:
	  G4ParallelGeometriesLimiterProcess* limiter = G4BiasingHelper::AddLimiterProcess(pmanager);

	  // -- attach the requested worlds to this process:
	  std::vector< G4String >& parallelWorlds = fParallelGeometriesForParticle[ particleName ];
	  for ( G4String world : parallelWorlds ) limiter->AddParallelWorld( world );
	}
      
    }

  
  // -- parallel geometries for particles in PDG ranges:
  G4int i = 0; // -- index for PDG range
  for ( G4int PDGlow : fPDGlowParallelGeometries )
    {
      G4int PDGhigh     = fPDGhighParallelGeometries[i];
      auto & geometries = fPDGrangeParallelGeometries[i];
      
      particleIterator->reset();

      while( (*particleIterator)() )
	{
	  G4ParticleDefinition*    particle = particleIterator->value();
	  G4int                 particlePDG = particle->GetPDGEncoding();
	  G4ProcessManager*        pmanager = particle->GetProcessManager();

	  if ( ( particlePDG >= PDGlow ) && ( particlePDG <= PDGhigh ) )
	    {
	      // -- §§ exclude particles from individual list ?
	      // -- insert biasing process for handling parallel geometries:
	      G4ParallelGeometriesLimiterProcess* limiter = G4BiasingHelper::AddLimiterProcess(pmanager);

	      // -- attached the requested worlds to this process:
	      for ( auto& geometry : geometries ) limiter->AddParallelWorld( geometry );
	    }
	}
      // -- increment index for next PDG range:
      i++;
    }


  // -- parallel geometries for all neutral / charged particles:
  particleIterator->reset();
  G4bool islAllNeutral = false;
  for(auto isln : fAllNeutralParallelGeometriesISL)
  { islAllNeutral |= isln; }
  G4bool islAllCharged = false;
  for(auto islc : fAllChargedParallelGeometriesISL)
  { islAllCharged |= islc; }

  while((*particleIterator)())
  {
    G4ParticleDefinition*    particle = particleIterator->value();
    G4ProcessManager*        pmanager = particle->GetProcessManager();
    if(particle->GetPDGCharge() == 0.)
    {
      // Neutral particle
      if(particle->IsShortLived() && !islAllNeutral) continue;
      G4ParallelGeometriesLimiterProcess* limiter = G4BiasingHelper::AddLimiterProcess(pmanager);
      G4int j = 0;
      for(G4String wNameN : fParallelGeometriesForCharged)
      {
        if(!(particle->IsShortLived()) || fAllNeutralParallelGeometriesISL[j])
        { limiter->AddParallelWorld(wNameN); }
        j++;
      }
    }
    else
    {
      // charged
      if(particle->IsShortLived() && !islAllCharged) continue;
      G4ParallelGeometriesLimiterProcess* limiter = G4BiasingHelper::AddLimiterProcess(pmanager);
      G4int j = 0;
      for(G4String wNameC : fParallelGeometriesForCharged)
      {
        if(!(particle->IsShortLived()) || fAllChargedParallelGeometriesISL[j])
        { limiter->AddParallelWorld(wNameC); }
        j++;
      }
    }
  }
  
}
