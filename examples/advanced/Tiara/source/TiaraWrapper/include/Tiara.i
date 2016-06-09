# $Id: Tiara.i,v 1.8 2004/12/08 15:37:14 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-07-01 $
# -------------------------------------------------------------------

%module Tiara
%{
#include "G4PVPlacement.hh"

#include "TiaraSim.hh"
#include "TiaraGeometry.hh"
#include "TiaraMaterials.hh"
#include "TiaraVSourceEnergyGenerator.hh"
#include "TiaraSampledEnergy.hh"
#include "TiaraDPSSampledEnergy.hh"
#include "TiaraFixedEnergyGenerator.hh"
#include "TiaraVDirectionGenerator.hh"
#include "TiaraIsotropicDirections.hh"
#include "TiaraPrimaryGeneratorAction.hh"
#include "TiaraCellScorer.hh"
#include "TiaraCellScorerStore.hh"
#include "TiaraPhysicsList.hh"
#include "LHEP_BIC_HP.hh"
#include "LHEP_BIC.hh"
#include "LHEP.hh"
#include "LHEP_PRECO.hh"
#include "LHEP_LEAD.hh"
#include "LHEP_PRECO_HP.hh"
#include "TiaraVisEventAction.hh"
#include "TiaraTimedEventAction.hh"
#include "TiaraMeasure.hh"
#include "TiaraVisManager.hh"
#include "TiaraTally.hh"
#include "TiaraRandom.hh"
%}

%include typemaps.i
%include std_string.i
%include std_vector.i

%template(hist_vec_dbl) std::vector< double >;
%template(particle_vec_string) std::vector< std::string >;


%include CLHEP.i
%import G4Kernel.i




%include TiaraDimensions.hh
%include TiaraMaterials.hh
%include TiaraGeometry.hh


%include TiaraMeasure.hh
%include TiaraTally.hh


%include G4VCellScorer.hh
%include TiaraCellScorer.hh

%include TiaraCellScorerStore.hh
%extend TiaraCellScorerStore {
  G4VCellScorerStore *GetG4VCellScorerStore() {
    return self;
  }
}

%include TiaraVSourceEnergyGenerator.hh
%include TiaraSampledEnergy.hh
%include TiaraDPSSampledEnergy.hh
%include TiaraFixedEnergyGenerator.hh

%include TiaraVDirectionGenerator.hh
%include TiaraIsotropicDirections.hh

%include TiaraPrimaryGeneratorAction.hh

%include TiaraSim.hh

%include TiaraVisEventAction.hh

%include TiaraTimedEventAction.hh

%include TiaraVisManager.hh

%include TiaraRandom.hh

%include LHEP_BIC_HP.hh
%template(LHEP_BIC_HP) TLHEP_BIC_HP<G4VModularPhysicsList>; 
%extend TLHEP_BIC_HP<G4VModularPhysicsList> {
  std::string  getName() const {
    return std::string("LHEP_BIC_HP");
  }
}

%include LHEP_BIC.hh
%template(LHEP_BIC) TLHEP_BIC<G4VModularPhysicsList>;
%extend TLHEP_BIC<G4VModularPhysicsList> {
  std::string  getName() const {
    return std::string("LHEP_BIC");
  }
}

%include LHEP.hh
%template(LHEP) TLHEP<G4VModularPhysicsList>;
%extend TLHEP<G4VModularPhysicsList> {
  std::string  getName() const {
    return std::string("LHEP");
  }
}

%include LHEP_PRECO.hh
%template(LHEP_PRECO) TLHEP_PRECO<G4VModularPhysicsList>;
%extend TLHEP_PRECO<G4VModularPhysicsList> {
  std::string  getName() const {
    return std::string("LHEP_PRECO");
  }
}

%include LHEP_LEAD.hh
%template(LHEP_LEAD) TLHEP_LEAD<G4VModularPhysicsList>;
%extend TLHEP_LEAD<G4VModularPhysicsList> {
  std::string  getName() const {
    return std::string("LHEP_LEAD");
  }
}

%include LHEP_PRECO_HP.hh
%template(LHEP_PRECO_HP) TLHEP_PRECO_HP<G4VModularPhysicsList>;
%extend TLHEP_PRECO_HP<G4VModularPhysicsList> {
  std::string  getName() const {
    return std::string("LHEP_PRECO_HP");
  }
}
