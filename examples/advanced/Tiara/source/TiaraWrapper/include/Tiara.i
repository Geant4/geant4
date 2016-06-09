# $Id: Tiara.i,v 1.7 2003/12/08 17:53:26 gcosmo Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-06-00-patch-01 $
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
#include "Physics_LHEP_PRECO.hh"
#include "Physics_LHEP_PRECO_HP.hh"
#include "Physics_LHEP_LEAD_HP.hh"
#include "Physics_CASCADE_HP.hh"
#include "LHEP_BIC.hh"
#include "LHEP_BIC_HP.hh"
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

%include TiaraPhysicsList.hh
%extend TiaraPhysicsList {
  std::string  getName() const {
    return std::string("Test12");
  }
}

%include Physics_LHEP_PRECO.hh
%extend Physics_LHEP_PRECO {
  std::string  getName() const {
    return std::string("LHEP_PRECO");
  }
}

%include Physics_LHEP_PRECO_HP.hh
%extend Physics_LHEP_PRECO_HP {
  std::string  getName() const {
    return std::string("LHEP_PRECO_HP");
  }
}

%include Physics_LHEP_LEAD_HP.hh
%extend Physics_LHEP_LEAD_HP {
  std::string  getName() const {
    return std::string("LHEP_LEAD_HP");
  }
}

%include Physics_CASCADE_HP.hh
%extend Physics_CASCADE_HP {
  std::string  getName() const {
    return std::string("CASCADE_HP");
  }
}

%include LHEP_BIC.hh 
%extend LHEP_BIC {
  std::string  getName() const {
    return std::string("LHEP_BIC");
  }
}

%include LHEP_BIC_HP.hh
%extend LHEP_BIC_HP {
  std::string  getName() const {
    return std::string("LHEP_BIC_HP");
  }
}


