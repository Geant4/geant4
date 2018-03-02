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
// $Id: G4ScoreQuantityMessenger.cc 108471 2018-02-15 14:12:06Z gcosmo $
//
// ---------------------------------------------------------------------
// Modifications
// 08-Oct-2010 T.Aso remove unit of G4PSPassageCellCurrent.
//  24-Mar-2011  T.Aso  Add StepChecker for debugging.
//  24-Mar-2011  T.Aso  Size and segmentation for replicated cylinder.
//  01-Jun-2012  T.Aso  Support weighted/dividedByArea options 
//                      in flatCurrent and flatFulx commands.
//  27-Mar-2013  T.Aso  Unit option in the kineticEnergy filter was 
//                     supported.
//                       
// ---------------------------------------------------------------------

#include "G4ScoreQuantityMessenger.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"

#include "G4PSCellCharge3D.hh"
#include "G4PSCellFlux3D.hh"
#include "G4PSCellFluxForCylinder3D.hh"
#include "G4PSPassageCellFlux3D.hh"
#include "G4PSPassageCellFluxForCylinder3D.hh"
#include "G4PSEnergyDeposit3D.hh"
#include "G4PSDoseDeposit3D.hh"
#include "G4PSDoseDepositForCylinder3D.hh"
#include "G4PSNofStep3D.hh"
#include "G4PSNofSecondary3D.hh"
//
#include "G4PSTrackLength3D.hh"
#include "G4PSPassageCellCurrent3D.hh"
#include "G4PSPassageTrackLength3D.hh"
#include "G4PSFlatSurfaceCurrent3D.hh"
#include "G4PSFlatSurfaceFlux3D.hh"
#include "G4PSSphereSurfaceCurrent3D.hh"
#include "G4PSSphereSurfaceFlux3D.hh"
#include "G4PSCylinderSurfaceCurrent3D.hh"
#include "G4PSCylinderSurfaceFlux3D.hh"
#include "G4PSNofCollision3D.hh"
#include "G4PSPopulation3D.hh"
#include "G4PSTrackCounter3D.hh"
#include "G4PSTermination3D.hh"
#include "G4PSMinKinEAtGeneration3D.hh"
//
// For debug purpose
#include "G4PSStepChecker3D.hh"

#include "G4SDChargedFilter.hh"
#include "G4SDNeutralFilter.hh"
#include "G4SDKineticEnergyFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDParticleWithEnergyFilter.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"

G4ScoreQuantityMessenger::G4ScoreQuantityMessenger(G4ScoringManager* SManager)
:fSMan(SManager)
{
  QuantityCommands();
  FilterCommands();
}

void G4ScoreQuantityMessenger::FilterCommands()
{
  G4UIparameter* param;

  //
  // Filter commands 
  filterDir = new G4UIdirectory("/score/filter/");
  filterDir->SetGuidance("  Scoring filter commands.");
  //
  fchargedCmd = new G4UIcmdWithAString("/score/filter/charged",this);
  fchargedCmd->SetGuidance("Charged particle filter.");
  fchargedCmd->SetParameterName("fname",false);
  //
  fneutralCmd = new G4UIcmdWithAString("/score/filter/neutral",this);
  fneutralCmd->SetGuidance("Neutral particle filter.");
  fneutralCmd->SetParameterName("fname",false);
  //
  fkinECmd = new G4UIcommand("/score/filter/kineticEnergy",this);
  fkinECmd->SetGuidance("Kinetic energy filter.");
  fkinECmd->SetGuidance("[usage] /score/filter/kineticEnergy fname Elow Ehigh unit");
  fkinECmd->SetGuidance("  fname     :(String) Filter Name ");
  fkinECmd->SetGuidance("  Elow      :(Double) Lower edge of kinetic energy");
  fkinECmd->SetGuidance("  Ehigh     :(Double) Higher edge of kinetic energy");
  fkinECmd->SetGuidance("  unit      :(String) unit of given kinetic energy");
  param = new G4UIparameter("fname",'s',false);
  fkinECmd->SetParameter(param);
  param = new G4UIparameter("elow",'d',true);
  param->SetDefaultValue("0.0");
  fkinECmd->SetParameter(param);
  param = new G4UIparameter("ehigh",'d',true);
  fkinECmd->SetParameter(param);
  G4String smax = DtoS(DBL_MAX);
  param->SetDefaultValue(smax);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("keV");
  fkinECmd->SetParameter(param);
  //
  fparticleCmd = new G4UIcommand("/score/filter/particle",this);
  fparticleCmd->SetGuidance("Particle filter.");
  fparticleCmd->SetGuidance("[usage] /score/filter/particle fname p0 .. pn");
  fparticleCmd->SetGuidance("  fname     :(String) Filter Name ");
  fparticleCmd->SetGuidance("  p0 .. pn  :(String) particle names");
  param = new G4UIparameter("fname",'s',false);
  fparticleCmd->SetParameter(param);
  param = new G4UIparameter("particlelist",'s',false);
  param->SetDefaultValue("");
  fparticleCmd->SetParameter(param);
  //
  //
  //
  fparticleKinECmd = new G4UIcommand("/score/filter/particleWithKineticEnergy",this);
  fparticleKinECmd->SetGuidance("Particle with kinetic energy filter.");
  fparticleKinECmd->SetGuidance("[usage] /score/filter/particleWithKineticEnergy fname Elow Ehigh unit p0 .. pn");
  fparticleKinECmd->SetGuidance("  fname     :(String) Filter Name ");
  fparticleKinECmd->SetGuidance("  Elow      :(Double) Lower edge of kinetic energy");
  fparticleKinECmd->SetGuidance("  Ehigh     :(Double) Higher edge of kinetic energy");
  fparticleKinECmd->SetGuidance("  unit      :(String) unit of given kinetic energy");
  fparticleKinECmd->SetGuidance("  p0 .. pn  :(String) particle names");
  param = new G4UIparameter("fname",'s',false);
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("elow",'d',false);
  param->SetDefaultValue("0.0");
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("ehigh",'d',true);
  param->SetDefaultValue(smax);
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("keV");
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("particlelist",'s',false);
  param->SetDefaultValue("");
  fparticleKinECmd->SetParameter(param);
  //
  //
}

G4ScoreQuantityMessenger::~G4ScoreQuantityMessenger()
{
    delete         quantityDir;
    delete         qTouchCmd;
    delete         qGetUnitCmd;
    delete         qSetUnitCmd;

    //
    delete    qCellChgCmd;
    delete    qCellFluxCmd;
    delete    qPassCellFluxCmd;
    delete    qeDepCmd;
    delete    qdoseDepCmd;
    delete    qnOfStepCmd;
    delete    qnOfSecondaryCmd;
    //
    delete          qTrackLengthCmd;
    delete          qPassCellCurrCmd;
    delete          qPassTrackLengthCmd;
    delete          qFlatSurfCurrCmd;
    delete          qFlatSurfFluxCmd;
//    delete          qSphereSurfCurrCmd;
//    delete          qSphereSurfFluxCmd;
//    delete          qCylSurfCurrCmd;
//    delete          qCylSurfFluxCmd;
    delete          qNofCollisionCmd;
    delete          qPopulationCmd;
    delete          qTrackCountCmd;
    delete          qTerminationCmd;
    delete          qMinKinEAtGeneCmd;
    //
    delete          qStepCheckerCmd;
    //
    delete   filterDir;
    delete   fchargedCmd;
    delete   fneutralCmd;
    delete   fkinECmd;
    delete   fparticleCmd;
    delete   fparticleKinECmd;
}

void G4ScoreQuantityMessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
      G4ExceptionDescription ed;
      //
      // Get Current Mesh
      //
      G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
      if(!mesh)
      {
        ed << "ERROR : No mesh is currently open. Open/create a mesh first. Command ignored.";
        command->CommandFailed(ed);
        return;
      }
      // Tokens
      G4TokenVec token;
      FillTokenVec(newVal,token);
      //
      // Commands for Current Mesh
      if(command==qTouchCmd) {
              mesh->SetCurrentPrimitiveScorer(newVal);
      } else if(command == qGetUnitCmd ){
          G4cout << "Unit:  "<< mesh->GetCurrentPSUnit() <<G4endl;
      } else if(command == qSetUnitCmd ){
          mesh->SetCurrentPSUnit(newVal);
      } else if(command== qCellChgCmd) {
          if ( CheckMeshPS(mesh,token[0],command) ){
              G4PSCellCharge3D* ps = new G4PSCellCharge3D(token[0]);
              ps->SetUnit(token[1]);
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qCellFluxCmd) {
          if ( CheckMeshPS(mesh,token[0],command) ){
            if( mesh->GetShape()==boxMesh ) {
              G4PSCellFlux3D* ps = new G4PSCellFlux3D(token[0]);
              ps->SetUnit(token[1]);
              mesh->SetPrimitiveScorer(ps);
            } else if( mesh->GetShape()==cylinderMesh ) {
              G4PSCellFluxForCylinder3D* ps = 
                new G4PSCellFluxForCylinder3D(token[0]);
              ps->SetUnit(token[1]);
              G4ThreeVector msize = mesh->GetSize(); // gevin in R Z N/A
              ps->SetCylinderSize(msize[0],msize[1]); // given in dr dz
              G4int nSeg[3];
              mesh->GetNumberOfSegments(nSeg);
              ps->SetNumberOfSegments(nSeg);
              mesh->SetPrimitiveScorer(ps);
            }
          }
      } else if(command== qPassCellFluxCmd) {
          if ( CheckMeshPS(mesh,token[0],command) ){
            if( mesh->GetShape()==boxMesh ) {
              G4PSPassageCellFlux3D* ps = new G4PSPassageCellFlux3D(token[0]);
              ps->SetUnit(token[1]);
              mesh->SetPrimitiveScorer(ps);
            } else if( mesh->GetShape()==cylinderMesh ) {
              G4PSPassageCellFluxForCylinder3D* ps = 
                new G4PSPassageCellFluxForCylinder3D(token[0]);
              ps->SetUnit(token[1]);
              G4ThreeVector msize = mesh->GetSize();  // gevin in R Z N/A
              ps->SetCylinderSize(msize[0],msize[1]); // given in dr dz
              G4int nSeg[3];
              mesh->GetNumberOfSegments(nSeg);
              ps->SetNumberOfSegments(nSeg);
              mesh->SetPrimitiveScorer(ps);
            }
          }
      } else if(command==qeDepCmd) {
          if ( CheckMeshPS(mesh,token[0],command) ){
              G4PSEnergyDeposit3D* ps =new G4PSEnergyDeposit3D(token[0]);
              ps->SetUnit(token[1]);
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qdoseDepCmd) {
          if ( CheckMeshPS(mesh,token[0],command) ){
            if( mesh->GetShape()==boxMesh ) {
              G4PSDoseDeposit3D* ps = new G4PSDoseDeposit3D(token[0]);
              ps->SetUnit(token[1]);
              mesh->SetPrimitiveScorer(ps);
            } else if( mesh->GetShape()==cylinderMesh ) {
              G4PSDoseDepositForCylinder3D* ps = 
                new G4PSDoseDepositForCylinder3D(token[0]);
              ps->SetUnit(token[1]);
              G4ThreeVector msize = mesh->GetSize(); // gevin in R Z N/A
              ps->SetCylinderSize(msize[0],msize[1]); // given in dr dz
              G4int nSeg[3];
              mesh->GetNumberOfSegments(nSeg);
              ps->SetNumberOfSegments(nSeg);
              mesh->SetPrimitiveScorer(ps);
            }
          }
      } else if(command== qnOfStepCmd) {
          if ( CheckMeshPS(mesh,token[0],command) ){
              G4PSNofStep3D* ps = new G4PSNofStep3D(token[0]);
              ps->SetBoundaryFlag(StoB(token[1]));
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qnOfSecondaryCmd) {
          if ( CheckMeshPS(mesh,token[0],command) ){
              G4PSNofSecondary3D* ps =new G4PSNofSecondary3D(token[0]);
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qTrackLengthCmd) {
          if ( CheckMeshPS(mesh,token[0],command) ){
              G4PSTrackLength3D* ps = new G4PSTrackLength3D(token[0]);
              ps->Weighted(StoB(token[1]));
              ps->MultiplyKineticEnergy(StoB(token[2]));
              ps->DivideByVelocity(StoB(token[3]));
              ps->SetUnit(token[4]);
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qPassCellCurrCmd){
          if( CheckMeshPS(mesh,token[0],command) ) {
              G4PSPassageCellCurrent* ps = new G4PSPassageCellCurrent3D(token[0]);
              ps->Weighted(StoB(token[1]));
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qPassTrackLengthCmd){
          if( CheckMeshPS(mesh,token[0],command) ) {
              G4PSPassageTrackLength* ps = new G4PSPassageTrackLength3D(token[0]);
              ps->Weighted(StoB(token[1]));
              ps->SetUnit(token[2]);
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qFlatSurfCurrCmd){
          if( CheckMeshPS(mesh,token[0],command)) {
              G4PSFlatSurfaceCurrent3D* ps = 
                new G4PSFlatSurfaceCurrent3D(token[0],StoI(token[1]));
              ps->Weighted(StoB(token[2]));
              ps->DivideByArea(StoB(token[3]));
              if ( StoB(token[3]) ){
                ps->SetUnit(token[4]);
              }else{
                ps->SetUnit("");
              }
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qFlatSurfFluxCmd){
          if( CheckMeshPS(mesh, token[0],command )) {
              G4PSFlatSurfaceFlux3D* ps = new G4PSFlatSurfaceFlux3D(token[0],StoI(token[1]));
              ps->Weighted(StoB(token[2]));
              ps->DivideByArea(StoB(token[3]));
              if ( StoB(token[3]) ){
                ps->SetUnit(token[4]);
              }else{
                ps->SetUnit("");
              }
              mesh->SetPrimitiveScorer(ps);
          }
//    } else if(command== qSphereSurfCurrCmd){
//        if( CheckMeshPS(mesh, token[0],command )) {
//            G4PSSphereSurfaceCurrent3D* ps = 
//              new G4PSSphereSurfaceCurrent3D(token[0],StoI(token[1]));
//            ps->Weighted(StoB(token[2]));
//            ps->DivideByArea(StoB(token[3]));
//            if ( StoB(token[3]) ){
//              ps->SetUnit(token[4]);
//            }else{
//              ps->SetUnit("");
//            }
//            mesh->SetPrimitiveScorer(ps);
//        }
//   } else if(command== qSphereSurfFluxCmd){
//        if( CheckMeshPS(mesh,token[0],command)) {
//            G4PSSphereSurfaceFlux3D* ps = new G4PSSphereSurfaceFlux3D(token[0], StoI(token[1]));
//            ps->Weighted(StoB(token[2]));
//            ps->DivideByArea(StoB(token[3]));
//            if ( StoB(token[3]) ){
//              ps->SetUnit(token[4]);
//            }else{
//              ps->SetUnit("");
//            }
//            mesh->SetPrimitiveScorer(ps);
//        }
//   } else if(command== qCylSurfCurrCmd){
//        if( CheckMeshPS(mesh, token[0],command ) ) {
//            G4PSCylinderSurfaceCurrent3D* ps = 
//              new G4PSCylinderSurfaceCurrent3D(token[0],StoI(token[1]));
//            ps->Weighted(StoB(token[2]));
//            ps->DivideByArea(StoB(token[3]));
//            if ( StoB(token[3]) ){
//              ps->SetUnit(token[4]);
//            }else{
//              ps->SetUnit("");
//            }
//            ps->SetUnit(token[4]);
//            mesh->SetPrimitiveScorer(ps);
//        }
//   } else if(command== qCylSurfFluxCmd){
//        if( CheckMeshPS(mesh, token[0],command ) {
//            G4PSCylinerSurfaceFlux3D* ps =new G4PSCylinderSurfaceFlux3D(token[0], StoI(token[1]));
//            ps->Weighted(StoB(token[2]));
//            ps->DivideByArea(StoB(token[3]));
//            if ( StoB(token[3]) ){
//              ps->SetUnit(token[4]);
//            }else{
//              ps->SetUnit("");
//            }
//            mesh->SetPrimitiveScorer(ps);
//        }
      } else if(command== qNofCollisionCmd){
          if( CheckMeshPS(mesh,token[0],command)) {
              G4PSNofCollision3D* ps =new G4PSNofCollision3D(token[0]); 
              ps->Weighted(StoB(token[1]));
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qPopulationCmd){
          if( CheckMeshPS(mesh,token[0],command) ) {
              G4PSPopulation3D* ps =new G4PSPopulation3D(token[0]); 
              ps->Weighted(StoB(token[1]));
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qTrackCountCmd){
          if( CheckMeshPS(mesh,token[0],command)) {
              G4PSTrackCounter3D* ps =new G4PSTrackCounter3D(token[0],StoI(token[1])); 
              ps->Weighted(StoB(token[2]));
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qTerminationCmd){
          if( CheckMeshPS(mesh,token[0],command)) {
              G4PSTermination3D* ps =new G4PSTermination3D(token[0]); 
              ps->Weighted(StoB(token[1]));
              mesh->SetPrimitiveScorer(ps);
          }

      } else if(command== qMinKinEAtGeneCmd){
          if( CheckMeshPS(mesh,token[0],command) ){
              G4PSMinKinEAtGeneration3D* ps =new G4PSMinKinEAtGeneration3D(token[0]); 
              ps->SetUnit(token[1]);
              mesh->SetPrimitiveScorer(ps);
          }
      } else if(command== qStepCheckerCmd){
          if( CheckMeshPS(mesh,token[0],command) ){
              G4PSStepChecker3D* ps =new G4PSStepChecker3D(token[0]); 
              mesh->SetPrimitiveScorer(ps);
          }

            //
            // Filters 
            // 
      }else if(command== fchargedCmd){
            if(!mesh->IsCurrentPrimitiveScorerNull()) {
              mesh->SetFilter(new G4SDChargedFilter(token[0])); 
            } else {
              ed << "WARNING[" << fchargedCmd->GetCommandPath()
                 << "] : Current quantity is not set. Set or touch a quantity first.";
              command->CommandFailed(ed);
            }
      }else if(command== fneutralCmd){
            if(!mesh->IsCurrentPrimitiveScorerNull()) {
              mesh->SetFilter(new G4SDNeutralFilter(token[0])); 
            } else {
              ed << "WARNING[" << fneutralCmd->GetCommandPath()
                 << "] : Current quantity is not set. Set or touch a quantity first.";
              command->CommandFailed(ed);
            }
      }else if(command== fkinECmd){
            if(!mesh->IsCurrentPrimitiveScorerNull()) {
              G4String& name = token[0];
              G4double elow  = StoD(token[1]);
              G4double ehigh = StoD(token[2]);
              G4double unitVal = G4UnitDefinition::GetValueOf(token[3]);
              mesh->SetFilter(new G4SDKineticEnergyFilter(name,elow*unitVal,ehigh*unitVal));
            } else {
              ed << "WARNING[" << fkinECmd->GetCommandPath()
                 << "] : Current quantity is not set. Set or touch a quantity first.";
              command->CommandFailed(ed);
            }
      }else if(command== fparticleKinECmd){
            if(!mesh->IsCurrentPrimitiveScorerNull()) {
              FParticleWithEnergyCommand(mesh,token); 
            } else {
              ed << "WARNING[" << fparticleKinECmd->GetCommandPath()
                 << "] : Current quantity is not set. Set or touch a quantity first.";
              command->CommandFailed(ed);
            }
      } else if(command==fparticleCmd) {
            if(!mesh->IsCurrentPrimitiveScorerNull()) {
              FParticleCommand(mesh,token);
            } else {
              ed << "WARNING[" << fparticleCmd->GetCommandPath()
                 << "] : Current quantity is not set. Set or touch a quantity first.";
              command->CommandFailed(ed);
            }
      }
}

G4String G4ScoreQuantityMessenger::GetCurrentValue(G4UIcommand * /*command*/)
{
  G4String val;

  return val;
}

void G4ScoreQuantityMessenger::FillTokenVec(G4String newValues, G4TokenVec& token){

    G4Tokenizer next(newValues);
    G4String val;
    while ( !(val = next()).isNull() ) { // Loop checking 12.18.2015 M.Asai
        token.push_back(val);
    }
}


void G4ScoreQuantityMessenger::FParticleCommand(G4VScoringMesh* mesh, G4TokenVec& token){
    //
    // Filter name
    G4String name = token[0];
    //
    // particle list
    std::vector<G4String> pnames;
    for ( G4int i = 1; i<(G4int)token.size(); i++){
        pnames.push_back(token[i]);
    }
    //
    // Attach Filter
    mesh->SetFilter(new G4SDParticleFilter(name,pnames));
}    

void G4ScoreQuantityMessenger::FParticleWithEnergyCommand(G4VScoringMesh* mesh,G4TokenVec& token){
    G4String& name = token[0];
    G4double  elow = StoD(token[1]);
    G4double  ehigh= StoD(token[2]);
    G4double  unitVal = G4UnitDefinition::GetValueOf(token[3]);
    G4SDParticleWithEnergyFilter* filter = 
        new G4SDParticleWithEnergyFilter(name,elow*unitVal,ehigh*unitVal);
    for ( G4int i = 4; i < (G4int)token.size(); i++){
        filter->add(token[i]);
    }
    mesh->SetFilter(filter);
}

G4bool G4ScoreQuantityMessenger::CheckMeshPS(G4VScoringMesh* mesh, G4String& psname,
                                             G4UIcommand* command){
    if(!mesh->FindPrimitiveScorer(psname)) {
        return true;
    } else {
        G4ExceptionDescription ed;
        ed << "WARNING[" << qTouchCmd->GetCommandPath()
           << "] : Quantity name, \"" << psname << "\", is already existing.";
        command->CommandFailed(ed);
        mesh->SetNullToCurrentPrimitiveScorer();
        return false;
    }
}
