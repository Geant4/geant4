#include "G4DosimPhysicsMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4DosimPhysics.hh"

G4DosimPhysicsMessenger::G4DosimPhysicsMessenger(PhysicsList* pMother): m_pMother(pMother)
{
  m_pChangeModel = new G4UIcmdWithAString("/process/change",this);
  m_pChangeModel->SetGuidance("Switch between low energy inelastic model (LEinelastic),");
  m_pChangeModel->SetGuidance("        precompound model (Precompound) and");
  m_pChangeModel->SetGuidance("        Mars below 5 GeV (Mars)");
  m_pChangeModel->SetParameterName("ModelName",true);
  m_pChangeModel->SetDefaultValue("LEinelastic");
  m_pChangeModel->AvailableForStates(Idle);
}

G4DosimPhysicsMessenger::~G4DosimPhysicsMessenger()
{
  delete m_pChangeModel;
}

void G4DosimPhysicsMessenger::SetNewValue(G4UIcommand* pCmd,G4String szValue)
{
  if(pCmd==m_pChangeModel){
    if(szValue=="LEinelastic")
      m_pMother->ChangeModel(LEinelastic);
    else if(szValue=="Precompound")
      m_pMother->ChangeModel(Precompound);
    else if(szValue=="Mars")
      m_pMother->ChangeModel(Mars);
    else if(szValue == "Chiral")
      m_pMother->ChangeModel(Chiral);
    else if(szValue=="new")
      m_pMother->ChangeModel(PrecompoundNew);
    else if(szValue=="all")
      m_pMother->ChangeModel(PrecompoundAll);
    else
      G4cout<<szValue<<" unknown model"<<G4endl;
  }
}
