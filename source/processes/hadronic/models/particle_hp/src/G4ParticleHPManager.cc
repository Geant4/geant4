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
// Class Description
// Manager of NetronHP
//
// 121031 First implementation done by T. Koi (SLAC/PPA)
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko, July-2023 Basic revision of particle HP classes
//
#include "G4ParticleHPManager.hh"

#include "G4Exception.hh"
#include "G4HadronicException.hh"
#include "G4ParticleHPMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicParameters.hh"
#include "G4ParticleHPThreadLocalManager.hh"
#include "G4SystemOfUnits.hh"

#include <zlib.h>
#include <fstream>

G4ParticleHPManager* G4ParticleHPManager::instance = nullptr;

G4ParticleHPManager::G4ParticleHPManager()
  : theMinEnergyDBRC(0.1 * CLHEP::eV),
    theMaxEnergyDBRC(210. * CLHEP::eV),
    theMaxEnergyDoppler(30. * CLHEP::keV)
{
  messenger = new G4ParticleHPMessenger(this);
  verboseLevel = G4HadronicParameters::Instance()->GetVerboseLevel();
  char* ss = std::getenv("NeutronHPNames");
  if (nullptr != ss) { CHECK_HP_NAMES = true; }
  ss = std::getenv("G4PHP_DO_NOT_CHECK_DIFF_COEFF_REPR");
  if (nullptr != ss) { PHP_CHECK = false; }
  ss = std::getenv("G4PHP_MULTIPLICITY_METHOD");
  if (nullptr != ss && "BetweenInts" == G4String(ss)) { PHP_USE_POISSON = false; }
  ss = std::getenv("G4ParticleHPDebug");
  if (nullptr != ss) { DEBUG = true; }

  // identify and check data path once - it should exist
  const char* nch = G4FindDataDir("G4NEUTRONHPDATA");
  if (nullptr == nch) {
    G4Exception("G4ParticleHPManager::G4ParticleHPManager()","hadhp01",
                FatalException, "G4NEUTRONXSDATA is not defined - check path");
  } else {
    fDataPath[0] = G4String(nch);
  }
  // path may be defined by two environment variables
  // it is not mandatory to access PHP data - path may be not defined
  const char* ttp = G4FindDataDir("G4PARTICLEHPDATA");
  G4String tendl = (nullptr == ttp) ? "" : G4String(ttp);
  const char* ssp = G4FindDataDir("G4PROTONHPDATA");
  fDataPath[1] = (nullptr == ssp) ? tendl + "/Proton" : G4String(ssp);

  ssp = G4FindDataDir("G4DEUTERONHPDATA");
  fDataPath[2] = (nullptr == ssp) ? tendl + "/Deuteron" : G4String(ssp);

  ssp = G4FindDataDir("G4TRITONHPDATA");
  fDataPath[3] = (nullptr == ssp) ? tendl + "/Triton" : G4String(ssp);

  ssp = G4FindDataDir("G4HE3HPDATA");
  fDataPath[4] = (nullptr == ssp) ? tendl + "/He3" : G4String(ssp);

  ssp = G4FindDataDir("G4ALPHAHPDATA");
  fDataPath[5] = (nullptr == ssp) ? tendl + "/Alpha" : G4String(ssp);
}

G4ParticleHPManager::~G4ParticleHPManager()
{
  delete messenger;
}

G4ParticleHPManager* G4ParticleHPManager::GetInstance()
{
  static G4ParticleHPManager manager;
  if (instance == nullptr) {
    instance = &manager;
  }
  return instance;
}

G4int G4ParticleHPManager::GetPHPIndex(const G4ParticleDefinition* part) const {
  G4int pdg = part->GetPDGEncoding();
  G4int idx;
  if (pdg == 2112) { idx = 0; }
  else if (pdg == 2212) { idx = 1; }
  else if (pdg == 1000010020) { idx = 2; }
  else if (pdg == 1000010030) { idx = 3; }
  else if (pdg == 1000020030) { idx = 4; }
  else if (pdg == 1000020040) { idx = 5; }
  else {
    idx = 0;
    G4ExceptionDescription ed;
    ed << "Particle " << part->GetParticleName()
       << " cannot be handled by the ParticleHP sub-library";
    G4Exception("G4ParticleHPManager::G4ParticleHPManager()","hadhp01",
                FatalException, ed, "");
  }
  return idx;
}

const G4String&
G4ParticleHPManager::GetParticleHPPath(const G4ParticleDefinition* part) const {
  return fDataPath[GetPHPIndex(part)];
}

void G4ParticleHPManager::OpenReactionWhiteBoard()
{
  G4ParticleHPThreadLocalManager::GetInstance()->OpenReactionWhiteBoard();
}

G4ParticleHPReactionWhiteBoard* G4ParticleHPManager::GetReactionWhiteBoard()
{
  return G4ParticleHPThreadLocalManager::GetInstance()->GetReactionWhiteBoard();
}

void G4ParticleHPManager::CloseReactionWhiteBoard()
{
  G4ParticleHPThreadLocalManager::GetInstance()->CloseReactionWhiteBoard();
}

void G4ParticleHPManager::GetDataStream(const G4String& filename, std::istringstream& iss)
{
  G4String* data = nullptr;
  G4String compfilename(filename);
  compfilename += ".z";
  auto in = new std::ifstream(compfilename, std::ios::binary | std::ios::ate);
  if (in->good()) {
    // Use the compressed file
    std::streamoff file_size = in->tellg();
    in->seekg(0, std::ios::beg);
    auto compdata = new Bytef[file_size];

    while (*in) {  // Loop checking, 11.05.2015, T. Koi
      in->read((char*)compdata, file_size);
    }

    auto complen = (uLongf)(file_size * 4);
    auto uncompdata = new Bytef[complen];

    while (Z_OK != uncompress(uncompdata, &complen, compdata, file_size))
    {  // Loop checking, 11.05.2015, T. Koi
      delete[] uncompdata;
      complen *= 2;
      uncompdata = new Bytef[complen];
    }
    delete[] compdata;
    // Now "complen" has uncomplessed size
    data = new G4String((char*)uncompdata, (G4long)complen);
    delete[] uncompdata;
  }
  else {
    // Use regular text file
    std::ifstream thefData(filename, std::ios::in | std::ios::ate);
    if (thefData.good()) {
      std::streamoff file_size = thefData.tellg();
      thefData.seekg(0, std::ios::beg);
      auto filedata = new char[file_size];
      while (thefData) {  // Loop checking, 11.05.2015, T. Koi
        thefData.read(filedata, file_size);
      }
      thefData.close();
      data = new G4String(filedata, file_size);
      delete[] filedata;
    }
    else {
      // found no data file
      // set error bit to the stream
      iss.setstate(std::ios::badbit);
    }
  }
  if (data != nullptr) {
    iss.str(*data);
    G4String id;
    iss >> id;
    if (id == "G4NDL") {
      // Register information of file
      G4String source;
      iss >> source;
      register_data_file(filename, source);
    }
    else {
      iss.seekg(0, std::ios::beg);
    }
  }
  in->close();
  delete in;
  delete data;
}

void G4ParticleHPManager::GetDataStream2(const G4String& filename, std::istringstream& iss)
{
  // Checking existance of data file

  G4String compfilename(filename);
  compfilename += ".z";
  auto in = new std::ifstream(compfilename, std::ios::binary | std::ios::ate);
  if (in->good()) {
    // Compressed file is exist
    in->close();
  }
  else {
    std::ifstream thefData(filename, std::ios::in | std::ios::ate);
    if (thefData.good()) {
      // Regular text file is exist
      thefData.close();
    }
    else {
      // found no data file
      // set error bit to the stream
      iss.setstate(std::ios::badbit);
    }
  }
  delete in;
}

void G4ParticleHPManager::SetVerboseLevel(G4int newValue)
{
  G4cout << "You are setting a new verbose level for Particle HP package." << G4endl;
  G4cout << "the new value will be used in whole of the Particle HP package, i.e., models and "
            "cross sections for Capture, Elastic, Fission and Inelastic interaction."
         << G4endl;
  verboseLevel = newValue;
}

void G4ParticleHPManager::register_data_file(const G4String& filename, const G4String& source)
{
  mDataEvaluation.insert(std::pair<G4String, G4String>(filename, source));
}

void G4ParticleHPManager::DumpDataSource()
{
  G4cout << "Data source of this Partile HP calculation are " << G4endl;
  for (const auto& it : mDataEvaluation) {
    G4cout << it.first << " " << it.second << G4endl;
  }
  G4cout << G4endl;
}

void G4ParticleHPManager::DumpSetting()
{
  if(isPrinted) { return; }
  G4cout << G4endl
         << "=======================================================" << G4endl
         << "======       ParticleHP Physics Parameters     ========" << G4endl
         << "=======================================================" << G4endl
         << " Use only photo-evaporation      " << USE_ONLY_PHOTONEVAPORATION << G4endl
         << " Skip missing isotopes           " << SKIP_MISSING_ISOTOPES << G4endl
         << " Neglect Doppler                 " << NEGLECT_DOPPLER << G4endl
         << " Do not adjust final state       " << DO_NOT_ADJUST_FINAL_STATE << G4endl
         << " Produce fission fragments       " << PRODUCE_FISSION_FRAGMENTS << G4endl
         << " Use WendtFissionModel           " << USE_WENDT_FISSION_MODEL << G4endl
         << " Use NRESP71Model                " << USE_NRESP71_MODEL << G4endl
         << " Use DBRC                        " << USE_DBRC << G4endl
         << " PHP use Poisson                 " << PHP_USE_POISSON << G4endl
         << " PHP check                       " << PHP_CHECK << G4endl
         << " CHECK HP NAMES                  " << CHECK_HP_NAMES << G4endl
         << " Enable DEBUG                    " << DEBUG << G4endl
         << "=======================================================" << G4endl << G4endl;
  isPrinted = true;
}
