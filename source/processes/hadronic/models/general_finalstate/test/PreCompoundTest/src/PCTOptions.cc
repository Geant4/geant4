#include "PCTOptions.hh"

#include "PCTTools.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include <locale>

void PCTOptions::AskForFermiBreakUpMode()
{
  this->PutEmptyLine(2);
  // Use Fermi (y/n)
  std::cout << "Fermi BreakUp Model\n"
	      << "-------------------\n";
  
  char fermi_yn = '*';
  do 
    {
      std::cout << "Do you want to use FermiBreakUp Model? (y/n): ";
      std::cin >> fermi_yn;
      fermi_yn = toupper(fermi_yn);
    }
  while (fermi_yn != 'Y' && fermi_yn != 'N');

  if (fermi_yn == 'Y') fermi = true;
  else fermi = false;

  return;
}

void PCTOptions::AskForEvaporationMode()
{
  this->PutEmptyLine(2);
  // EVAP MODE (standard or GEM)
  do
    {
      G4String aMode;
      std::cout << "Use Standard (std) Evaporation or GEM model: ";
      std::cin >> aMode;
      std::transform(aMode.begin(),aMode.end(),
                     aMode.begin(),
                     toupper);
      if (aMode == "STANDARD" || aMode == "STD") theEvapMode = standard_mode;
      else if (aMode == "GEM") theEvapMode = GEM_mode;
      else theEvapMode = no_mode;
    }
  while (theEvapMode != standard_mode && theEvapMode != GEM_mode);
  return;
}

void PCTOptions::AskForPreeqEmissionMode()
{
  this->PutEmptyLine(2);
  // PRECOMPOUND EMISSION MODE (standard or HETC)
  do
    {
      G4String aMode;
      std::cout << "Use Standard (std) or HETC precompound emission model: ";
      std::cin >> aMode;
      std::transform(aMode.begin(),aMode.end(),
                     aMode.begin(),
                     toupper);
      if (aMode == "STANDARD" || aMode == "STD") thePreeqEmissionMode = standard_emission_mode;
      else if (aMode == "HETC") thePreeqEmissionMode = HETC_emission_mode;
      else thePreeqEmissionMode = no_emission_mode;
    }
  while (thePreeqEmissionMode != standard_emission_mode && 
	 thePreeqEmissionMode != HETC_emission_mode);
  return;
}

void PCTOptions::AskForPreeqTransitionMode()
{
  this->PutEmptyLine(2);
  // PRECOMPOUND EMISSION MODE (standard or GNASH)
  do
    {
      G4String aMode;
      std::cout << "Use Standard (std) or GNASH precompound transition model: ";
      std::cin >> aMode;
      std::transform(aMode.begin(),aMode.end(),
                     aMode.begin(),
                     toupper);
      if (aMode == "STANDARD" || aMode == "STD") thePreeqTransitionMode = standard_transition_mode;
      else if (aMode == "GNASH") thePreeqTransitionMode = GNASH_transition_mode;
      else thePreeqTransitionMode = no_transition_mode;
    }
  while (thePreeqTransitionMode != standard_transition_mode && 
	 thePreeqTransitionMode != GNASH_transition_mode);
  return;
}


void PCTOptions::AskForTarget()
{
  this->PutEmptyLine(2);
  std::cout << "Target Nucleus\n"
	      << "--------------\n";
  do 
    {
      G4String z;
      std::cout << "Z = ";
      std::cin >> z;
      tZ = StrToInt(z);
    }
  while (tZ < 0);
  do 
    {
      G4String a;
      std::cout << "A (0 for natural elemet) = ";
      std::cin >> a;
      tA = StrToInt(a);
    }
  while (tA < 0 || (tA != 0 && tA < tZ) );
  if (tA > 0) natural = false;
  else 
    {
      natural = true;
      G4int isotopeA = 0;
      G4double porcentaje = 0.0;
      G4double norma = 0.0;
      G4int numisotopes = 0;
      G4String snisotopes;
      do 
	{
	  std::cout << "\tNumber of isotopes: ";
	  std::cin >> snisotopes;
	  numisotopes = StrToInt(snisotopes);
	}
      while (numisotopes < 1);
      material.reserve(numisotopes);
      for (G4int isotop = 0; isotop < numisotopes; isotop++)
	{
	  std::cout << "\tIsotope n. " << isotop + 1 << '\n';
	  std::cout << "\t\tA = ";
	  std::cin >> isotopeA;
	  std::cout << "\t\tabundance = ";
	  std::cin >> porcentaje;
	  norma += porcentaje;
	  material.push_back(std::make_pair(isotopeA,norma));
	}
      if (norma != 1.0 && norma != 100.0)
	{
	  std::cout << "Composition is not normalized to the unit. "
		      << "I can manage this properly, but you should\n"
		      << "check that the data you have input is correct.\n";
	}
      G4double prev = 0.0;
      for (std::vector< std::pair<G4int,G4double> >::iterator it = material.begin();
	   it != material.end(); it++)
	{
	  std::cout << "A = " << it->first << "\t(" << it->second - prev << ")\n";
	  prev = it->second;
	}
    }
  return;
}


void PCTOptions::AskForProjectile()
{
  this->PutEmptyLine(2);
  std::cout << "Projetile\n"
	      << "---------\n";
  pA = 1;
  pZ = 1;
  G4String T;
  do 
    {
      std::cout << "Projectile Kinetic Energy (in MeV) : ";
      std::cin >> T;
      pKineticEnergy = StrToDouble(T)*MeV;
    }
  while (pKineticEnergy < 0.0);
  char dirtype = '*';
  do
    {
      std::cout << "Shoot projetile in (X), (Y), (Z) axe or (R)andom direction: ";
      std::cin  >> dirtype; 
      dirtype = toupper(dirtype);
    }
  while (dirtype != 'X' && dirtype != 'Y' && dirtype != 'Z' && dirtype != 'R');
  if (dirtype == 'X')
    pdirection = XaxisDirection;
  else if (dirtype == 'Y')
    pdirection = YaxisDirection;
  else if (dirtype == 'Z')
    pdirection = ZaxisDirection;
  else 
    pdirection = RandomDirection;
  return;
}


void PCTOptions::AskForINC()
{
  this->PutEmptyLine(2);
  std::cout << "Intra nuclear cascade\n"
	      << "---------------------\n";
  
  char inc_yn = '*';
  do 
    {
      std::cout << "Do you want to use INC? (y/n): ";
      std::cin >> inc_yn;
      inc_yn = toupper(inc_yn);
    }
  while (inc_yn != 'Y' && inc_yn != 'N');

  if (inc_yn == 'Y') inc = true;
  else inc = false;

  return;
}

void PCTOptions::AskForExcitons()
{
  this->PutEmptyLine(2);
  std::cout << "Excitons\n"
	      << "--------\n";
  particles = 2;      
  holes = 1;
  charged = 1;

  char autoexciton = '*';
  do 
    {
      std::cout << "Do you want automatic exciton configuration? (y/n): ";
      std::cin >> autoexciton;
      autoexciton = toupper(autoexciton);
    }
  while (autoexciton != 'Y' && autoexciton != 'N');

  //    if (G4UniformRand() > G4double(Z0)/G4double(A0)) MyCharged++;

  if (autoexciton == 'N')
    {
      G4String tmp;
      do 
	{
	  std::cout << "Number of Particles: ";
	  std::cin >> tmp;
	  particles = StrToInt(tmp);
	}
      while (particles < 0);
      do 
	{
	  std::cout << "Number of Holes: ";
	  std::cin >> tmp;
	  holes = StrToInt(tmp);
	}
      while (holes < 0);
      do 
	{
	  std::cout << "Number of Charged Particles ( < " << particles << " ): ";
	  std::cin >> tmp;
	  charged = StrToInt(tmp);
	}
      while (charged > particles);
    }
  else 
    {
      std::cout << "The initial configuration is:\n\tparticles =========> " << particles
                <<                              "\n\tholes =============> " << holes
                <<                              "\n\tcharged particles => " << charged
                << '\n';
    }
  return;
}

void PCTOptions::AskForIterations()
{
  this->PutEmptyLine(2);
  iterations = -1;
  G4String tmp;
  do 
    {
      std::cout << "Number of Iterations: ";
      std::cin >> tmp;
      iterations = StrToInt(tmp);
    }
  while (iterations < 0);
  return;
}

void PCTOptions::GenerateFilename()
{
  
  G4int anA = tA;
  if (natural) anA = material.begin()->first;
  filename = (const_cast<G4IonTable*>(G4ParticleTable::GetParticleTable()->GetIonTable()))->
    GetIon(int(tZ),int(anA))->GetParticleName();
  G4String::size_type idx = filename.find('[');
  if (idx != G4String::npos) 
    {
      filename = filename.substr(0,idx);
    }
  if ( natural ) // remove the isotope number and put a 0
    {
      idx = filename.find_first_of("0123456789");
      if (idx != G4String::npos)
	{
	  filename = filename.substr(0,idx);
	  filename += '0';
	}
    }
  filename += "-";
  filename += DoubleToStr(pKineticEnergy/MeV);
  filename += "MeV-";

  if (this->UsingINC())
    filename += "INC-";

  if (this->UsingFermi())
    filename += "FERMI-";

  if (thePreeqEmissionMode == HETC_emission_mode)
    filename += "HETC";
  else 
    filename += "STD";

  if (thePreeqTransitionMode == GNASH_transition_mode)
    filename += "+GNASH";
  else 
    filename += "+STD";

  if (theEvapMode == GEM_mode) 
    filename += "+GEM";
  else 
    filename += "+STD";
  
  if (pdirection == XaxisDirection)
    filename += "-X";
  else if (pdirection == YaxisDirection)
    filename += "-Y";
  else if (pdirection == ZaxisDirection)
    filename += "-Z";
  else 
    filename += "-R";
  
  // Add the initial exciton configuration
  filename += "-" + IntToStr(particles) + "p" + 
    IntToStr(holes) + "h" + 
    IntToStr(charged) + "c";

  // Add extension
  filename += ".dat.gz";

  return;
}


void PCTOptions::Initialize()
{
  this->AskForFermiBreakUpMode();
  this->AskForPreeqEmissionMode();
  this->AskForPreeqTransitionMode();
  this->AskForEvaporationMode();
  this->AskForTarget();
  this->AskForProjectile();
  this->AskForINC();
  if ( !this->UsingINC() )
    this->AskForExcitons();
  this->AskForIterations();
  this->GenerateFilename();
  return;
}

