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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLParticleType.hh"
#include "G4INCLConfig.hh"
#include "G4INCLParticleSpecies.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  const G4int Config::randomSeedMin = 1;
  const G4int Config::randomSeedMax = ((1<<30)-1)+(1<<30); // 2^31-1

  Config::Config()
  {
    init();
  }

  Config::Config(G4int /*A*/, G4int /*Z*/, G4INCL::ParticleSpecies proj, G4double projectileE)
  {
    init();
    projectileSpecies = proj;
    projectileKineticEnergy = projectileE;
  }

  // NOT used in Geant4 mode
#if defined(HAS_BOOST_PROGRAM_OPTIONS) && !defined(INCLXX_IN_GEANT4_MODE)
  Config::Config(G4int argc, char *argv[], G4bool isFullRun) :
    runOptDesc("Run options"),
    hiddenOptDesc("Hidden options"),
    genericOptDesc("Generic options"),
    physicsOptDesc("Physics options"),
    naturalTarget(false)
  {
    const std::string suggestHelpMsg("You might want to run `INCLCascade --help' to get a help message.\n");

    // Define the names of the de-excitation models
    const std::string theNoneName = "none";
#ifdef INCL_DEEXCITATION_SMM
    const std::string theSMMName = "SMM";
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    const std::string theGEMINIXXName = "GEMINIXX";
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
    const std::string theABLAv3pName = "ABLAv3p";
#endif
#ifdef INCL_DEEXCITATION_ABLA07
    const std::string theABLA07Name = "ABLA07";
#endif

    // Define the default de-excitation model, in decreasing order of priority
    std::string defaultDeExcitationModel = theNoneName;
#ifdef INCL_DEEXCITATION_SMM
    defaultDeExcitationModel = theSMMName;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    defaultDeExcitationModel = theGEMINIXXName;
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
    defaultDeExcitationModel = theABLAv3pName;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
    defaultDeExcitationModel = theABLA07Name;
#endif

    const std::string listSeparator = "\n  \t";
    deExcitationModelList =
      listSeparator + theNoneName
#ifdef INCL_DEEXCITATION_ABLA07
      + listSeparator + theABLA07Name
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
      + listSeparator + theABLAv3pName
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
      + listSeparator + theGEMINIXXName
#endif
#ifdef INCL_DEEXCITATION_SMM
      + listSeparator + theSMMName
#endif
      ;

    // Append " (default)" to the name of the default model
    size_t defaultModelIndex = deExcitationModelList.find(defaultDeExcitationModel);
    if(defaultModelIndex!=std::string::npos) {
      deExcitationModelList = deExcitationModelList.substr(0, defaultModelIndex+defaultDeExcitationModel.size())
        + " (default)"
        + deExcitationModelList.substr(defaultModelIndex+defaultDeExcitationModel.size(), std::string::npos);
    }

    // Spell out the G4bool values
    std::cout << std::boolalpha;

    try {

      // Hidden options
      hiddenOptDesc.add_options()
        ("input-file", po::value<std::string>(&inputFileName), "input file")
        ("impact-parameter", po::value<G4double>(&impactParameter)->default_value(-1.), "impact parameter")
        ;

      // Generic options
      std::stringstream verbosityDescription;
      verbosityDescription << "set verbosity level:\n"
        << " 0: \tquiet, suppress all output messages\n"
        << " " << InfoMsg << ": \tminimal logging\n"
        << " " << FatalMsg << ": \tlog fatal error messages as well\n"
        << " " << ErrorMsg << ": \tlog error messages as well\n"
        << " " << WarningMsg << ": \tlog warning messages as well\n"
        << " " << DebugMsg << ": \tlog debug messages as well\n"
        << " " << DataBlockMsg << ": \tlog data-block messages as well";

      genericOptDesc.add_options()
        ("help,h", "produce this help message")
        ("version", "print version string and exit")
        ;

      // Run-specific options
      std::stringstream randomSeed1Description, randomSeed2Description;
      randomSeed1Description << "first seed for the random-number generator (between "
        << randomSeedMin << "and " << randomSeedMax << ")";
      randomSeed2Description << "second seed for the random-number generator (between "
        << randomSeedMin << "and " << randomSeedMax << ")";

      runOptDesc.add_options()
        ("title", po::value<std::string>(&title)->default_value("INCL default run title"), "run title")
        ("output,o", po::value<std::string>(&outputFileRoot), "root for generating output file names. File-specific suffixes (.root, .out, etc.) will be appended to this root. Defaults to the input file name, if given; otherwise, defaults to a string composed of the explicitly specified options and of a customisable suffix, if provided using the -s option")
        ("suffix,s", po::value<std::string>(&fileSuffix), "suffix to be appended to generated output file names")
        ("logfile,l", po::value<std::string>(&logFileName), "log file name. Defaults to `<output_root>.log'. Use `-' if you want to redirect logging to stdout")
        ("number-shots,N", po::value<G4int>(&nShots), "* number of shots")
        ("target,t", po::value<std::string>(&targetString), "* target nuclide. Can be specified as Fe56, 56Fe, Fe-56, 56-Fe, Fe_56, 56_Fe or Fe. If the mass number is omitted, natural target composition is assumed.")
        ("projectile,p", po::value<std::string>(&projectileString), "* projectile name:\n"
         "  \tproton, p\n"
         "  \tneutron, n\n"
         "  \tpi+, piplus, pion+, pionplus\n"
         "  \tpi0, pizero, pion0, pionzero\n"
         "  \tpi-, piminus, pion-, pionminus\n"
         "  \td, t, a, deuteron, triton, alpha\n"
         "  \tHe-4, He4, 4He (and so on)")
        ("energy,E", po::value<G4double>(&projectileKineticEnergy), "* total kinetic energy of the projectile, in MeV")
        ("verbose-event", po::value<G4int>(&verboseEvent)->default_value(-1), "request verbose logging for the specified event only")
        ("random-seed-1", po::value<G4int>(&randomSeed1)->default_value(666), randomSeed1Description.str().c_str())
        ("random-seed-2", po::value<G4int>(&randomSeed2)->default_value(777), randomSeed2Description.str().c_str())
#ifdef INCL_ROOT_USE
        ("root-selection", po::value<std::string>(&rootSelectionString)->default_value(""), "ROOT selection for abridged output ROOT tree. For example: \"A==1 && Z==0 && theta<3\" selects only events where a neutron is scattered in the forward direction.")
#endif
        ("inclxx-datafile-path", po::value<std::string>(&INCLXXDataFilePath)->default_value("../data/"),
         "path to the INCL++ data files")
#ifdef INCL_DEEXCITATION_ABLA07
        ("abla07-datafile-path", po::value<std::string>(&abla07DataFilePath)->default_value("../de-excitation/abla07/upstream/tables/"),
         "path to the ABLA07 data files")
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
        ("ablav3p-cxx-datafile-path", po::value<std::string>(&ablav3pCxxDataFilePath)->default_value("../de-excitation/ablaxx/data/G4ABLA3.0/"),
         "path to the ABLAv3p data files")
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
        ("geminixx-datafile-path", po::value<std::string>(&geminixxDataFilePath)->default_value("../de-excitation/geminixx/upstream/"),
         "path to the GEMINI++ data files")
#endif
        ("verbosity,v", po::value<G4int>(&verbosity)->default_value(4), verbosityDescription.str().c_str())
        ;

      // Physics options
      physicsOptDesc.add_options()
        ("de-excitation,d", po::value<std::string>(&deExcitationString)->default_value(defaultDeExcitationModel.c_str()), ("which de-excitation model to use:" + deExcitationModelList).c_str())
#ifdef INCL_DEEXCITATION_FERMI_BREAKUP
        ("max-mass-fermi-breakup", po::value<G4int>(&maxMassFermiBreakUp)->default_value(18), "Maximum remnant mass for Fermi break-up. Default: 18.")
#endif
        ("pauli", po::value<std::string>(&pauliString)->default_value("strict-statistical"), "Pauli-blocking algorithm:\n"
         "  \tstrict-statistical (default)\n"
         "  \tstrict\n"
         "  \tstatistical\n"
         "  \tglobal\n"
         "  \tnone")
        ("cdpp", po::value<G4bool>(&CDPP)->default_value(true), "whether to apply CDPP after collisions:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("coulomb", po::value<std::string>(&coulombString)->default_value("non-relativistic"), "Coulomb-distortion algorithm:\n  \tnon-relativistic (default)\n  \tnone")
        ("potential", po::value<std::string>(&potentialString)->default_value("isospin-energy"), "nucleon potential:\n  \tisospin-energy-smooth\n  \tisospin-energy (default)\n  \tisospin\n  \tconstant")
        ("pion-potential", po::value<G4bool>(&pionPotential)->default_value("true"), "whether to use a pion potential:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("local-energy-BB", po::value<std::string>(&localEnergyBBString)->default_value("first-collision"), "local energy in baryon-baryon collisions:\n  \talways\n  \tfirst-collision (default)\n  \tnever")
        ("local-energy-pi", po::value<std::string>(&localEnergyPiString)->default_value("first-collision"), "local energy in pi-N collisions and in delta decays:\n  \talways\n  \tfirst-collision (default)\n  \tnever")
        ("cluster-algorithm", po::value<std::string>(&clusterAlgorithmString)->default_value("intercomparison"), "clustering algorithm for production of composites:\n  \tintercomparison (default)\n  \tnone")
        ("cluster-max-mass", po::value<G4int>(&clusterMaxMass)->default_value(8), "maximum mass of produced composites:\n  \tminimum 2\n  \tmaximum 12")
        ("back-to-spectator", po::value<G4bool>(&backToSpectator)->default_value("true"), "whether to use back-to-spectator:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("use-real-masses", po::value<G4bool>(&useRealMasses)->default_value("true"), "whether to use real masses for the outgoing particle energies:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("separation-energies", po::value<std::string>(&separationEnergyString)->default_value("INCL"), "how to assign the separation energies of the INCL nucleus:\n  \tINCL (default)\n  \treal\n  \treal-light")
        ("fermi-momentum", po::value<std::string>(&fermiMomentumString)->default_value("constant"), "how to assign the Fermi momentum of the INCL nucleus:\n  \tconstant (default)\n  \tconstant-light\n  \tmass-dependent")
        ("cutNN", po::value<G4double>(&cutNN)->default_value(1910.), "minimum CM energy for nucleon-nucleon collisions, in MeV. Default: 1910.")
        ("rp-correlation", po::value<G4double>(&rpCorrelationCoefficient)->default_value(1.), "correlation coefficient for the r-p correlation. Default: 1 (full correlation).")
        ("rp-correlation-p", po::value<G4double>(&rpCorrelationCoefficientProton)->default_value(1.), "correlation coefficient for the proton r-p correlation. Overrides the value specified using the rp-correlation option. Default: 1 (full correlation).")
        ("rp-correlation-n", po::value<G4double>(&rpCorrelationCoefficientNeutron)->default_value(1.), "correlation coefficient for the neutron r-p correlation. Overrides the value specified using the rp-correlation option. Default: 1 (full correlation).")
        ("neutron-skin-thickness", po::value<G4double>(&neutronSkinThickness)->default_value(0.), "thickness of the neutron skin, in fm. Default: 0.")
        ("neutron-skin-additional-diffuseness", po::value<G4double>(&neutronSkinAdditionalDiffuseness)->default_value(0.), "additional diffuseness of the neutron density distribution (with respect to the proton diffuseness), in fm. Default: 0.")
        ("refraction", po::value<G4bool>(&refraction)->default_value(false), "whether to use refraction when particles are transmitted. Default: false.")
        ;

      // Select options allowed on the command line
      po::options_description cmdLineOptions;
      cmdLineOptions.add(hiddenOptDesc).add(genericOptDesc).add(runOptDesc).add(physicsOptDesc);

      // Select options allowed in config files
      po::options_description configFileOptions;
      configFileOptions.add(runOptDesc).add(physicsOptDesc);

      // Select visible options
      po::options_description visibleOptions;
      visibleOptions.add(genericOptDesc).add(runOptDesc).add(physicsOptDesc);

      // Declare input-file as a positional option (if we just provide a file
      // name on the command line, it should be interpreted as an input-file
      // option).
      po::positional_options_description p;
      p.add("input-file", 1);

      // Disable guessing of option names
      G4int cmdstyle =
        po::command_line_style::default_style &
        ~po::command_line_style::allow_guessing;

      // Result of the option processing
      po::store(po::command_line_parser(argc, argv).
          style(cmdstyle).
          options(cmdLineOptions).positional(p).run(), variablesMap);
      po::notify(variablesMap);

      // If an input file was specified, merge the options with the command-line
      // options.
      if(variablesMap.count("input-file")) {
        std::ifstream inputFileStream(inputFileName.c_str());
        if(!inputFileStream) {
          std::cerr << "Cannot open input file: " << inputFileName << std::endl;
          std::exit(EXIT_FAILURE);
        } else {
          // Merge options from the input file
          po::parsed_options parsedOptions = po::parse_config_file(inputFileStream, configFileOptions, true);

          // Make sure that the unhandled options are all "*-datafile-path"
          std::vector<std::string> unhandledOptions =
            po::collect_unrecognized(parsedOptions.options, po::exclude_positional);
          G4bool ignoreNext = false;
          const std::string match = "-datafile-path";
          for(std::vector<std::string>::const_iterator i=unhandledOptions.begin(), e=unhandledOptions.end(); i!=e; ++i) {
            if(ignoreNext) {
              ignoreNext=false;
              continue;
            }
            if(i->rfind(match) == i->length()-match.length()) {
              std::cout << "Ignoring unrecognized option " << *i << std::endl;
              ignoreNext = true;
            } else {
              std::cerr << "Error: unrecognized option " << *i << std::endl;
              std::cerr << suggestHelpMsg;
              std::exit(EXIT_FAILURE);
            }
          }

          // Store the option values in the variablesMap
          po::store(parsedOptions, variablesMap);
          po::notify(variablesMap);
        }
        inputFileStream.close();
      }

      // Process the options from the user-specific config file ~/.inclxxrc
      std::string configFileName;
      const char * const configFileVar = getenv("INCLXXRC");
      if(configFileVar)
        configFileName = configFileVar;
      else {
        const char * const homeDirectoryPointer = getenv("HOME");
        if(homeDirectoryPointer) { // Check if we can find the home directory
          std::string homeDirectory(homeDirectoryPointer);
          configFileName = homeDirectory + "/.inclxxrc";
        } else {
          std::cerr << "Could not determine the user's home directory. "
            << "Are you running Linux, Unix or BSD?"<< std::endl;
          std::exit(EXIT_FAILURE);
        }
      }

      std::ifstream configFileStream(configFileName.c_str());
      std::cout << "Reading config file " << configFileName << std::endl;
      if(!configFileStream) {
        std::cout << "INCL++ config file " << configFileName
          << " not found. Continuing the run regardless."
          << std::endl;
      } else {
        // Merge options from the input file
        po::parsed_options parsedOptions = po::parse_config_file(configFileStream, configFileOptions, true);
        po::store(parsedOptions, variablesMap);

        // Make sure that the unhandled options are all "*-datafile-path"
        std::vector<std::string> unhandledOptions =
          po::collect_unrecognized(parsedOptions.options, po::exclude_positional);
        G4bool ignoreNext = false;
        const std::string match = "-datafile-path";
        for(std::vector<std::string>::const_iterator i=unhandledOptions.begin(), e=unhandledOptions.end(); i!=e; ++i) {
          if(ignoreNext) {
            ignoreNext=false;
            continue;
          }
          if(i->rfind(match) == i->length()-match.length()) {
            std::cout << "Ignoring unrecognized option " << *i << std::endl;
            ignoreNext = true;
          } else {
            std::cerr << "Error: unrecognized option " << *i << std::endl;
            std::cerr << suggestHelpMsg;
            std::exit(EXIT_FAILURE);
          }
        }

        // Store the option values in the variablesMap
        po::store(parsedOptions, variablesMap);
        po::notify(variablesMap);
      }
      configFileStream.close();

      /* *******************
       * Process the options
       * *******************/

      // -h/--help: print the help message and exit successfully
      if(variablesMap.count("help")) {
        std::cout
          << "Usage: INCLCascade [options] <input_file>" << std::endl
          << std::endl << "Options marked with a * are compulsory, i.e. they must be provided either on\nthe command line or in the input file." << std::endl
          << visibleOptions << std::endl;
        std::exit(EXIT_SUCCESS);
      }

      // --version: print the version string and exit successfully
      if(variablesMap.count("version")) {
        std::cout <<"INCL++ version " << getVersionString() << std::endl;
        std::exit(EXIT_SUCCESS);
      }

      // Check if the required options are present
      if(isFullRun) {
        std::string missingOption("");
        if(!variablesMap.count("number-shots"))
          missingOption = "number-shots";
        else if(!variablesMap.count("target"))
          missingOption = "target";
        else if(!variablesMap.count("projectile"))
          missingOption = "projectile";
        else if(!variablesMap.count("energy"))
          missingOption = "energy";
        if(!missingOption.empty()) {
          std::cerr << "Required option " << missingOption << " is missing." << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      } else {
        std::cout <<"Not performing a full run. This had better be a test..." << std::endl;
      }

      // -p/--projectile: projectile species
      projectileSpecies = ParticleSpecies(projectileString);
      if(projectileSpecies.theType == G4INCL::UnknownParticle && isFullRun) {
        std::cerr << "Error: unrecognized particle type " << projectileString << std::endl;
        std::cerr << suggestHelpMsg;
        std::exit(EXIT_FAILURE);
      }

      // -t/--target: target species
      if(variablesMap.count("target")) {
        targetSpecies = ParticleSpecies(targetString);
        if(targetSpecies.theType!=Composite) {
          std::cerr << "Unrecognized target. You specified: " << targetString << std::endl
            << "  The target nuclide must be specified in one of the following forms:" << std::endl
            << "    Fe56, 56Fe, Fe-56, 56-Fe, Fe_56, 56_Fe, Fe" << std::endl
            << "  You can also use IUPAC element names (such as Uuh)." << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
        if(targetSpecies.theA==0)
          naturalTarget = true;
      }

      // --pauli
      if(variablesMap.count("pauli")) {
        std::string pauliNorm = pauliString;
        std::transform(pauliNorm.begin(), pauliNorm.end(), pauliNorm.begin(), ::tolower);
        if(pauliNorm=="statistical")
          pauliType = StatisticalPauli;
        else if(pauliNorm=="strict")
          pauliType = StrictPauli;
        else if(pauliNorm=="strict-statistical")
          pauliType = StrictStatisticalPauli;
        else if(pauliNorm=="global")
          pauliType = GlobalPauli;
        else if(pauliNorm=="none")
          pauliType = NoPauli;
        else {
          std::cerr << "Unrecognized Pauli-blocking algorithm. Must be one of:" << std::endl
            << "  strict-statistical (default)" << std::endl
            << "  strict" << std::endl
            << "  statistical" << std::endl
            << "  global" << std::endl
            << "  none" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --coulomb
      if(variablesMap.count("coulomb")) {
        std::string coulombNorm = coulombString;
        std::transform(coulombNorm.begin(), coulombNorm.end(), coulombNorm.begin(), ::tolower);
        if(coulombNorm=="non-relativistic")
          coulombType = NonRelativisticCoulomb;
        else if(coulombNorm=="none")
          coulombType = NoCoulomb;
        else {
          std::cerr << "Unrecognized Coulomb-distortion algorithm. Must be one of:" << std::endl
            << "  non-relativistic (default)" << std::endl
            << "  none" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --potential
      if(variablesMap.count("potential")) {
        std::string potentialNorm = potentialString;
        std::transform(potentialNorm.begin(), potentialNorm.end(), potentialNorm.begin(), ::tolower);
        if(potentialNorm=="isospin-energy-smooth") {
          potentialType = IsospinEnergySmoothPotential;
        } else if(potentialNorm=="isospin-energy") {
          potentialType = IsospinEnergyPotential;
        } else if(potentialNorm=="isospin")
          potentialType = IsospinPotential;
        else if(potentialNorm=="constant")
          potentialType = ConstantPotential;
        else {
          std::cerr << "Unrecognized potential type. Must be one of:" << std::endl
            << "  isospin-energy-smooth" << std::endl
            << "  isospin-energy (default)" << std::endl
            << "  isospin" << std::endl
            << "  constant" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --local-energy-BB
      if(variablesMap.count("local-energy-BB")) {
        std::string localEnergyBBNorm = localEnergyBBString;
        std::transform(localEnergyBBNorm.begin(), localEnergyBBNorm.end(), localEnergyBBNorm.begin(), ::tolower);
        if(localEnergyBBNorm=="always") {
          localEnergyBBType = AlwaysLocalEnergy;
        } else if(localEnergyBBNorm=="first-collision")
          localEnergyBBType = FirstCollisionLocalEnergy;
        else if(localEnergyBBNorm=="never")
          localEnergyBBType = NeverLocalEnergy;
        else {
          std::cerr << "Unrecognized local-energy-BB type. Must be one of:" << std::endl
            << "  always" << std::endl
            << "  first-collision (default)" << std::endl
            << "  never" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --local-energy-pi
      if(variablesMap.count("local-energy-pi")) {
        std::string localEnergyPiNorm = localEnergyPiString;
        std::transform(localEnergyPiNorm.begin(), localEnergyPiNorm.end(), localEnergyPiNorm.begin(), ::tolower);
        if(localEnergyPiNorm=="always") {
          localEnergyPiType = AlwaysLocalEnergy;
        } else if(localEnergyPiNorm=="first-collision")
          localEnergyPiType = FirstCollisionLocalEnergy;
        else if(localEnergyPiNorm=="never")
          localEnergyPiType = NeverLocalEnergy;
        else {
          std::cerr << "Unrecognized local-energy-pi type. Must be one of:" << std::endl
            << "  always" << std::endl
            << "  first-collision" << std::endl
            << "  never (default)" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // -d/--de-excitation
      if(variablesMap.count("de-excitation")) {
        std::string deExcitationNorm = deExcitationString;
        std::transform(deExcitationNorm.begin(),
            deExcitationNorm.end(),
            deExcitationNorm.begin(), ::tolower);
        if(deExcitationNorm=="none")
          deExcitationType = DeExcitationNone;
#ifdef INCL_DEEXCITATION_ABLAXX
        else if(deExcitationNorm=="ablav3p")
          deExcitationType = DeExcitationABLAv3p;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
        else if(deExcitationNorm=="abla07")
          deExcitationType = DeExcitationABLA07;
#endif
#ifdef INCL_DEEXCITATION_SMM
        else if(deExcitationNorm=="smm")
          deExcitationType = DeExcitationSMM;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
        else if(deExcitationNorm=="geminixx")
          deExcitationType = DeExcitationGEMINIXX;
#endif
        else {
          std::cerr << "Unrecognized de-excitation model. "
            << "Must be one of:" << std::endl
            << deExcitationModelList << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      } else {
        deExcitationType = DeExcitationNone;
      }

      // --cluster-algorithm
      if(variablesMap.count("cluster-algorithm")) {
        std::string clusterAlgorithmNorm = clusterAlgorithmString;
        std::transform(clusterAlgorithmNorm.begin(),
            clusterAlgorithmNorm.end(),
            clusterAlgorithmNorm.begin(), ::tolower);
        if(clusterAlgorithmNorm=="none")
          clusterAlgorithmType = NoClusterAlgorithm;
        else if(clusterAlgorithmNorm=="intercomparison")
          clusterAlgorithmType = IntercomparisonClusterAlgorithm;
        else {
          std::cerr << "Unrecognized cluster algorithm. "
            << "Must be one of:" << std::endl
            << "  intercomparison (default)" << std::endl
            << "  none" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      } else {
        clusterAlgorithmType = IntercomparisonClusterAlgorithm;
      }

      // --cluster-max-mass
      if(variablesMap.count("cluster-max-mass") && clusterMaxMass < 2 && clusterMaxMass > 12) {
        std::cerr << "Maximum cluster mass outside the allowed range. Must be between 2 and 12 (included)"
          << std::endl
          << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
      }

      // --separation-energies
      if(variablesMap.count("separation-energies")) {
        std::string separationEnergyNorm = separationEnergyString;
        std::transform(separationEnergyNorm.begin(),
            separationEnergyNorm.end(),
            separationEnergyNorm.begin(), ::tolower);
        if(separationEnergyNorm=="incl")
          separationEnergyType = INCLSeparationEnergy;
        else if(separationEnergyNorm=="real")
          separationEnergyType = RealSeparationEnergy;
        else if(separationEnergyNorm=="real-light")
          separationEnergyType = RealForLightSeparationEnergy;
        else {
          std::cerr << "Unrecognized separation-energies option. "
            << "Must be one of:" << std::endl
            << "  INCL (default)" << std::endl
            << "  real" << std::endl
            << "  real-light" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      } else {
        separationEnergyType = INCLSeparationEnergy;
      }

      // --fermi-momentum
      if(variablesMap.count("fermi-momentum")) {
        std::string fermiMomentumNorm = fermiMomentumString;
        std::transform(fermiMomentumNorm.begin(),
            fermiMomentumNorm.end(),
            fermiMomentumNorm.begin(), ::tolower);
        if(fermiMomentumNorm=="constant")
          fermiMomentumType = ConstantFermiMomentum;
        else if(fermiMomentumNorm=="constant-light")
          fermiMomentumType = ConstantLightFermiMomentum;
        else if(fermiMomentumNorm=="mass-dependent")
          fermiMomentumType = MassDependentFermiMomentum;
        else {
          std::cerr << "Unrecognized fermi-momentum option. "
            << "Must be one of:" << std::endl
            << "  constant (default)" << std::endl
            << "  constant-light" << std::endl
            << "  mass-dependent" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      } else {
        fermiMomentumType = ConstantFermiMomentum;
      }

      // --rp-correlation / --rp-correlation-p / --rp-correlation-n
      if(variablesMap.count("rp-correlation")) {
        if(!variablesMap.count("rp-correlation-p") || variablesMap.find("rp-correlation-p")->second.defaulted())
          rpCorrelationCoefficientProton = rpCorrelationCoefficient;
        if(!variablesMap.count("rp-correlation-n") || variablesMap.find("rp-correlation-n")->second.defaulted())
          rpCorrelationCoefficientNeutron = rpCorrelationCoefficient;
      }

      // -s/--suffix
      if(!variablesMap.count("suffix")) {
        // update the value in the variables_map
        variablesMap.insert(std::make_pair("suffix", po::variable_value(boost::any(fileSuffix), false)));
      }

      // --output: construct a reasonable output file root if not specified
      if(!variablesMap.count("output") && isFullRun) {
        std::stringstream outputFileRootStream;
        // If an input file was specified, use its name as the output file root
        if(variablesMap.count("input-file"))
          outputFileRootStream << inputFileName << fileSuffix;
        else {
          outputFileRootStream.precision(0);
          outputFileRootStream.setf(std::ios::fixed, std::ios::floatfield);
          outputFileRootStream <<
            ParticleTable::getShortName(projectileSpecies) << "_" <<
            ParticleTable::getShortName(targetSpecies) << "_" <<
            projectileKineticEnergy;
          outputFileRootStream.precision(2);

          // Append suffixes to the output file root for each explicitly specified CLI option
          typedef po::variables_map::const_iterator BPOVMIter;
          for(BPOVMIter i=variablesMap.begin(), e=variablesMap.end(); i!=e; ++i) {
            std::string const &name = i->first;
            // Only process CLI options
            if(name!="projectile"
               && name!="target"
               && name!="energy"
               && name!="number-shots"
               && name!="random-seed-1"
               && name!="random-seed-2"
               && name!="verbosity"
               && name!="verbose-event"
               && name!="suffix"
#ifdef INCL_ROOT_USE
               && name!="root-selection"
#endif
               && name!="inclxx-datafile-path"
#ifdef INCL_DEEXCITATION_ABLA07
               && name!="abla07-datafile-path"
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
               && name!="ablav3p-cxx-datafile-path"
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
               && name!="geminixx-datafile-path"
#endif
              ) {
              po::variable_value v = i->second;
              if(!v.defaulted()) {
                const std::type_info &type = v.value().type();
                if(type==typeid(std::string))
                  outputFileRootStream << "_" << name << "=" << v.as<std::string>();
                else if(type==typeid(G4float))
                  outputFileRootStream << "_" << name << "=" << v.as<G4float>();
                else if(type==typeid(G4double))
                  outputFileRootStream << "_" << name << "=" << v.as<G4double>();
                else if(type==typeid(G4int))
                  outputFileRootStream << "_" << name << "=" << v.as<G4int>();
                else if(type==typeid(G4bool))
                  outputFileRootStream << "_" << name << "=" << v.as<G4bool>();
              }
            }
          }

          outputFileRootStream << fileSuffix;
        }

        // update the variable
        outputFileRoot = outputFileRootStream.str();
        // update the value in the variables_map
        variablesMap.insert(std::make_pair("output", po::variable_value(boost::any(outputFileRoot), false)));
      }

      // -l/--logfile
      if(!variablesMap.count("logfile")) {
        // update the variable
        logFileName = outputFileRoot + ".log";
        // update the value in the variables_map
        variablesMap.insert(std::make_pair("logfile", po::variable_value(boost::any(logFileName), false)));
      }

      // --random-seed-1 and 2
      if(!variablesMap.count("random-seed-1")) {
        if(randomSeed1<randomSeedMin || randomSeed1>randomSeedMax) {
          std::cerr << "Invalid value for random-seed-1. "
            << "Allowed range: [" << randomSeedMin << ", " << randomSeedMax << "]." << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }
      if(!variablesMap.count("random-seed-2")) {
        if(randomSeed2<randomSeedMin || randomSeed2>randomSeedMax) {
          std::cerr << "Invalid value for random-seed-2. "
            << "Allowed range: [" << randomSeedMin << ", " << randomSeedMax << "]." << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

    }
    catch(std::exception& e)
    {
      std::cerr << e.what() << "\n";
      std::cerr << suggestHelpMsg;
      std::exit(EXIT_FAILURE);
    }

  }
#else
    Config::Config(G4int /*argc*/, char * /*argv*/ [], G4bool /*isFullRun*/)
    {
      init();
    }
#endif

  Config::~Config()
  {}

  void Config::init() {
      verbosity = 1;
      inputFileName = "";
      title = "INCL default run title";
      nShots = 1000;
      naturalTarget = false;
      projectileString = "proton";
      projectileSpecies = G4INCL::Proton;
      projectileKineticEnergy = 1000.0;
      verboseEvent = -1;
      randomSeed1 = 666;
      randomSeed2 = 777;
      pauliString = "strict-statistical";
      pauliType = StrictStatisticalPauli;
      CDPP = true;
      coulombString = "non-relativistic";
      coulombType = NonRelativisticCoulomb;
      potentialString = "isospin-energy";
      potentialType = IsospinEnergyPotential;
      pionPotential = true;
      localEnergyBBString = "first-collision";
      localEnergyBBType = FirstCollisionLocalEnergy;
      localEnergyPiString = "first-collision";
      localEnergyPiType = FirstCollisionLocalEnergy;
      deExcitationString = "none";
      deExcitationType = DeExcitationNone;
      clusterAlgorithmString = "intercomparison";
      clusterAlgorithmType = IntercomparisonClusterAlgorithm;
      clusterMaxMass = 8;
      backToSpectator = true;
      useRealMasses = true;
      impactParameter = -1.;
      separationEnergyString = "INCL";
      separationEnergyType = INCLSeparationEnergy;
      fermiMomentumString = "constant";
      fermiMomentumType = ConstantFermiMomentum;
      cutNN = 1910.;
#ifdef INCL_DEEXCITATION_FERMI_BREAKUP
      maxMassFermiBreakUp = 18;
#endif
      rpCorrelationCoefficient = 1.;
      rpCorrelationCoefficientProton = 1.;
      rpCorrelationCoefficientNeutron = 1.;
      neutronSkinThickness = 0.;
      neutronSkinAdditionalDiffuseness = 0.;
      refraction=false;
  }

  std::string Config::summary() {
    std::stringstream message;
    message << "INCL++ version " << getVersionString() << std::endl;
    if(projectileSpecies.theType != Composite)
      message << "Projectile: " << ParticleTable::getName(projectileSpecies) << std::endl;
    else
      message << "Projectile: composite, A=" << projectileSpecies.theA << ", Z=" << projectileSpecies.theZ << std::endl;
    message << "  energy = " << projectileKineticEnergy << std::endl;
    if(targetSpecies.theA>0)
      message << "Target: A = " << targetSpecies.theA << " Z = " << targetSpecies.theZ << std::endl;
    else
      message << "Target: natural isotopic composition, Z = " << targetSpecies.theZ << std::endl;
    message << "Number of requested shots = " << nShots << std::endl;
    return message.str();
  }

#if defined(HAS_BOOST_PROGRAM_OPTIONS) && !defined(INCLXX_IN_GEANT4_MODE)
  std::string const Config::echo() const {
    std::stringstream ss;
    ss << "###########################\n"
      << "### Start of input echo ###\n"
      << "###########################\n\n"
      << "# You may re-use this snippet of the log file as an input file!\n"
      << "# Options marked with a * are compulsory.\n"
      << "\n### Run options\n" << echoOptionsDescription(runOptDesc)
      << "\n### Physics options\n" << echoOptionsDescription(physicsOptDesc)
      << "\n# the projectile nuclide was parsed as Z=" << projectileSpecies.theZ
      << ", A=" << projectileSpecies.theA
      << "\n# the target nuclide was parsed as Z=" << targetSpecies.theZ;
    if(targetSpecies.theA>0)
      ss << ", A=" << targetSpecies.theA;
    else
      ss << ", natural target";
    ss << "\n\n#########################\n"
      << "### End of input echo ###\n"
      << "#########################" << std::endl;

    return ss.str();
  }

  std::string Config::echoOptionsDescription(const po::options_description &aDesc) const {
    typedef std::vector< boost::shared_ptr< po::option_description > > OptVector;
    typedef std::vector< boost::shared_ptr< po::option_description > >::const_iterator OptIter;

    std::stringstream ss;
    ss << std::boolalpha;
    OptVector const &anOptVect = aDesc.options();
    for(OptIter opt=anOptVect.begin(), e=anOptVect.end(); opt!=e; ++opt) {
      std::string description = (*opt)->description();
      String::wrap(description);
      String::replaceAll(description, "\n", "\n# ");
      ss << "\n# " << description << std::endl;
      const std::string &name = (*opt)->long_name();
      ss << name << " = ";
      po::variable_value const &value = variablesMap.find(name)->second;
      std::type_info const &type = value.value().type();
      if(type == typeid(std::string)) {
        const std::string svalue = value.as<std::string>();
        if(svalue.empty())
          ss << "\"\"";
        else
          ss << svalue;
      } else if(type == typeid(G4int))
        ss << value.as<G4int>();
      else if(type == typeid(G4float))
        ss << value.as<G4float>();
      else if(type == typeid(G4double))
        ss << value.as<G4double>();
      else if(type == typeid(G4bool))
        ss << value.as<G4bool>();
      ss << '\n';
    }
    return ss.str();
  }
#endif

}
