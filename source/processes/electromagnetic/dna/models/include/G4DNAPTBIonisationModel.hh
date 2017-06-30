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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)
//

#ifndef G4DNAPTBIONISATIONMODEL_h
#define G4DNAPTBIONISATIONMODEL_h 1

#include "G4VDNAModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4DNACrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4DNAGenericIonsManager.hh"

#include "G4LogLogInterpolation.hh"

#include "G4DNAPTBIonisationStructure.hh"
#include "G4DNAPTBAugerModel.hh"
#include "G4NistManager.hh"

/*!
 * \brief The G4DNAPTBIonisationModel class
 * Implements the PTB ionisation model.
 */
class G4DNAPTBIonisationModel : public G4VDNAModel
{

public:
    /*!
     * \brief G4DNAPTBIonisationModel
     * Constructor
     * \param applyToMaterial
     * \param p
     * \param nam
     * \param isAuger
     */
    G4DNAPTBIonisationModel(const G4String &applyToMaterial = "all",
                            const G4ParticleDefinition* p = 0,
                            const G4String &nam = "DNAPTBIonisationModel",
                            const G4bool isAuger=true);

    /*!
     * \brief ~G4DNAPTBIonisationModel
     * Destructor
     */
    virtual ~G4DNAPTBIonisationModel();

    /*!
     * \brief Initialise
     * Method called once at the beginning of the simulation. It is used to setup the list of the materials managed by the model
     * and the energy limits. All the materials are setup but only a part of them can be activated by the user through the constructor.
     */
    virtual void Initialise(const G4ParticleDefinition* particle, const G4DataVector& = *(new G4DataVector()), G4ParticleChangeForGamma* fpChangeForGamme=nullptr);

    /*!
     * \brief CrossSectionPerVolume
     * Mandatory for every model the CrossSectionPerVolume method is in charge of returning the
     * cross section value corresponding to the material, particle and energy current values.
     * \param material
     * \param materialName
     * \param p
     * \param ekin
     * \param emin
     * \param emax
     * \return the cross section value
     */
    virtual G4double CrossSectionPerVolume(const G4Material* material,
                                           const G4String& materialName,
                                           const G4ParticleDefinition* p,
                                           G4double ekin,
                                           G4double emin,
                                           G4double emax);

    /*!
     * \brief SampleSecondaries
     * If the model is selected for the ModelInterface then SampleSecondaries will be called.
     * The method sets the characteristics of the particles implied with the physical process after the ModelInterface (energy, momentum...).
     * This method is mandatory for every model.
     * \param materialName
     * \param particleChangeForGamma
     * \param tmin
     * \param tmax
     */
    virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                   const G4MaterialCutsCouple*,
                                   const G4String& materialName,
                                   const G4DynamicParticle*,
                                   G4ParticleChangeForGamma *particleChangeForGamma,
                                   G4double tmin,
                                   G4double tmax);

protected:

private:

    G4DNAPTBAugerModel*      fDNAPTBAugerModel; ///< PTB Auger model instanciated in the constructor and deleted in the destructor of the class

    G4int verboseLevel; ///< verbose level

    G4DNAPTBIonisationStructure ptbStructure; /*!< ptbStructure class which contains the shell binding energies */

    typedef std::map<G4String, std::map<G4String, std::map<double, std::map<double, std::map<double, double> > > > > TriDimensionMap;
    TriDimensionMap diffCrossSectionData;
    TriDimensionMap fEnergySecondaryData;
    std::map<G4String, std::map<G4String, std::vector<double> > > fTMapWithVec;
    typedef std::map<G4String, std::map<G4String, std::map<double, std::vector<double> > > > VecMap;
    VecMap fEMapWithVector;
    typedef std::map<G4String, std::map<G4String, std::map<double, std::map<double, std::vector<double> > > > > VecMapWithShell;
    VecMapWithShell fProbaShellMap;

    G4double RandomizeEjectedElectronEnergy(G4ParticleDefinition * aParticleDefinition, G4double incomingParticleEnergy, G4int shell, const G4String& materialName);
    double DifferentialCrossSection(G4ParticleDefinition * aParticleDefinition, G4double k, G4double energyTransfer, G4int shell, const G4String &materialName);

    /*!
     * \brief RandomizeEjectedElectronEnergyFromCumulated
     * Uses the cumulated tables to find the energy of the ejected particle (electron)
     * \param particleDefinition
     * \param k
     * \param shell
     * \param materialName
     * \return the ejected electron energy
     */
    G4double RandomizeEjectedElectronEnergyFromCumulated(G4ParticleDefinition *particleDefinition, G4double k, G4int shell, const G4String& materialName);

    /*!
     * \brief RandomizeEjectedElectronDirection
     * Method to calculate the ejected electron direction
     * \param aParticleDefinition
     * \param incomingParticleEnergy
     * \param outgoingParticleEnergy
     * \param cosTheta
     * \param phi
     */
    void RandomizeEjectedElectronDirection(G4ParticleDefinition * aParticleDefinition, G4double incomingParticleEnergy, G4double
                                           outgoingParticleEnergy, G4double & cosTheta, G4double & phi );
    /*!
     * \brief ReadDiffCSFile
     * Method to read the differential cross section files.
     * \param materialName
     * \param particleName
     * \param file
     * \param scaleFactor
     */
    void ReadDiffCSFile(const G4String &materialName, const G4String &particleName, const G4String &file, const G4double scaleFactor);

    /*!
     * \brief QuadInterpolator
     * \param e11
     * \param e12
     * \param e21
     * \param e22
     * \param xs11
     * \param xs12
     * \param xs21
     * \param xs22
     * \param t1
     * \param t2
     * \param t
     * \param e
     * \return the interpolated value
     */
    G4double QuadInterpolator(G4double e11, G4double e12, G4double e21, G4double e22, G4double xs11, G4double xs12, G4double xs21, G4double xs22, G4double t1, G4double t2, G4double t, G4double e);
    /*!
     * \brief LogLogInterpolate
     * \param e1
     * \param e2
     * \param e
     * \param xs1
     * \param xs2
     * \return the interpolate value
     */
    G4double LogLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);

    // copy constructor and hide assignment operator
    G4DNAPTBIonisationModel(const  G4DNAPTBIonisationModel&); // prevent copy-construction
    G4DNAPTBIonisationModel & operator=(const  G4DNAPTBIonisationModel &right); // prevent assignement
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
