#ifndef G24ActionInitialization_HPP
#define G24ActionInitialization_HPP
#include"G4VUserActionInitialization.hh"
#include"G24PrimaryGeneratorAction.hpp"

class G24ActionInitialization: public G4VUserActionInitialization
{
    public:
    G24ActionInitialization();
    virtual ~G24ActionInitialization();
    virtual void Build() const override;
    G4double fbeta;
    G4double mass_c2;

};

#endif