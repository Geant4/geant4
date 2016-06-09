# Geant4MacroUtilities - this module defines helper macros and functions
#
# GEANT4_ADD_FEATURE(NAME DESCRIPTION)
#   Use this macro to add a Geant4 specific feature NAME, assumed to be a 
#   boolean, to the list of enabled/disabled features, together with a short 
#   DESCRIPTION.
#
# GEANT4_PRINT_ENABLED_FEATURES()
#   Prints list of enabled Geant4 features and their description only. Just a
#   simplified version of that in FeatureSummary.



macro(GEANT4_ADD_FEATURE _var _description)
    if(${_var})
        set_property(GLOBAL APPEND PROPERTY GEANT4_ENABLED_FEATURES ${_var})
    else(${_var})
        set_property(GLOBAL APPEND PROPERTY GEANT4_DISABLED_FEATURES ${_var})
    endif(${_var})

    set_property(GLOBAL PROPERTY ${_var}_DESCRIPTION "${_description}")
endmacro(GEANT4_ADD_FEATURE)


macro(GEANT4_PRINT_ENABLED_FEATURES)
    set(_currentFeatureText "The following Geant4 features are enabled:")
    get_property(_enabledFeatures GLOBAL PROPERTY GEANT4_ENABLED_FEATURES)

    foreach(_feature ${_enabledFeatures})
        set(_currentFeatureText "${_currentFeatureText}\n${_feature}")

        get_property(_desc GLOBAL PROPERTY ${_feature}_DESCRIPTION)

        if(_desc)
            set(_currentFeatureText "${_currentFeatureText}: ${_desc}")
            set(_desc NOTFOUND)
        endif(_desc)
    endforeach(_feature)

    message(STATUS "${_currentFeatureText}\n")
endmacro(GEANT4_PRINT_ENABLED_FEATURES)
