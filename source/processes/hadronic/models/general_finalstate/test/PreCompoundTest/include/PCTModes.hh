#ifndef PCTModes_hh
#define PCTModes_hh 1

typedef enum evap_modes
  {
    no_mode = 0,
    standard_mode,
    GEM_mode
  } PCTEVAPModes;

typedef enum preeq_emission_modes
  {
    no_emission_mode = 0,
    standard_emission_mode,
    HETC_emission_mode
  } PCTPREEQEmissionModes;

typedef enum preeq_transition_modes
  {
    no_transition_mode = 0,
    standard_transition_mode,
    GNASH_transition_mode
  } PCTPREEQTransitionModes;

#endif // PCTModes_hh
