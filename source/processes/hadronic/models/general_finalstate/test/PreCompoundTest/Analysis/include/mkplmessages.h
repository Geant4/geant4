#ifndef mkplmessages_h
#define mkplmessages_h

enum MKPL_MSG
  {
    ENERGY_SELECTOR_MSG        = 1001,
    ENERGY_SELECTOR_DELETED    = 1,
    ENERGY_SELECTOR_DE         = 5,
    ENERGY_SELECTOR_DA         = 6,
    ENERGY_SELECTOR_DD         = 7,
    ENERGY_SELECTOR_DDA        = 8,
    ANGLE_SELECTOR_MSG         = 2001,
    ANGLE_SELECTOR_DELETED     = 1,
    ANGLE_SELECTOR_ONE         = 2,
    ANGLE_SELECTOR_ALL         = 3,
    ANGLE_SELECTOR_CLOSE       = 4,
    RANGE_SELECTOR_MSG         = 3001,
    RANGE_SELECTOR_DELETED     = 1,
    RANGE_SELECTOR_ONE         = 2,
    RANGE_SELECTOR_CLOSE       = 4,
    COMPARISON_SELECTOR_MSG    = 4001,
    COMPARISON_SELECTOR_DELETED= 1,
    COMPARISON_SELECTOR_COMP   = 2,
    COMPARISON_SELECTOR_TEST   = 4
  };

inline Int_t mkpl_mkmsg(MKPL_MSG msg, MKPL_MSG submsg) 
{
  return Int_t(msg << 8) + submsg; 
}

#endif
