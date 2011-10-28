#ifndef G4EXCEPTIONORIGIN_HH
#define G4EXCEPTIONORIGIN_HH

#define __Exception_Origin__ \
    G4String exceptionOrigin(__FILE__); \
    exceptionOrigin += " ("; \
    std::ostringstream os; \
    os << __LINE__; \
    exceptionOrigin += os.str(); \
    exceptionOrigin += ")";

#endif // G4EXCEPTIONORIGIN_HH
