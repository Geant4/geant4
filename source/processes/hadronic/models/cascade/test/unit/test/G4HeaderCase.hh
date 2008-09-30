#ifndef G4HeaderCase_HH
#define G4HeaderCase_HH

#include "../G4Unit.hh"

testSuite(HeaderSuite)

class CountLog: public BERT::Runner {
    unsigned m_count;
public:
    CountLog(): m_count(0) {}
    unsigned count() { return m_count; }
    virtual void enter(BERT::SuiteBase&) {}
    virtual void exit(BERT::SuiteBase&) {}
    virtual void check(BERT::CaseBase&) { ++m_count; }
};

// this case is implemented in a header file, which would be included by 2 source files
// the case should be counted only once
testCase(OnceCase, HeaderSuite)
{
    using namespace BERT;
    CountLog cl;
    suiteInstance<HeaderSuite>().run(cl);
    assertTrue(equalValueInfo(1u, cl.count()));
}

#endif // G4HeaderCase_HH
