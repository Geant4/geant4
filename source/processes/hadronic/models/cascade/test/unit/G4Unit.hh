// A simple unit test aid for C++
// Based on unit-- by Tsong Chong

#ifndef G4UNIT_HH
#define G4UNIT_HH

#include <iostream>
#include <string>
#include <sstream>

namespace BERT {

class SuiteBase;
class CaseBase;
// log entire test process, providing simple control
class Runner {
public:
    virtual ~Runner() {}
    // enter() & exit() are called for every test suite
    virtual void enter(SuiteBase&) =0;
    virtual void exit(SuiteBase&) =0;
    // check() is called for every test case
    virtual void check(CaseBase&) =0;
};


// base class of tests
class TestBase {
    // tests will chain up as a linked list to avoid dependence on
    // dynamic memory allocation init
    virtual TestBase*& nextPointer() = 0;
    friend class SuiteBase;
public:
    virtual ~TestBase() {}
    virtual void run(Runner&) =0;
    virtual const char* caption() const =0;
};


struct SuiteState {
    bool added;
    TestBase* tests;
    TestBase* nextPointer;
};
// abstract class
// test suite defination
class SuiteBase: public TestBase {
protected:
    SuiteState* m_state;
private:
    // header of linked list
    TestBase*& tests()
    {
        return m_state->tests;
    }
    virtual TestBase*& nextPointer()
    {
        return m_state->nextPointer;
    }
public:
    void setState(SuiteState* statep)
    {
        m_state = statep;
    }
    virtual void run(Runner& log)
    {
        log.enter(*this);
        for (TestBase* p = tests(); p != 0; p = p->nextPointer()) {
            p->run(log);
        }
        log.exit(*this);
    }
    bool add(TestBase& test)
    {
        test.nextPointer() = tests();
        tests() = &test;
        return true;
    }
};

// implement singleton work
template <typename ThisSuite>
extern inline ThisSuite& suiteInstance()
{
    static SuiteState state = { false, 0, 0 };
    static ThisSuite s_inst;
    s_inst.setState(&state);
    return s_inst;
}
// abstract class
template <class ThisSuite, class UpperSuite>
class Suite: public SuiteBase {
    typedef Suite Self;
public:
    bool addToUpper()
    {
        // add will be performed only once
        return m_state->added = m_state->added || suiteInstance<UpperSuite>().add(*this);
    }
};

// the top most test suite
// every test case should be attached to Root directly or indirectly
// it's add to upper is not called, OR there will be a INFINITE LOOP
class Root: public Suite<Root, Root> {
public:
    virtual const char* caption() const
    {
        return "all tests";
    }
};


struct CaseState {
    bool added;
    TestBase* nextPointer;
};
// test work implementation
// abstract class
class CaseBase: public TestBase {
protected:
    CaseState* m_state;
private:
    virtual TestBase*& nextPointer()
    {
        return m_state->nextPointer;
    }
public:
    virtual void test() =0;
    virtual const char* file() const =0;
    virtual int line() const =0;
public:
    void setState(CaseState* statep)
    {
        m_state = statep;
    }
    virtual void run(Runner& log)
    {
        log.check(*this);
    }
};

// singleton work
template <typename ThisCase>
extern inline ThisCase& caseInstance()
{
    static CaseState state = { false, 0 };
    static ThisCase s_inst;
    s_inst.setState(&state);
    return s_inst;
}
// abstract class
// see Suite
template <class ThisCase, class UpperSuite>
class Case: public CaseBase {
private:
    typedef Case Self;
public:
    bool addToUpper()
    {
        return m_state->added = m_state->added || suiteInstance<UpperSuite>().add(*this);
    }
};


class AssertionInfo {
private:
    bool m_successful;
    bool m_hasMessage;
    std::string m_message;

    std::string m_rawMessage;
public:
    AssertionInfo(bool successp)
        : m_successful(successp), m_hasMessage(false) {}
    AssertionInfo(bool successp, const std::string& messagep)
        : m_successful(successp), m_hasMessage(true), m_message(messagep) {}
    bool successful() const { return m_successful; }
    operator bool () const { return successful(); }
    bool operator ! () const { return !successful(); }

    bool hasMessage() const { return m_hasMessage; }
    const std::string& message() const { return m_message; }

    const std::string& rawMessage() const { return m_rawMessage; }
    void rawMessage(const std::string& rawMess) { m_rawMessage = rawMess; }
};

class AssertionChecker;
// singleton, 7.1.2 - 4
extern inline AssertionChecker*& currentChecker()
{
    static AssertionChecker* s_ptr;
    return s_ptr;
}
class AssertionChecker {
    typedef AssertionChecker Self;
public:
    static void reg(Self& checker) { currentChecker() = &checker; }
    static Self& get() { return *currentChecker(); }

    virtual void before(const char* file, int line) =0;
    virtual void after() =0;
    virtual void check(const AssertionInfo& info) =0;
};
// make sure before() and after() are called
struct CheckerGuard {
    CheckerGuard(const char* file, int line)
        { AssertionChecker::get().before(file, line); }
    ~CheckerGuard() { AssertionChecker::get().after(); }
};

template <typename T1, typename T2>
inline bool equalValue(const T1& expected, const T2& actual)
{
    return expected == actual;
}
template <typename T1, typename T2>
inline AssertionInfo equalValueInfo(const T1& expected, const T2& actual)
{
    if (equalValue(expected, actual)) return AssertionInfo(true);

    std::ostringstream oo;
    oo << "expected <" << expected << ">, but was <" << actual << ">";
    return AssertionInfo(false, oo.str());
}

template <typename T1, typename T2, typename T3>
inline bool closeValue(const T1& expected, const T2& actual, const T3& precision)
{
    return !(expected + precision < actual)
        && !(actual + precision < expected);
}
template <typename T1, typename T2, typename T3>
inline AssertionInfo closeValueInfo(const T1& expected, const T2& actual, const T3& precision)
{
    if (closeValue(expected, actual, precision)) return AssertionInfo(true);

    std::ostringstream oo;
    oo << "expected <" << expected << ">, but was <" << actual << ">, not close enough with precision <" << precision << ">";
    return AssertionInfo(false, oo.str());
}

} // namespace BERT

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

// wrapper code for Suite<>
// a test suite should be a Suite<>
#define testSuitePrefix(thisSuite, upperSuite, prefix)                        \
class thisSuite: public BERT::Suite<thisSuite, upperSuite> {            \
    public: virtual const char* caption() const { return #thisSuite; }        \
};                                                                            \
const bool prefix##thisSuite                                                  \
    = BERT::suiteInstance<thisSuite>().addToUpper(); // 3.6.2 foot note 31

#define subSuite(thisSuite, upperSuite)                                       \
testSuitePrefix(thisSuite, upperSuite, added_)

#define testSuite(thisSuite) subSuite(thisSuite, BERT::Root)

// wrapper code for Case<>
// a test case should be a Case<>
#define testCasePrefix(thisCase, upperSuite, prefix)                          \
class thisCase: public BERT::Case<thisCase, upperSuite> {               \
public:                                                                       \
    virtual const char* caption() const { return #thisCase; }                 \
private:                                                                      \
    virtual const char* file() const { return __FILE__; }                     \
    virtual int line() const { return __LINE__; }                             \
    virtual void test();                                                      \
};                                                                            \
const bool prefix##thisCase                                                   \
    = BERT::caseInstance<thisCase>().addToUpper();                      \
inline void thisCase::test()

#define testCase(thisCase, upperSuite)                                        \
testCasePrefix(thisCase, upperSuite, added_)


// test assertion
// param: expression to test
// true for pass, false for failure
#define assertTrue(opt_expression)                                            \
{                                                                             \
    BERT::CheckerGuard opt_tempGard(__FILE__, __LINE__);                \
    BERT::AssertionInfo opt_info(opt_expression);                       \
    opt_info.rawMessage("Assertion failed: <" #opt_expression ">");           \
    BERT::AssertionChecker::get().check(opt_info);                      \
}

// test assertion
// assume failure when executed
// param: error message
#define assertNoArrive(opt_message)                                           \
{                                                                             \
    BERT::CheckerGuard opt_tempGard(__FILE__, __LINE__);                \
    BERT::AssertionInfo opt_info(false, (opt_message));                 \
    opt_info.rawMessage(opt_message);                                         \
    BERT::AssertionChecker::get().check(opt_info);                      \
}

#endif // G4UNIT_HH
