#ifndef ProtectedMemory_h
#define ProtectedMemory_h 1

#ifdef __cplusplus
extern "C" {
#endif

void BuildProtectedMemory(size_t totalHeap);
void setuphandlers();

void addThread();
void starttracer();
void finishtracer();
void StartDetection();
void FinishDetection();
void PhaseChange();

#ifdef __cplusplus
}
#endif

#endif
