#ifndef STORAGEMANAGER_HH
#define STORAGEMANAGER_HH

#include <string>
#include <fstream>

namespace AIDA {
   class IAnalysisFactory;
   class ITree;
   class IHistogramFactory;
   class IHistogram1D;
   class IHistogram2D;
}


class StorageManager {

 public:
   static StorageManager* Instance(std::string fileNameBase);
   static void Destroy(); 

   AIDA::IHistogramFactory* GetIHistogramFactory();
   std::ostream& GetOutputFileStream();

 protected:
   virtual ~StorageManager();
   StorageManager(std::string fileNameBase);

   StorageManager(const StorageManager& only);
   const StorageManager& operator=(const StorageManager& only);

 private:
   static StorageManager* instance;

   AIDA::IHistogramFactory* histogramFactory;
   AIDA::IAnalysisFactory* analysisFactory;
   AIDA::ITree* tree;

   std::ofstream os;
};

#endif // STORAGEMANAGER_HH
