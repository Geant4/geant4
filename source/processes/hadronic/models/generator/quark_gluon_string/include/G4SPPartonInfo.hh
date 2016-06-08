#ifndef G4SPPartonInfo_h
#define G4SPPartonInfo_h

class G4SPPartonInfo
{
  public:
    G4SPPartonInfo(G4int diq, G4int q, G4double prob) 
    { diQuarkPDGCode = diq; quarkPDGCode = q; probability = prob; }
    G4int GetQuark() const {return quarkPDGCode;}
    G4int GetDiQuark() const {return diQuarkPDGCode;}
    G4double GetProbability() const {return probability;}      
    G4bool operator == (const G4SPPartonInfo & aInfo) const
    {return this == &aInfo;}
  private:      
    G4int quarkPDGCode;
    G4int diQuarkPDGCode;
    G4double probability;
};
    
#endif
