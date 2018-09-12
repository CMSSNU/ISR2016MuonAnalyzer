#ifndef ISR2016MuonAnalyzer_h
#define ISR2016MuonAnalyzer_h

#include "AnalyzerCore.h"
class ISR2016MuonAnalyzer : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  ISR2016MuonAnalyzer();
  ~ISR2016MuonAnalyzer();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( ISR2016MuonAnalyzer, 1);
};
#endif
