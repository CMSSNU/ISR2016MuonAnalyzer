#ifndef ISR2016MuonAnalyzer_h
#define ISR2016MuonAnalyzer_h

#include "AnalyzerCore.h"
#include "TProfile2D.h"
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
  void FillProfile2D(TString histname,double x,double y,double z,double w,double xmin,double xmax,int nbinsx,double ymin,double ymax,int nbinsy,TString lable="",TString labely="");
  void FillHists(TString prefix,TString suffix,double dimass,double dipt,double met,int nvtx,double l1pt,double l2pt,double l1eta,double l2eta,double ww);
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
