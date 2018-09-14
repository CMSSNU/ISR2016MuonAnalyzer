// $Id: ISR2016MuonAnalyzer.cc 1 2018-8-6 21:58 hsseo $
/***************************************************************************
 * @Project: LQISR2016MuonAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author HyonSan Seo       <hyon.san.seo@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ISR2016MuonAnalyzer.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ISR2016MuonAnalyzer);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ISR2016MuonAnalyzer::ISR2016MuonAnalyzer() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ISR2016MuonAnalyzer");
  
  Message("In ISR2016MuonAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //MakeCleverHistograms(sighist_mm,"DiMuon");


}


void ISR2016MuonAnalyzer::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  //Message("Making clever hists for Z ->ll test code", INFO);
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:


  return;
}


void ISR2016MuonAnalyzer::ExecuteEvents()throw( LQError ){
  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;   

  // Get generator information for gen/reco ratio correction (only for Drell-Yan MC).
  // For this, we needs full-phase space generator level information.
  // So, You should run this with FLATCAT.
  if(k_sample_name.Contains("DY")){  //only for Drell-Yan MC
    std::vector<snu::KTruth> truthcol =  eventbase->GetTruth();   //get truth particles
    TLorentzVector gendy,genl1,genl2;   //gendy: Z/gamma* fourvector, genl1: muon fourvector, genl2: anti-muon fourvector
    int truthsize=truthcol.size();

    //loop for collect dimuon Drell-Yan product
    for(int i=0;i<truthsize;i++){
      snu::KTruth truth=truthcol.at(i);
      if(truth.GenStatus()!=1) continue;  //stable-particle-requirement
      if(truth.PdgId()==13&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl1+=truth;
	gendy+=truth;
      }else if(truth.PdgId()==-13&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl2+=truth;
	gendy+=truth;
      }else if(truth.PdgId()==22){   //collect photons from DY muons by FSR
	int imother=truth.IndexMother();
	if(imother>=truthsize) continue;
	snu::KTruth mother=truthcol.at(truth.IndexMother());
	if(abs(mother.PdgId())==13&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
	  gendy+=truth;
	}
      }
    }
    
    //Fill generator level hists, if we find both of muon and anti-muon
    if(genl1!=TLorentzVector(0,0,0,0)&&genl2!=TLorentzVector(0,0,0,0)){  
      double genptlep1 = genl1.Pt();
      double genptlep2 = genl2.Pt();
      double genetalep1 = genl1.Eta();
      double genetalep2 = genl2.Eta();
      double gendipt = gendy.Pt();
      double gendimass = gendy.M();  
      
      cout<<"before filling gen2D"<<endl;
      ((TProfile2D*)maphist2D["genprofileM"])->Fill(gendimass,gendipt,gendimass,weight);
      ((TProfile2D*)maphist2D["genprofilePt"])->Fill(gendimass,gendipt,gendipt,weight);
      FillHist("gendipt",gendipt,weight,0.,100.,100);
      FillHist("gendimass",gendimass,weight,0.,400.,200);
      if(genptlep1>genptlep2){
	FillHist("genl1pt",genptlep1,weight,0.,100.,50);
	FillHist("genl2pt",genptlep2,weight,0.,100.,50);
	FillHist("genl1eta",genetalep1,weight,-4.,4.,80);
	FillHist("genl2eta",genetalep2,weight,-4.,4.,80);
      }else{
	FillHist("genl1pt",genptlep2,weight,0.,100.,50);
	FillHist("genl2pt",genptlep1,weight,0.,100.,50);
	FillHist("genl1eta",genetalep2,weight,-4.,4.,80);
	FillHist("genl2eta",genetalep1,weight,-4.,4.,80);
      }    
      const float massrange[]={40,60,60,80,80,100,100,200,200,350};
      for(int im=0;im<5;im++){
	if(gendimass>=massrange[2*im]&&gendimass<massrange[2*im+1]){
	  FillHist(Form("gendipt_m%d",im),gendipt,weight,0.,100.,100);
	  FillHist(Form("gendimass_m%d",im),gendimass,weight,0.,400.,200);
	  break;
	}
      }
    }
  }
  
  // use singlemuon trigger
  TString single_trig1 = "HLT_IsoMu24_v";
  TString single_trig2 = "HLT_IsoTkMu24_v";  
  vector<TString> trignames;
  trignames.push_back(single_trig1);
  trignames.push_back(single_trig2);
  if(!isData) weight*=WeightByTrigger(trignames,TargetLumi); // change MC luminocity from 1pb^-1(default) to data luminocity

  FillCutFlow("NoCut", weight);

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  FillCutFlow("BasicEventCut", weight);    

  std::vector<snu::KMuon> muons =GetMuons("MUON_POG_TIGHT",true); //get muon collection
  std::vector<snu::KElectron> electrons =  GetElectrons(false,false, "ELECTRON_NOCUT"); //get electron collection

  //select dimuon events
  if(muons.size()!=2) return;
  FillCutFlow("dimuonCut",weight);

  //check whether this event pass single muon trigger
  if(!PassTriggerOR(trignames)) return;
  //check whether one of dimuon fire trigger
  if(!(muons[0].TriggerMatched(single_trig1)||muons[1].TriggerMatched(single_trig1)||muons[0].TriggerMatched(single_trig2)||muons[1].TriggerMatched(single_trig2))) return;
  FillCutFlow("triggerCut",weight);  
  
  //apply scale factors
  if(!isData){
    weight*=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);  //apply pileup reweight
    weight*=mcdata_correction->TriggerScaleFactor(electrons, muons, trignames[0],0);  //apply trigger scale factor
    weight*=mcdata_correction->MuonTrackingEffScaleFactor(muons); //apply muon tracking efficiency scale factor
    weight*=mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muons, 0);  //apply muon ID scale factor
    weight*=mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muons, 0);  //apply muon isolation scale factor
  }
  FillCutFlow("AfterScaleFactor",weight);

  //select opposite sign events
  if(muons.at(0).Charge()==muons.at(1).Charge()) return;
  FillCutFlow("osCut",weight);
  
  //relative isolation cut
  if(muons.at(0).RelIso04() > 0.15) return;
  if(muons.at(1).RelIso04() > 0.15) return;
  FillCutFlow("RisoCut",weight);

  
  CorrectMuonMomentum(muons);   //apply rochester correction
  CorrectedMETRochester(muons);  //update MET after rochester correction

  double met=eventbase->GetEvent().MET();
  double ptlep1 = muons[0].Pt();
  double ptlep2 = muons[1].Pt();
  double etalep1 = muons[0].Eta();
  double etalep2 = muons[1].Eta();
  double dipt = (muons[0]+muons[1]).Pt();
  double dimass = (muons[0]+muons[1]).M();  

  //dilepton mass cut
  if(dimass<40||dimass>350) return;
  FillCutFlow("MassCut",weight);

  //lepton pT&eta cut
  if(!(((ptlep1>27)||(ptlep2>27))&&(ptlep1>15)&&(ptlep2>15)&&(fabs(etalep1)<2.4)&&(fabs(etalep2)<2.4))) return;
  FillCutFlow("PtEtaCut",weight);

  //dilepton pT cut
  if(dipt>100) return;
  FillCutFlow("DileptonPtCut",weight);

  //missing E_T cut
  //if(met>35) return;
  //FillCutFlow("METCut",weight);

  //mark 'DY -> tau tau' events
  bool mcfromtau = (muons[0].MCFromTau()||muons[1].MCFromTau());
  TString prefix="";
  if(mcfromtau&&k_sample_name.Contains("DY")) prefix="tau_";
  
  //fill hists
  cout<<"before filling 2D"<<endl;
  ((TProfile2D*)maphist2D["profileM"])->Fill(dimass,dipt,dimass,weight);
  ((TProfile2D*)maphist2D["profilePt"])->Fill(dimass,dipt,dipt,weight);
  FillHist(prefix+"dipt",dipt,weight,0.,100.,100);
  FillHist(prefix+"dimass",dimass,weight,0.,400.,200);
  FillHist(prefix+"met",met,weight,0.,100.,50);
  FillHist(prefix+"nvtx",eventbase->GetEvent().nVertices(),weight,0.,50.,50);
  if(ptlep1>ptlep2){
    FillHist(prefix+"l1pt",ptlep1,weight,0.,100.,50);
    FillHist(prefix+"l2pt",ptlep2,weight,0.,100.,50);
    FillHist(prefix+"l1eta",etalep1,weight,-4.,4.,80);
    FillHist(prefix+"l2eta",etalep2,weight,-4.,4.,80);
  }else{
    FillHist(prefix+"l1pt",ptlep2,weight,0.,100.,50);
    FillHist(prefix+"l2pt",ptlep1,weight,0.,100.,50);
    FillHist(prefix+"l1eta",etalep2,weight,-4.,4.,80);
    FillHist(prefix+"l2eta",etalep1,weight,-4.,4.,80);
  }    
  const float massrange[]={40,60,60,80,80,100,100,200,200,350};
  for(int im=0;im<5;im++){
    if(dimass>=massrange[2*im]&&dimass<massrange[2*im+1]){
      FillHist(Form("%sdipt_m%d",prefix.Data(),im),dipt,weight,0.,100.,100);
      FillHist(Form("%sdimass_m%d",prefix.Data(),im),dimass,weight,0.,400.,200);
      FillHist(Form("%smet_m%d",prefix.Data(),im),met,weight,0.,100.,50);
      FillHist(Form("%snvtx_m%d",prefix.Data(),im),eventbase->GetEvent().nVertices(),weight,0.,50.,50);
      if(ptlep1>ptlep2){
	FillHist(Form("%sl1pt_m%d",prefix.Data(),im),ptlep1,weight,0.,100.,50);
	FillHist(Form("%sl2pt_m%d",prefix.Data(),im),ptlep2,weight,0.,100.,50);
	FillHist(Form("%sl1eta_m%d",prefix.Data(),im),etalep1,weight,-4.,4.,80);
	FillHist(Form("%sl2eta_m%d",prefix.Data(),im),etalep2,weight,-4.,4.,80);
      }else{
	FillHist(Form("%sl1pt_m%d",prefix.Data(),im),ptlep2,weight,0.,100.,50);
	FillHist(Form("%sl2pt_m%d",prefix.Data(),im),ptlep1,weight,0.,100.,50);
	FillHist(Form("%sl1eta_m%d",prefix.Data(),im),etalep2,weight,-4.,4.,80);
	FillHist(Form("%sl2eta_m%d",prefix.Data(),im),etalep1,weight,-4.,4.,80);
      }    
      break;
    }
  }

  return;
}// End of execute event loop
  


void ISR2016MuonAnalyzer::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ISR2016MuonAnalyzer::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

ISR2016MuonAnalyzer::~ISR2016MuonAnalyzer() {
  
  Message("In ISR2016MuonAnalyzer Destructor" , INFO);
  
}


void ISR2016MuonAnalyzer::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ISR2016MuonAnalyzer::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ISR2016MuonAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/
  if(k_sample_name.Contains("DY")){
    maphist2D["genprofileM"]=new TProfile2D("genprofileM","genprofileM",200,0,400,200,0,400);
    maphist2D["genprofilePt"]=new TProfile2D("genprofilePt","genprofilePt",200,0,400,200,0,400);
  }
  maphist2D["profileM"]=new TProfile2D("profileM","profileM",200,0,400,200,0,400);
  maphist2D["profilePt"]=new TProfile2D("profilePt","profilePt",200,0,400,200,0,400);
}


void ISR2016MuonAnalyzer::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



