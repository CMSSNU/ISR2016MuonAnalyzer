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

  // for SingleMuon trigger (not default)
  TString single_trig1 = "HLT_IsoMu24_v";
  TString single_trig2 = "HLT_IsoTkMu24_v";  
  vector<TString> trignames_single;
  trignames_single.push_back(single_trig1);
  trignames_single.push_back(single_trig2);

  // for DoubleMuon trigger (default)
  vector<TString> trignames_double;
  trignames_double.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  trignames_double.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  trignames_double.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  trignames_double.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  double lumiweight_single=1,lumiweight_double=1;
  if(!isData){
    lumiweight_single*=WeightByTrigger(trignames_single,TargetLumi); // change MC luminocity from 1pb^-1(default) to data luminocity
    lumiweight_double*=WeightByTrigger(trignames_double,TargetLumi);
  }
  weightGen=weight*lumiweight_double;

  FillCutFlow("NoCut", weightGen);
  // Get generator information for gen/reco ratio correction (only for Drell-Yan MC).
  // For this, we needs full-phase space generator level information.
  // So, You should run this with FLATCAT.
  bool ishardFSR=false;
  if(k_sample_name.Contains("DY")){  //only for Drell-Yan MC
    std::vector<snu::KTruth> truthcol =  eventbase->GetTruth();   //get truth particles
    TLorentzVector gendy,genl1,genl2,genl1pre,genl2pre;   //gendy: Z/gamma* fourvector, genl1: muon fourvector, genl2: anti-muon fourvector, *pre: before FSR
    int truthsize=truthcol.size();

    //loop for collect dimuon Drell-Yan product
    for(int i=0;i<truthsize;i++){
      snu::KTruth truth=truthcol.at(i);
      if(truth.GenStatus()!=1) continue;  //stable-particle-requirement
      if(truth.PdgId()==13&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl1+=truth;
      }else if(truth.PdgId()==-13&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl2+=truth;
      }else if(truth.PdgId()==22){   //collect photons from DY muons by FSR
	int imother=truth.IndexMother();
	if(imother>=truthsize) continue;
	snu::KTruth mother=truthcol.at(truth.IndexMother());
	if(mother.PdgId()==13&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
	  genl1pre+=truth;
	}else if(mother.PdgId()==-13&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
	  genl2pre+=truth;
	}
      }
    }
    genl1pre+=genl1;
    genl2pre+=genl2;
    gendy=genl1pre+genl2pre;

    //Fill generator level hists, if we find both of muon and anti-muon
    if(genl1!=TLorentzVector(0,0,0,0)&&genl2!=TLorentzVector(0,0,0,0)){  
      double gendimass = gendy.M();  
      double gendipt = gendy.Pt();
      double genl1pt = genl1.Pt();
      double genl2pt = genl2.Pt();
      double genl1eta = genl1.Eta();
      double genl2eta = genl2.Eta();
      double genl1mass = genl1.M();
      double genl2mass = genl2.M();
      if(genl1pt<genl2pt){
	double temppt=genl2pt; 
	double tempeta=genl2eta;
	double tempmass=genl2mass;
	genl2pt=genl1pt;
	genl2eta=genl1eta;
	genl2mass=genl1mass;
	genl1pt=temppt;
	genl1eta=tempeta;
	genl1mass=tempmass;
      }
      FillHists("gen","",gendimass,gendipt,-99999,-99999,genl1pt,genl2pt,genl1eta,genl2eta,weightGen);
      if(gendy.M()-(genl1+genl2).M()>5){
	ishardFSR=true;
	FillHists("gen","_hardFSR",gendimass,gendipt,-99999,-99999,genl1pt,genl2pt,genl1eta,genl2eta,weightGen);
      }
      ///////for unfolding///////
      ptPreFSR.push_back(genl1pre.Pt());
      ptPreFSR.push_back(genl2pre.Pt());
      ptPreFSR.push_back(gendipt);
      mPreFSR.push_back(genl1pre.M());
      mPreFSR.push_back(genl2pre.M());
      mPreFSR.push_back(gendimass);
      ptGen.push_back(genl1pt);
      ptGen.push_back(genl2pt);
      ptGen.push_back((genl1+genl2).Pt());
      mGen.push_back(0.000511);
      mGen.push_back(0.000511);
      mGen.push_back((genl1+genl2).M());
      if((genl1pre.Pt()>20||genl2pre.Pt()>20)&&(genl1pre.Pt()>10&&genl2pre.Pt()>10)&&fabs(genl1pre.Eta())<2.4&&fabs(genl2pre.Eta())<2.4){
	isfiducialPreFSR = 1;
      }
      if(genl1pt>20&&genl2pt>10&&fabs(genl1eta)<2.4&&fabs(genl2eta)<2.4){
	FillHists("gen","_fiducial",gendimass,gendipt,-99999,-99999,genl1pt,genl2pt,genl1eta,genl2eta,weight*lumiweight_double);
	isfiducialGen = 1;
      }
    }
  }  

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  FillCutFlow("BasicEventCut", weightGen);    

  std::vector<snu::KMuon> muons =GetMuons("MUON_POG_TIGHT",true); //get muon collection
  std::vector<snu::KElectron> electrons =  GetElectrons(false,false, "ELECTRON_NOCUT"); //get electron collection

  //select dimuon events
  if(muons.size()!=2) return;
  //if(electrons.size()>0) return;
  FillCutFlow("dimuonCut",weightGen);

  //check whether this event pass SingleMuon trigger
  bool passtrigger_single=true;
  if(!PassTriggerOR(trignames_single)) passtrigger_single=false;
  //check whether one of dimuon fire trigger
  if(!(muons[0].TriggerMatched(single_trig1)||muons[1].TriggerMatched(single_trig1)||muons[0].TriggerMatched(single_trig2)||muons[1].TriggerMatched(single_trig2))) passtrigger_single=false;
  
  //check whether this event pass DoubleMuon trigger
  bool passtrigger_double=PassTriggerOR(trignames_double);
  //check whether dimuon fire trigger
  bool matchtrigger_double=false;
  for(int i=0;i<4;i++){
    if(muons[0].TriggerMatched(trignames_double[i])&&muons[1].TriggerMatched(trignames_double[i])) matchtrigger_double=true;
  }
  passtrigger_double&=matchtrigger_double;

  //scale factors
  double PUreweight=1,PUreweight_up=1,PUreweight_down=1;
  double IDSF=1,IDSF_up=1,IDSF_down=1;
  double ISOSF=1,ISOSF_up=1,ISOSF_down=1;
  double triggerSF_single=1;
  if(!isData){
    PUreweight=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);  //pileup reweight
    PUreweight_up=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),1);  //pileup reweight sys
    PUreweight_down=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),-1);  //pileup reweight sys

    IDSF=mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muons, 0);  //muon ID scale factor
    IDSF_up=mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muons, 1);  //muon ID scale factor sys
    IDSF_down=mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muons, -1);  //muon ID scale factor sys

    ISOSF=mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muons, 0);  //muon isolation scale factor
    ISOSF_up=mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muons, 1);  //muon isolation scale factor sys
    ISOSF_down=mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muons, -1);  //muon isolation scale factor sys

    triggerSF_single=mcdata_correction->TriggerScaleFactor(electrons, muons, trignames_single[0],0);  //trigger scale factor for SingleMuon trigger
    //triggerSF_double=mcdata_correction->TriggerScaleFactor(electrons, muons, trignames_double[0],0);  //no scale factor for double muon trigger yet
  }
  weightTotal=weight*lumiweight_double*PUreweight*IDSF*ISOSF;
  FillCutFlow("AfterScaleFactor",weightTotal);

  //select opposite sign events
  if(muons.at(0).Charge()==muons.at(1).Charge()) return;
  FillCutFlow("osCut",weightTotal);
  
  //relative isolation cut
  if(muons.at(0).RelIso04() > 0.15) return;
  if(muons.at(1).RelIso04() > 0.15) return;
  FillCutFlow("RisoCut",weightTotal);

  
  CorrectMuonMomentum(muons);   //apply rochester correction
  CorrectedMETRochester(muons);  //update MET after rochester correction
  
  double dimass = (muons[0]+muons[1]).M();  
  double dipt = (muons[0]+muons[1]).Pt();
  double dieta = (muons[0]+muons[1]).Eta();
  double met=eventbase->GetEvent().MET();
  int nvtx=eventbase->GetEvent().nVertices();
  double l1pt = muons[0].Pt();
  double l2pt = muons[1].Pt();
  double l1eta = muons[0].Eta();
  double l2eta = muons[1].Eta();
  double l1mass = muons[0].M();
  double l2mass = muons[1].M();

  //dilepton mass cut
  //if(dimass<40||dimass>350) return;
  if(dimass<15) return;
  FillCutFlow("MassCut",weightTotal);

  //missing E_T cut
  //if(met>35) return;
  //FillCutFlow("METCut",weightTotal);

  //mark 'DY -> tau tau' events
  bool mcfromtau = (muons[0].MCFromTau()||muons[1].MCFromTau());
  DYtautau=mcfromtau;
  TString prefix="";
  if(mcfromtau&&k_sample_name.Contains("DY")) prefix="tau_";

  if(l1pt<l2pt){
    double temppt=l2pt; 
    double tempeta=l2eta;
    double tempmass=l2mass;
    l2pt=l1pt;
    l2eta=l1eta;
    l2mass=l1mass;
    l1pt=temppt;
    l1eta=tempeta;
    l1mass=tempmass;
  }

  //fill hists
  if(passtrigger_single){
    //lepton pT&eta cut
    if(l1pt>27&&l2pt>15&&(fabs(l1eta)<2.4)&&(fabs(l2eta)<2.4)){
      FillHists(prefix,"_SingleMuon",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_single*PUreweight*IDSF*ISOSF*triggerSF_single);
      if(ishardFSR) FillHists(prefix,"_SingleMuon_hardFSR",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_single*PUreweight*IDSF*ISOSF*triggerSF_single);
    }
  }
  if(passtrigger_double){
    //////////////Default/////////////////////
    if(l1pt>20&&l2pt>10&&(fabs(l1eta)<2.4)&&(fabs(l2eta)<2.4)){
      FillCutFlow("PtEtaCut",weightTotal);
      FillHists(prefix,"",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal); 
      if(ishardFSR) FillHists(prefix,"_hardFSR",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal);
     
      ///////////for unfolding tree////////////// 
      ispassRec=1;
      ptRec.push_back(l1pt);
      ptRec.push_back(l2pt);
      ptRec.push_back(dipt);
      etaRec.push_back(l1eta);
      etaRec.push_back(l2eta);
      etaRec.push_back(dieta);
      mRec.push_back(l1mass);
      mRec.push_back(l2mass);
      mRec.push_back(dimass);

      //////////////systematics/////////////////
      if(!isData){
	FillHists(prefix,"_sys_PUreweight_up",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_double*PUreweight_up*IDSF*ISOSF);
	FillHists(prefix,"_sys_PUreweight_down",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_double*PUreweight_down*IDSF*ISOSF);
	FillHists(prefix,"_sys_IDSF_up",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_double*PUreweight*IDSF_up*ISOSF);
	FillHists(prefix,"_sys_IDSF_down",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_double*PUreweight*IDSF_down*ISOSF);
	FillHists(prefix,"_sys_ISOSF_up",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_double*PUreweight*IDSF*ISOSF_up);
	FillHists(prefix,"_sys_ISOSF_down",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_double*PUreweight*IDSF*ISOSF_down);
	//scale variation
	vector<Float_t> scaleweight=eventbase->GetEvent().ScaleWeights();
	int imax=scaleweight.size();
	for(int i=0;i<imax;i++) FillHists(prefix,Form("_sys_scaleweight%d",i),dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_double*PUreweight*IDSF*ISOSF*scaleweight.at(i));
	//pdf sys
	//vector<Float_t> pdfweight=eventbase->GetEvent().PdfWeights();
	//int imax=pdfweight.size();
	//for(int i=0;i<imax;i++) FillHists(prefix,Form("_pdfweight%d",i),dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weight*lumiweight_double*PUreweight*IDSF*ISOSF*pdfweight.at(i));
      }
    }
    /////////For cut optimzation
    if(l1pt>20&&l2pt>15&&(fabs(l1eta)<2.4)&&(fabs(l2eta)<2.4)){
      FillHists(prefix,"_pt2015",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal);
      if(ishardFSR) FillHists(prefix,"_pt2015_hardFSR",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal);
    }
    if(l1pt>20&&l2pt>20&&(fabs(l1eta)<2.4)&&(fabs(l2eta)<2.4)){
      FillHists(prefix,"_pt2020",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal);
      if(ishardFSR) FillHists(prefix,"_pt2020_hardFSR",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal);
    }
    if(l1pt>20&&l2pt>10&&(fabs(l1eta)<1)&&(fabs(l2eta)<1)){
      FillHists(prefix,"_eta1",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal);
      if(ishardFSR) FillHists(prefix,"_eta1_hardFSR",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal);
    }
    if(l1pt>20&&l2pt>20&&(fabs(l1eta)<1)&&(fabs(l2eta)<1)){
      FillHists(prefix,"_pt2020_eta1",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal);
      if(ishardFSR) FillHists(prefix,"_pt2020_eta1_hardFSR",dimass,dipt,met,nvtx,l1pt,l2pt,l1eta,l2eta,weightTotal);
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

  DeclareVariable(ptPreFSR,"ptPreFSR","tree"); 
  DeclareVariable(mPreFSR,"mPreFSR","tree");
  DeclareVariable(etaRec,"etaRec","tree"); 
  DeclareVariable(ptRec,"ptRec","tree"); 
  DeclareVariable(mRec,"mRec","tree");
  DeclareVariable(ptGen,"ptGen","tree"); 
  DeclareVariable(etaGen,"etaGen","tree"); 
  DeclareVariable(mGen,"mGen","tree"); 
 
  DeclareVariable(ispassRec,"ispassRec","tree"); 
  DeclareVariable(isfiducialGen,"isfiducialGen","tree"); 
  DeclareVariable(isfiducialPreFSR,"isfiducialPreFSR","tree"); 

  DeclareVariable(weightGen,"weightGen","tree"); 
  DeclareVariable(weightTotal,"weightTotal","tree"); 
  DeclareVariable(DYtautau,"DYtautau","tree"); 

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

  isfiducialGen = 0;
  isfiducialPreFSR = 0;
  ispassRec = 0;
  DYtautau = 0;

  weightGen = 1.;
  weightTotal = 1.;

  etaRec.clear();
  ptRec.clear();
  mRec.clear();

  etaGen.clear();
  ptGen.clear();
  mGen.clear();

  ptPreFSR.clear();
  mPreFSR.clear();

}

void ISR2016MuonAnalyzer::FillProfile2D(TString histname, double x, double y, double z, double w, double xmin, double xmax, int nbinsx, double ymin, double ymax, int nbinsy , TString label, TString labely){
  m_logger << DEBUG << "FillProfile2D : " << histname << LQLogger::endmsg;
  if(GetHist2D(histname)) ((TProfile2D*)GetHist2D(histname))->Fill(x,y,z,w);
  else{
    if (nbinsx < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the profile" << LQLogger::endmsg;
    maphist2D[histname]=new TProfile2D(histname,histname,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
    ((TProfile2D*)maphist2D[histname])->Sumw2();
    if(GetHist2D(histname)) GetHist2D(histname)->GetXaxis()->SetTitle(label);
    if(GetHist2D(histname)) GetHist2D(histname)->GetYaxis()->SetTitle(labely);
    if(GetHist2D(histname)) ((TProfile2D*)GetHist2D(histname))->Fill(x,y,z,w);
  }
}

void ISR2016MuonAnalyzer::FillHists(TString prefix,TString suffix,double dimass,double dipt,double met,int nvtx,double l1pt,double l2pt,double l1eta,double l2eta,double ww){
  FillProfile2D(prefix+"profileM"+suffix,dimass,dipt,dimass,ww,0.,400.,400,0.,400.,400);
  FillProfile2D(prefix+"profilePt"+suffix,dimass,dipt,dipt,ww,0.,400.,400,0.,400.,400);
  if(l1pt<l2pt){
    double temppt=l2pt;
    l2pt=l1pt;
    l1pt=temppt;
    double tempeta=l2eta;
    l2eta=l1eta;
    l1eta=tempeta;
  }
  if(dipt<100){;
    if(dipt!=-99999) FillHist(prefix+"dipt"+suffix,dipt,ww,0.,100.,100);
    if(dimass!=-99999) FillHist(prefix+"dimass"+suffix,dimass,ww,0.,500.,500);
    if(met!=-99999) FillHist(prefix+"met"+suffix,met,ww,0.,100.,50);
    if(nvtx!=-99999) FillHist(prefix+"nvtx"+suffix,nvtx,ww,0.,50.,50);
    if(l1pt!=-99999) FillHist(prefix+"l1pt"+suffix,l1pt,ww,0.,100.,50);
    if(l2pt!=-99999) FillHist(prefix+"l2pt"+suffix,l2pt,ww,0.,100.,50);
    if(l1eta!=-99999) FillHist(prefix+"l1eta"+suffix,l1eta,ww,-4.,4.,80);
    if(l2eta!=-99999) FillHist(prefix+"l2eta"+suffix,l2eta,ww,-4.,4.,80);
    const float massrange[]={18,28,28,38,38,50,50,60,60,80,80,100,100,200,200,350};
    const int massbinnum=sizeof(massrange)/sizeof(float);
    for(int im=0;im<massbinnum;im++){
      if(dimass>=massrange[2*im]&&dimass<massrange[2*im+1]){
	if(dipt!=-99999) FillHist(Form("%sdipt_m%d%s",prefix.Data(),im,suffix.Data()),dipt,ww,0.,100.,100);
	if(dimass!=-99999) FillHist(Form("%sdimass_m%d%s",prefix.Data(),im,suffix.Data()),dimass,ww,0.,500.,500);
	if(met!=-99999) FillHist(Form("%smet_m%d%s",prefix.Data(),im,suffix.Data()),met,ww,0.,100.,50);
	if(nvtx!=-99999) FillHist(Form("%snvtx_m%d%s",prefix.Data(),im,suffix.Data()),nvtx,ww,0.,50.,50);
	if(l1pt!=-99999) FillHist(Form("%sl1pt_m%d%s",prefix.Data(),im,suffix.Data()),l1pt,ww,0.,100.,50);
	if(l2pt!=-99999) FillHist(Form("%sl2pt_m%d%s",prefix.Data(),im,suffix.Data()),l2pt,ww,0.,100.,50);
	if(l1eta!=-99999) FillHist(Form("%sl1eta_m%d%s",prefix.Data(),im,suffix.Data()),l1eta,ww,-4.,4.,80);
	if(l2eta!=-99999) FillHist(Form("%sl2eta_m%d%s",prefix.Data(),im,suffix.Data()),l2eta,ww,-4.,4.,80);
	break;
      }
    }
  }
}
