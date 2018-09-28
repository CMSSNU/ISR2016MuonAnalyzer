#include<vector>
TString hsseopath="/data2/CAT_SKTreeOutput/JobOutPut/hsseo/LQanalyzer/data/output/CAT/ISR2016MuonAnalyzer/periodBtoH/";

////////////////////////////////////////////////
////////////////////TProfile2D functions////////
///////////////////////////////////////////////
bool TProfile2DAdd(TProfile2D* profile1,TProfile2D* profile2,double c=1.){
  // Performs the operation: profile1 += c*profile2                                                                                
  // delete buffer if it is there since it will become invalid                                                                     
  if (profile1->fBuffer) profile1->BufferEmpty(1);

  // Check profile compatibility                                                                                                   
  Int_t nx = profile1->GetNbinsX();
  Int_t ny = profile1->GetNbinsY();
  Int_t nz = profile1->GetNbinsZ();

  if ( nx != profile1->GetNbinsX() ||  nx != profile2->GetNbinsX() ||
       ny != profile1->GetNbinsY() ||  ny != profile2->GetNbinsY() ||
       nz != profile1->GetNbinsZ() ||  nz != profile2->GetNbinsZ() ) {
    cout<<"Error::TProfile2DAdd::Attempt to add profiles with different number of bins"<<endl;
    return kFALSE;
  }

  // Add statistics                                                                                                                 
  //Double_t ac = TMath::Abs(c);
  profile1->SetEntries(profile1->GetEntries() +c*profile2->GetEntries());
  Double_t s1[TH1::kNstat], s2[TH1::kNstat];
  Int_t i;
  for (i=0;i<TH1::kNstat;i++) {s1[i] = s2[i] = 0;}
  profile1->GetStats(s1);
  profile2->GetStats(s2);
  for (i=0;i<TH1::kNstat;i++) {
    if (i == 1) s1[i] += c*c*s2[i];
    else        s1[i] += c*s2[i];
  }
  profile1->PutStats(s1);
  
  // create sumw2 per bin if not set                                                                                                
  if (profile1->GetBinSumw2()->fN == 0) profile1->Sumw2();
  if (profile2->GetBinSumw2()->fN == 0) profile2->Sumw2();

  // Make the loop over the bins to calculate the Addition                                                                          
  Double_t *p1sumwz = profile1->GetArray();    Double_t *p2sumwz = profile2->GetArray();
  Double_t *p1sumwz2 = profile1->GetSumw2()->GetArray();   Double_t *p2sumwz2 = profile2->GetSumw2()->GetArray();
  //Double_t *p1sumw = profile1->GetBinArray();    Double_t *p2sumw = profile2->GetBinArray();
  Double_t *p1sumw2 = profile1->GetBinSumw2()->GetArray();   Double_t *p2sumw2 = profile2->GetBinSumw2()->GetArray();
  for (int bin =0;bin< profile1->fN;bin++) {
    p1sumwz[bin]+=c*p2sumwz[bin];
    p1sumwz2[bin]+=c*p2sumwz2[bin];
    double tempnewsumw2=p1sumw2[bin]+c*c*p2sumw2[bin];
    profile1->SetBinEntries(bin,profile1->GetBinEntries(bin)+c*profile2->GetBinEntries(bin));
    p1sumw2[bin]=tempnewsumw2;
  }
  return kTRUE;
}
void TProfile2DGetStats(TProfile2D* profile,double* stats,double xmin,double xmax,double ymin,double ymax){
  const int statnum=4;
  int xminbin=profile->GetXaxis()->FindBin(xmin);
  int xmaxbin=profile->GetXaxis()->FindBin(xmax);
  int yminbin=profile->GetYaxis()->FindBin(ymin);
  int ymaxbin=profile->GetYaxis()->FindBin(ymax);
  for(int i=0;i<statnum;i++) stats[i]=0;
  double* Asumw2=profile->GetBinSumw2()->GetArray();
  double* Asumwz=profile->GetArray();
  double* Asumwz2=profile->GetSumw2()->GetArray();
  for(int ixbin=xminbin;ixbin<xmaxbin;ixbin++){
    for(int iybin=yminbin;iybin<ymaxbin;iybin++){
      int ibin=profile->GetBin(ixbin,iybin);
      stats[0]+=profile->GetBinEntries(ibin);
      stats[1]+=Asumw2[ibin];
      stats[2]+=Asumwz[ibin];
      stats[3]+=Asumwz2[ibin];
    }
  }
}
double TProfile2DGetSumwz(TProfile2D* profile,double xmin,double xmax,double ymin,double ymax){
  double stats[4];
  TProfile2DGetStats(profile,stats,xmin,xmax,ymin,ymax);
  return stats[2];
}
double TProfile2DGetSumw(TProfile2D* profile,double xmin,double xmax,double ymin,double ymax){
  double stats[4];
  TProfile2DGetStats(profile,stats,xmin,xmax,ymin,ymax);
  return stats[0];
}
double TProfile2DGetMean(TProfile2D* profile,double xmin,double xmax,double ymin,double ymax){
  double stats[4];
  TProfile2DGetStats(profile,stats,xmin,xmax,ymin,ymax);
  return stats[2]/stats[0];
}
double TProfile2DGetMeanError(TProfile2D* profile,double xmin, double xmax,double ymin,double ymax){
  double stats[4];
  TProfile2DGetStats(profile,stats,xmin,xmax,ymin,ymax);
  return sqrt((stats[3]/stats[0]-pow(stats[2]/stats[0],2))/(pow(stats[0],2)/stats[1]));
}

///////////////////////////////////////////////////////////
////////////////////Histogram plotting functions////////
////////////////////////////////////////////////////////
TH1* get_hist_data_bgsub(TString dirname, TString histname,TString tdir="Hists/"){
  TSystemDirectory dir(dirname,dirname);
  vector<TH1*> hists;
  TH1* hdata;
  TString stream="DoubleMuon";
  if(histname.Contains("SingleMuon")) stream="SingleMuon";
  for(int i=0;i<dir.GetListOfFiles()->GetSize();i++){
    TSystemFile *file=(TSystemFile*)dir.GetListOfFiles()->At(i);
    if(file->IsDirectory()) continue;
    bool isdata=0;
    if(strstr(file->GetName(),"data")!=NULL) isdata=1;
    if(isdata&&strstr(file->GetName(),stream.Data())==NULL) continue;
    TFile tfile(TString(file->GetTitle())+TString(file->GetName()));
    TH1* hist=(TH1*)tfile.Get(tdir+histname);
    if(!hist) continue;
    hist->GetXaxis()->SetTitle(hist->GetName());
    TString title=file->GetName();
    TString subtitle=title(title.Index("_")+1,title.Index("cat")-title.Index("_")-2);
    hist->SetTitle(subtitle);
    hist->SetDirectory(0);
    if(isdata) hdata=hist;
    else if(strstr(file->GetName(),"DY")!=NULL){
      TH1* tauhist=(TH1*)tfile.Get(tdir+"tau_"+histname);
      if(!tauhist) continue;
      tauhist->GetXaxis()->SetTitle(hist->GetName());
      tauhist->SetTitle("tautau");
      tauhist->SetDirectory(0);
      hists.push_back(tauhist);
    }
    else hists.push_back(hist);
  }
  for(int i=0;i<hists.size();i++){
    TH1* hist=hists.at(i);
    if(strcmp(hdata->ClassName(),"TProfile2D")==0) TProfile2DAdd((TProfile2D*)hdata,(TProfile2D*)hist,-1.);
    else hdata->Add(hist,-1.);
  }
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(0.6);
  return hdata;
}
TH1* get_hist_mc_dy(TString dirname, TString histname,TString tdir="Hists/"){
  TSystemDirectory dir(dirname,dirname);
  TH1* hdy;
  for(int i=0;i<dir.GetListOfFiles()->GetSize();i++){
    TSystemFile *file=(TSystemFile*)dir.GetListOfFiles()->At(i);
    if(file->IsDirectory()) continue;
    if(strstr(file->GetName(),"DY")!=NULL){
      TFile tfile(TString(file->GetTitle())+TString(file->GetName()));
      TH1* hist=(TH1*)tfile.Get(tdir+histname);
      if(!hist) continue;
      hist->GetXaxis()->SetTitle(hist->GetName());
      TString title=file->GetName();
      TString subtitle=title(title.Index("_")+1,title.Index("cat")-title.Index("hsseo")-2);
      hist->SetTitle(subtitle);
      hist->SetDirectory(0);
      if(!hdy) hdy=hist;
      else hdy->Add(hist,1);
    }
  }
  hdy->SetTitle("DY");
  return hdy;
}

///////////////////////////////////////////////////////////
////////////////////graph functions///////////////////////
//////////////////////////////////////////////////////////
TGraphErrors* get_graph(TString dirname,void *histfunction,bool getgenlevel=0,TString suffix=""){
  //const double massrange[]={40,60,60,80,80,100,100,200,200,350};
  const double massrange[]={40,50,50,60,60,70,70,80,80,100,100,130,130,200,200,350};
  TString prefix="";
  if(getgenlevel) prefix="gen";
  TProfile2D* profilePt=(TProfile2D*)(*histfunction)(dirname,prefix+"profilePt"+suffix,"Hists2D/");
  TProfile2D* profileM=(TProfile2D*)(*histfunction)(dirname,prefix+"profileM"+suffix,"Hists2D/");
  TGraphErrors* graph=new TGraphErrors;
  for(int im=0;im<sizeof(massrange)/sizeof(double)/2;im++){
    graph->SetPoint(im,pow(TProfile2DGetMean(profileM,massrange[2*im],massrange[2*im+1],0,100),2),TProfile2DGetMean(profilePt,massrange[2*im],massrange[2*im+1],0,100));
    graph->SetPointError(im,2*TProfile2DGetMeanError(profileM,massrange[2*im],massrange[2*im+1],0,100),TProfile2DGetMeanError(profilePt,massrange[2*im],massrange[2*im+1],0,100));
  }
  return graph;
}

TGraphErrors* get_graph_data_corrected(TString dirname,TString suffix=""){
  TGraphErrors *mc_reco=get_graph(dirname,get_hist_mc_dy,0,suffix);
  TGraphErrors *mc_gen=get_graph(dirname,get_hist_mc_dy,1,"");
  TGraphErrors *data_reco=get_graph(dirname,get_hist_data_bgsub,0,suffix);
  TGraphErrors *graph=new TGraphErrors;
  
  for(int i=0;i<mc_reco->GetN();i++){
    double x=data_reco->GetX()[i]*mc_gen->GetX()[i]/mc_reco->GetX()[i];
    double y=data_reco->GetY()[i]*mc_gen->GetY()[i]/mc_reco->GetY()[i];
    graph->SetPoint(i,x,y);
    graph->SetPointError(i,sqrt(pow(data_reco->GetErrorX(i)*mc_gen->GetX()[i]/mc_reco->GetX()[i],2)+pow(data_reco->GetX()[i]*mc_gen->GetErrorX(i)/mc_reco->GetX()[i],2)+pow(data_reco->GetX()[i]*mc_gen->GetX()[i]/mc_reco->GetX()[i]/mc_reco->GetX()[i]*mc_reco->GetErrorX(i),2)),data_reco->GetErrorY(i));
  }
  return graph;
}
TGraphErrors* get_graph_old(TString dirname,void *histfunction,bool getgenlevel=0){
  const int immax=5;
  TGraphErrors* graph=new TGraphErrors;
  TString prefix="";
  if(getgenlevel) prefix="gen";
  for(int im=0;im<immax;im++){
    TH1D* pthist=(TH1D*)(*histfunction)(dirname,Form("%sdipt_m%d",prefix.Data(),im));
    TH1D* mhist=(TH1D*)(*histfunction)(dirname,Form("%sdimass_m%d",prefix.Data(),im));
    graph->SetPoint(im,mhist->GetMean()*mhist->GetMean(),pthist->GetMean());
    graph->SetPointError(im,2*mhist->GetMeanError(),pthist->GetMeanError());
  }
  return graph;
}
TGraphErrors* get_graph_data_corrected_old(TString dirname){
  TGraphErrors *mc_reco=get_graph_old(dirname,get_hist_mc_dy);
  TGraphErrors *mc_gen=get_graph_old(dirname,get_hist_mc_dy,1);
  TGraphErrors *data_reco=get_graph_old(dirname,get_hist_data_bgsub);
  TGraphErrors *graph=new TGraphErrors;
  
  for(int i=0;i<mc_reco->GetN();i++){
    double x=data_reco->GetX()[i]*mc_gen->GetX()[i]/mc_reco->GetX()[i];
    double y=data_reco->GetY()[i]*mc_gen->GetY()[i]/mc_reco->GetY()[i];
    graph->SetPoint(i,x,y);
    graph->SetPointError(i,sqrt(pow(data_reco->GetErrorX(i)*mc_gen->GetX()[i]/mc_reco->GetX()[i],2)+pow(data_reco->GetX()[i]*mc_gen->GetErrorX(i)/mc_reco->GetX()[i],2)+pow(data_reco->GetX()[i]*mc_gen->GetX()[i]/mc_reco->GetX()[i]/mc_reco->GetX()[i]*mc_reco->GetErrorX(i),2)),data_reco->GetErrorY(i));
  }
  return graph;
}
////////////////////////////////////////////////////////////
////////////////TCanvas plotting functions/////////////////
//////////////////////////////////////////////////////////
TCanvas* plot_compare_stack(TString dirname, TString histname, TString histdirectory="Hists/"){
  TCanvas* c=new TCanvas;
  TSystemDirectory dir(dirname,dirname);
  THStack* stack=new THStack(histname,histname);
  vector<TH1D*> hists;
  TH1D* hdata;
  TString stream="DoubleMuon";
  if(histname.Contains("SingleMuon")) stream="SingleMuon";

  //loop for files in directory
  for(int i=0;i<dir.GetListOfFiles()->GetSize();i++){
    TSystemFile *file=(TSystemFile*)dir.GetListOfFiles()->At(i);
    if(file->IsDirectory()) continue;
    bool isdata=0;
    if(strstr(file->GetName(),"data")!=NULL) isdata=1;
    if(isdata&&strstr(file->GetName(),stream.Data())==NULL) continue;
    TFile tfile(TString(file->GetTitle())+TString(file->GetName()));  
    TH1D* hist=(TH1D*)tfile.Get(histdirectory+histname);
    if(hist){
      hist->GetXaxis()->SetTitle(hist->GetName());
      TString title=file->GetName();
      TString subtitle=title(title.Index("_")+1,title.Index("cat")-title.Index("_")-2);
      hist->SetTitle(subtitle);
      hist->SetDirectory(0);
      if(isdata) hdata=hist;
      else if(strstr(file->GetName(),"DY")!=NULL){
	bool isexist=false;
	for(int j=0;j<hists.size();j++){
	  if(strstr(hists.at(j)->GetTitle(),"DY")!=NULL){
	    isexist=true;
	    hists.at(j)->Add(hist);
	    break;
	  }
	}
	if(!isexist) hists.push_back(hist);
      }else hists.push_back(hist);
    }
    
    TH1D* tauhist=(TH1D*)tfile.Get(histdirectory+"tau_"+histname);
    if(tauhist){
      tauhist->GetXaxis()->SetTitle(tauhist->GetName());
      tauhist->SetTitle("tautau");
      tauhist->SetDirectory(0);
      
      bool isexist=false;
      for(int j=0;j<hists.size();j++){
	if(strcmp(hists.at(j)->GetTitle(),"tautau")==0){
	  isexist=true;
	  hists.at(j)->Add(tauhist);
	  break;
	}
      }
      if(!isexist) hists.push_back(tauhist);      
    }
  }
  for(int i=0;i<hists.size();i++){
    for(int j=i+1;j<hists.size();j++){
      if(hists.at(i)->Integral()<hists.at(j)->Integral()){
	TH1D* temp=hists.at(j);
	hists.erase(hists.begin()+j);
	hists.insert(hists.begin()+j,hists.at(i));
	hists.erase(hists.begin()+i);
	hists.insert(hists.begin()+i,temp);
      }
    }
  }
  for(int i=0;i<hists.size();i++){
    TH1D* hist=hists.at(hists.size()-i-1);
    hist->SetLineColor(hists.size()-i+1);
    hist->SetFillColor(hists.size()-i+1);
    hist->SetFillStyle(1001);
    stack->Add(hist,"HIST e");
  }
  c->Divide(1,2);
  c->cd(1);
  gPad->SetPad(0,0.3,1,1);
  stack->Draw();
  stack->GetXaxis()->SetTitle(histname);
  stack->SetMaximum(stack->GetMaximum()>hdata->GetMaximum()?stack->GetMaximum()*1.05:hdata->GetMaximum()*1.05);
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(0.4);
  hdata->Draw("same e");
  gPad->BuildLegend(0.65,0.5,0.88,0.88);
  gPad->Modified();
  
  c->cd(2);
  gPad->SetPad(0,0,1,0.3);
  gPad->SetTopMargin(0.05);
  TH1D* hmc=hists.at(0)->Clone("hmc");
  hmc->SetDirectory(0);
  for(int i=1;i<hists.size();i++){
    hmc->Add(hists.at(i));
  }
  TH1D* ratio=hdata->Clone("ratio");
  ratio->Divide(hmc);
  ratio->SetDirectory(0);
  ratio->SetStats(0);
  ratio->SetTitle("");
  ratio->Draw();
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->GetYaxis()->SetLabelSize(0.1);
  gPad->SetGridy();

  return c;
}
TCanvas* plot_final(TString dirname,TString suffix=""){
  TCanvas* c=new TCanvas;
  TGraphErrors *mc_reco=get_graph(dirname,get_hist_mc_dy,0,suffix);
  TGraphErrors *mc_gen=get_graph(dirname,get_hist_mc_dy,1,"");
  TGraphErrors *data_reco=get_graph(dirname,get_hist_data_bgsub,0,suffix);
  TGraphErrors *data_cor=get_graph_data_corrected(dirname);
  
  mc_reco->SetTitle("mc_reco");
  mc_reco->Draw();
  c->SetLogx();  
  mc_reco->GetYaxis()->SetRangeUser(10,30);
  mc_reco->GetYaxis()->SetTitle("<p_{T}(DY)>");
  mc_reco->GetXaxis()->SetRangeUser(100,500000);
  mc_reco->GetXaxis()->SetTitle("M^{2}(DY)");
  mc_reco->SetFillColor(0);
  mc_reco->SetLineWidth(2);
  mc_gen->SetTitle("mc_gen");
  mc_gen->Draw("same");
  mc_gen->SetLineColor(4);
  mc_gen->SetFillColor(0);
  mc_gen->SetLineWidth(2);
  data_reco->SetTitle("data_reco");
  data_reco->Draw("same");
  data_reco->SetLineColor(2);
  data_reco->SetFillColor(0);
  data_reco->SetLineWidth(2);
  data_cor->SetTitle("data_cor");
  data_cor->Draw("same");
  data_cor->SetLineColor(6);
  data_cor->SetFillColor(0);
  data_cor->SetLineWidth(2);
  c->BuildLegend(0.6,0.12,0.88,0.33);
  mc_reco->SetTitle("");
  return c;
}

//////////////////////////////////////////////////////////////
/////////////////Save all plots/////////////////////////////////
/////////////////////////////////////////////////////////////
void SaveAll(TString dirname,TString suffix="",TString option="",TString format="pdf"){
  TString histnames[]={"dimass","dipt","met","nvtx","l1pt","l2pt","l1eta","l2eta"};
  for(int i=0;i<sizeof(histnames)/sizeof(TString);i++){
    TCanvas *c=plot_compare_stack(dirname,histnames[i]+suffix);
    c->SaveAs(histnames[i]+suffix+"."+format);
  }
  if(option.Contains("MORE")||option.Contains("more")){
    for(int im=0;im<5;im++){
      for(int i=0;i<sizeof(histnames)/sizeof(TString);i++){
	TCanvas *c=plot_compare_stack(dirname,Form("%s_m%d%s",histnames[i].Data(),im,suffix.Data()));
	c->SaveAs(Form("%s_m%d%s.%s",histnames[i].Data(),im,suffix.Data(),format.Data()));
      }
    }
  }
  TCanvas* c=plot_final(dirname,suffix);
  c->SaveAs("graph"+suffix+"."+format);
}
