#include<vector>
//const double massrange[]={40,60,60,80,80,100,100,200,200,350};
//const double massrange[]={40,50,50,60,60,70,70,80,80,100,100,130,130,200,200,350};
//const double massrange[]={20,30,30,40,40,50,50,60,60,70,70,80,80,100,100,130,130,200,200,350};
//const double massrange[]={20,22,22,24,24,26,26,28,28,30,30,32,32,34,34,36,36,38,38,40,40,42,42,44,44,46,46,48,48,50,50,52,52,56,56,60,60,64,64,68,68,72,72,76,76,80,80,84,84,88,88,92,92,96,96,100,100,130,130,160,160,200,200,350};
const double massrange[]={18,28,28,38,38,50,50,60,60,80,80,100,100,200,200,350};
const int massbinnum=sizeof(massrange)/sizeof(double)/2;
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
  int xminbin=profile->GetXaxis()->FindBin(xmin);
  int xmaxbin=profile->GetXaxis()->FindBin(xmax);
  int yminbin=profile->GetYaxis()->FindBin(ymin);
  int ymaxbin=profile->GetYaxis()->FindBin(ymax);
  const int statnum=4;
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
TH1D* TProfile2DGetHistX(TProfile2D* profile,double ymin,double ymax){
  int yminbin=profile->GetYaxis()->FindBin(ymin);
  int ymaxbin=profile->GetYaxis()->FindBin(ymax);
  TH2D* hist2d=profile->ProjectionXY("hist2d","b");
  TH1D* hist1d=hist2d->ProjectionX("hist1d",yminbin,ymaxbin-1);
  TH1D* hist=(TH1D*)hist1d->Clone(profile->GetName()+TString("_px"));
  delete hist2d;delete hist1d;
  return hist;
}
TH1D* TProfile2DGetHistY(TProfile2D* profile,double xmin,double xmax){
  int xminbin=profile->GetXaxis()->FindBin(xmin);
  int xmaxbin=profile->GetXaxis()->FindBin(xmax);
  TH2D* hist2d=profile->ProjectionXY("hist2d","b");
  TH1D* hist1d=hist2d->ProjectionY("hist1d",xminbin,xmaxbin-1);
  TH1D* hist=(TH1D*)hist1d->Clone(profile->GetName()+TString("_py"));
  delete hist2d;delete hist1d;
  return hist;
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
TGraphErrors* TGraphMultiply(TGraphErrors* graph1,TGraphErrors* graph2){
  if(graph1->GetN()!=graph2->GetN()){
    cout<<"Error::TGraphMultiply::Different number of points"<<endl;
    return NULL;
  }
  TGraphErrors *graph=new TGraphErrors;
  for(int i=0;i<graph1->GetN();i++){
    double x=graph1->GetX()[i]*graph2->GetX()[i];
    double y=graph1->GetY()[i]*graph2->GetY()[i];
    graph->SetPoint(i,x,y);
    graph->SetPointError(i,sqrt(pow(graph1->GetErrorX(i)*graph2->GetX()[i],2)+pow(graph1->GetX()[i]*graph2->GetErrorX(i),2)),sqrt(pow(graph1->GetErrorY(i)*graph2->GetY()[i],2)+pow(graph1->GetY()[i]*graph2->GetErrorY(i),2)));
  }
  return graph;
}  
TGraphErrors* TGraphDivide(TGraphErrors* graph1,TGraphErrors* graph2){
  if(graph1->GetN()!=graph2->GetN()){
    cout<<"Error::TGraphDivide::Different number of points"<<endl;
    return NULL;
  }
  TGraphErrors *graph=new TGraphErrors;
  for(int i=0;i<graph1->GetN();i++){
    double x=graph1->GetX()[i]/graph2->GetX()[i];
    double y=graph1->GetY()[i]/graph2->GetY()[i];
    graph->SetPoint(i,x,y);
    graph->SetPointError(i,sqrt(pow(graph1->GetErrorX(i)/graph2->GetX()[i],2)+pow(graph1->GetX()[i]/graph2->GetX()[i]/graph2->GetX()[i]*graph2->GetErrorX(i),2)),sqrt(pow(graph1->GetErrorY(i)/graph2->GetY()[i],2)+pow(graph1->GetY()[i]/graph2->GetY()[i]/graph2->GetY()[i]*graph2->GetErrorY(i),2)));
  }
  return graph;
}
TH1D* TGraphGetHistX(TGraphErrors* graph){
  int npoint=graph->GetN();
  TH1D* hist=new TH1D("graph_x","graph_x",npoint,0,npoint);
  for(int i=0;i<npoint;i++){
    hist->SetBinContent(i+1,graph->GetX()[i]);
    hist->SetBinError(i+1,graph->GetErrorX(i));
  }
  return hist;
}
TH1D* TGraphGetHistY(TGraphErrors* graph){
  int npoint=graph->GetN();
  TH1D* hist=new TH1D("graph_y","graph_y",npoint,0,npoint);
  for(int i=0;i<npoint;i++){
    hist->SetBinContent(i+1,graph->GetY()[i]);
    hist->SetBinError(i+1,graph->GetErrorY(i));
  }
  return hist;
}
TGraphErrors* get_graph(TString dirname,void *histfunction,bool getgenlevel=0,TString suffix=""){
  const double diptcut=400;
  TString prefix="";
  if(getgenlevel) prefix="gen";
  TProfile2D* profilePt=(TProfile2D*)(*histfunction)(dirname,prefix+"profilePt"+suffix,"Hists2D/");
  TProfile2D* profileM=(TProfile2D*)(*histfunction)(dirname,prefix+"profileM"+suffix,"Hists2D/");
  TGraphErrors* graph=new TGraphErrors;
  for(int im=0;im<sizeof(massrange)/sizeof(double)/2;im++){
    graph->SetPoint(im,pow(TProfile2DGetMean(profileM,massrange[2*im],massrange[2*im+1],0,diptcut),2),TProfile2DGetMean(profilePt,massrange[2*im],massrange[2*im+1],0,diptcut));
    graph->SetPointError(im,2*TProfile2DGetMeanError(profileM,massrange[2*im],massrange[2*im+1],0,diptcut),TProfile2DGetMeanError(profilePt,massrange[2*im],massrange[2*im+1],0,diptcut));
  }
  return graph;
}
TGraphErrors* get_graph_data_corrected(TString dirname,TString suffix=""){
  TGraphErrors* mc_reco=get_graph(dirname,get_hist_mc_dy,0,suffix);
  TGraphErrors* mc_gen=get_graph(dirname,get_hist_mc_dy,1,"");
  TGraphErrors* data_reco=get_graph(dirname,get_hist_data_bgsub,0,suffix);
  TGraphErrors* ratio=TGraphDivide(mc_gen,mc_reco);
  TGraphErrors* graph=TGraphMultiply(data_reco,ratio);
  delete mc_reco;delete mc_gen;delete data_reco;delete ratio;
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
TCanvas* plot_compare_stack(TString dirname, TString histname, TString histdirectory="Hists/",TString axistitle="",int rebin=0,double xmin=0,double xmax=0,double logymin=0){
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
      if(rebin) hist->Rebin(rebin);
      if(xmin||xmax) hist->GetXaxis()->SetRangeUser(xmin,xmax);
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
      if(rebin) tauhist->Rebin(rebin);
      if(xmin||xmax) tauhist->GetXaxis()->SetRangeUser(xmin,xmax);
      
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
  TH1D *zm,*zt,*wjet,*diboson,*tt;
  for(int i=0;i<hists.size();i++){
    TH1D* hist=hists.at(i);
    if(strstr(hist->GetTitle(),"DY")){
      hist->SetTitle("#gamma*/Z#rightarrow#mu^{+}#mu^{-}");
      zm=hist;
    }else if(strstr(hist->GetTitle(),"tautau")){
      hist->SetTitle("#gamma*/Z#rightarrow#tau^{+}#tau^{-}");
      zt=hist;
    }else if(strstr(hist->GetTitle(),"TT")){
      hist->SetTitle("t#bar{t}");
      tt=hist;
    }else if(strstr(hist->GetTitle(),"WJet")){
      hist->SetTitle("W");
      wjet=hist;
    }else if(strstr(hist->GetTitle(),"WW")){
      hist->SetTitle("WW, WZ, ZZ");
      diboson=hist;
      cout<<"diboson="<<hist<<hist->GetTitle()<<endl;
    }else if(strstr(hist->GetTitle(),"WZ")){
      for(int j=0;j<hists.size();j++){
	if(strstr(hists.at(j)->GetTitle(),"WW")){
	  hists.at(j)->Add(hist);
	  cout<<hists.at(j)<<"+="<<hist<<hist->GetTitle()<<endl;
	}
      }
      hists.erase(hists.begin()+i);
      i--;
    }else if(strstr(hist->GetTitle(),"ZZ")){
      for(int j=0;j<hists.size();j++){
	if(strstr(hists.at(j)->GetTitle(),"WW")){
	  hists.at(j)->Add(hist);
	  cout<<hists.at(j)<<"+="<<hist<<hist->GetTitle()<<endl;
	}
      }
      hists.erase(hists.begin()+i);
      i--;
    }
  }
  /*
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
    stack->Add(hists.at(i),"HIST e");
  }
  */
  stack->Add(wjet,"HIST");wjet->SetFillColor(4);wjet->SetLineColor(4);
  stack->Add(tt,"HIST");tt->SetFillColor(kGreen-5);tt->SetLineColor(kGreen-5);
  stack->Add(diboson,"HIST");diboson->SetFillColor(kOrange+8);diboson->SetLineColor(kOrange+8);
  stack->Add(zt,"HIST");zt->SetFillColor(kOrange-5);zt->SetLineColor(kOrange-5);
  stack->Add(zm,"HIST");zm->SetFillColor(kOrange);zm->SetLineColor(kOrange);

  c->Divide(1,2);
  c->cd(1);
  gPad->SetPad(0,0.3,1,1);
  stack->Draw();
  if(xmin||xmax) stack->GetXaxis()->SetRangeUser(xmin,xmax);
  if(axistitle!="") stack->SetTitle("");
  stack->GetYaxis()->SetTitle("Events");stack->GetYaxis()->SetTitleSize(0.06);stack->GetYaxis()->SetTitleOffset(0.8);
  stack->GetYaxis()->SetLabelSize(0.05);
  stack->GetXaxis()->SetTitle("");
  stack->GetXaxis()->SetLabelSize(0);
  if(logymin){stack->SetMinimum(logymin);gPad->SetLogy();}
  stack->SetMaximum(stack->GetMaximum()>hdata->GetMaximum()?stack->GetMaximum()*1.05:hdata->GetMaximum()*1.05);
  
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(0.4);
  hdata->Draw("same e");
  TLegend* legend=new TLegend(0.65,0.5,0.88,0.88);
  legend->AddEntry(hdata,"Data","lp");
  legend->AddEntry(zm,"","f");
  legend->AddEntry(zt,"","f");
  legend->AddEntry(diboson,"","f");
  legend->AddEntry(tt,"","f");
  legend->AddEntry(wjet,"","f");
  legend->Draw();
  legend->SetBorderSize(0);
  gPad->Modified();
  
  c->cd(2);
  gPad->SetPad(0,0,1,0.35);
  gPad->SetTopMargin(0.00);
  gPad->SetBottomMargin(0.26);
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
  ratio->GetYaxis()->SetRangeUser(0.52,1.48);
  ratio->GetYaxis()->SetLabelSize(0.1);ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTitle("Data / MC");ratio->GetYaxis()->SetTitleSize(0.1);ratio->GetYaxis()->SetTitleOffset(0.35);
  ratio->GetXaxis()->SetLabelSize(0.1);
  if(axistitle!="") ratio->GetXaxis()->SetTitle(axistitle);ratio->GetXaxis()->SetTitleSize(0.1);ratio->GetXaxis()->SetTitleOffset(1);
  TLine* line=new TLine(ratio->GetXaxis()->GetBinLowEdge(ratio->GetXaxis()->GetFirst()),1,ratio->GetXaxis()->GetBinLowEdge(ratio->GetXaxis()->GetLast()+1),1); line->SetLineColor(2); line->Draw();
  gPad->SetGridy();

  return c;
}
TCanvas* plot_compare(TH1* hist1,TH1* hist2){
  TCanvas* c=new TCanvas;
  c->Divide(1,2);
  c->cd(1);
  gPad->SetPad(0,0.3,1,1);
  hist1->Draw();
  hist1->SetMaximum(hist1->GetMaximum()>hist2->GetMaximum()?hist1->GetMaximum()*1.05:hist2->GetMaximum()*1.05);

  hist2->Draw("same e");
  gPad->BuildLegend(0.55,0.5,0.88,0.88);
  gPad->Modified();
  
  c->cd(2);
  gPad->SetPad(0,0,1,0.3);
  gPad->SetTopMargin(0.05);
  TH1D* ratio=hist1->Clone("ratio");
  ratio->Divide(hist2);
  ratio->SetDirectory(0);
  ratio->SetStats(0);
  ratio->SetTitle("");
  ratio->Draw();
  ratio->GetYaxis()->SetLabelSize(0.1);
  ratio->GetYaxis()->UnZoom();
  ratio->GetXaxis()->SetLabelSize(0.1);
  gPad->SetGridy();
  return c;
}  
TCanvas* plot_final(TString dirname,TString suffix=""){
  TCanvas* c=new TCanvas;
  TGraphErrors *mc_reco=get_graph(dirname,get_hist_mc_dy,0,suffix);
  TGraphErrors *mc_gen=get_graph(dirname,get_hist_mc_dy,1,"");
  TGraphErrors *data_reco=get_graph(dirname,get_hist_data_bgsub,0,suffix);
  TGraphErrors *data_cor=get_graph_data_corrected(dirname,suffix);
  
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
  mc_reco->SetTitle("graph");
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
    for(int im=0;im<7;im++){
      for(int i=0;i<sizeof(histnames)/sizeof(TString);i++){
	TCanvas *c=plot_compare_stack(dirname,Form("%s_m%d%s",histnames[i].Data(),im,suffix.Data()));
	c->SaveAs(Form("%s_m%d%s.%s",histnames[i].Data(),im,suffix.Data(),format.Data()));
      }
    }
  }
  TCanvas* c=plot_final(dirname,suffix);
  c->SaveAs("graph"+suffix+"."+format);
}




//////////////////////////////
//////////test 20181001/////
///////////////////////////
void test(TString prefix="",bool save=0){
  TProfile2D* genprofilePt=get_hist_mc_dy(hsseopath,"genprofilePt","Hists2D/");
  TProfile2D* profilePt=get_hist_mc_dy(hsseopath,"profilePt"+prefix,"Hists2D/");
  TProfile2D* profilePt_hardFSR=get_hist_mc_dy(hsseopath,"profilePt"+prefix+"_hardFSR","Hists2D/");
  TH1D* genhistm=TProfile2DGetHistX(genprofilePt,0,100);
  TH1D* histm=TProfile2DGetHistX(profilePt,0,100);
  TH1D* histm_hardFSR=TProfile2DGetHistX(profilePt_hardFSR,0,100);
  genhistm->SetStats(0);
  genhistm->GetXaxis()->SetTitle("mass");
  genhistm->SetNameTitle("genmass","genmass");

  histm->SetLineColor(2);
  histm->SetStats(0);
  histm->GetXaxis()->SetTitle("mass");
  histm->SetNameTitle("recomass","recomass");

  histm_hardFSR->SetLineColor(3);
  histm_hardFSR->SetStats(0);
  histm_hardFSR->GetXaxis()->SetTitle("mass");
  histm_hardFSR->SetNameTitle("recomass_hardFSR","recomass_hardFSR");
  plot_compare((TH1*)histm->Clone(),genhistm);
  plot_compare(histm_hardFSR,(TH1*)histm->Clone());

  for(int i=0;i<sizeof(massrange)/sizeof(double)/2;i++){
    TH1D* genhistpt=TProfile2DGetHistY(genprofilePt,massrange[2*i],massrange[2*i+1]);
    TH1D* histpt=TProfile2DGetHistY(profilePt,massrange[2*i],massrange[2*i+1]);
    TH1D* histpt_hardFSR=TProfile2DGetHistY(profilePt_hardFSR,massrange[2*i],massrange[2*i+1]);
    TString stringmassrange=Form("[%d,%d]",(int)massrange[2*i],(int)massrange[2*i+1]);
    genhistpt->SetStats(0);
    genhistpt->GetXaxis()->SetTitle("pT");
    genhistpt->GetXaxis()->SetRangeUser(0,100);
    genhistpt->SetNameTitle("genpt"+stringmassrange,"genpt"+stringmassrange);

    histpt->SetLineColor(2);
    histpt->SetStats(0);
    histpt->GetXaxis()->SetTitle("pT");
    histpt->GetXaxis()->SetRangeUser(0,100);
    histpt->SetNameTitle("recopt"+stringmassrange,"recopt"+stringmassrange);

    histpt_hardFSR->SetLineColor(3);
    histpt_hardFSR->SetStats(0);
    histpt_hardFSR->GetXaxis()->SetTitle("pT");
    histpt_hardFSR->GetXaxis()->SetRangeUser(0,100);
    histpt_hardFSR->SetNameTitle("recopt_hardFSR"+stringmassrange,"recopt_hardFSR"+stringmassrange);
    plot_compare((TH1*)histpt->Clone(),genhistpt);
    plot_compare(histpt_hardFSR,(TH1*)histpt->Clone());
  }

  TGraphErrors* graph_reco=get_graph(hsseopath,get_hist_mc_dy,0,prefix);
  TGraphErrors* graph_gen=get_graph(hsseopath,get_hist_mc_dy,1,"");
  TGraphErrors* ratio=TGraphDivide(graph_gen,graph_reco);
  TH1D* histx=TGraphGetHistX(ratio);
  TH1D* histy=TGraphGetHistY(ratio);
  histx->SetNameTitle("Rm","Rm");
  histx->SetStats(0);
  histy->SetNameTitle("RpT","RpT");
  histy->SetStats(0);
  for(int i=0;i<histx->GetNbinsX();i++){
    histx->GetXaxis()->SetBinLabel(i+1,Form("[%d,%d]",(int)massrange[2*i],(int)massrange[2*i+1]));
    histy->GetXaxis()->SetBinLabel(i+1,Form("[%d,%d]",(int)massrange[2*i],(int)massrange[2*i+1]));
  }
  TCanvas* c2=new TCanvas;
  histx->Draw();
  TCanvas* c3=new TCanvas;
  histy->Draw();
  plot_final(hsseopath,prefix);
  if(save){
    TSeqCollection* list=gROOT->GetListOfCanvases();
    for(int i=0;i<list->GetSize();i++){
	cout<<i<<" "<<list->At(i)->GetTitle()<<endl;
	TCanvas* canvas=list->At(i);
	canvas->cd();
	canvas->Update();
	canvas->Modified();
	TPaveText* title=canvas->cd(1)->GetPrimitive("title");
	canvas->SaveAs(TString(title->GetLine(0)->GetTitle())+".png");
	canvas->cd(1)->SetLogy();
	canvas->SaveAs(TString(title->GetLine(0)->GetTitle())+"_log.png");
    }    
  }
}
  
