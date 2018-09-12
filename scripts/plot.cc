#include<vector>
TString AnalyzerName="ISR2016MuonAnalyzer";
TString mypath="/data2/CAT_SKTreeOutput/JobOutPut/hsseo/LQanalyzer/data/output/CAT/ISR2016MuonAnalyzer/periodBtoH/";

//Make comparison plot between data and MC
TCanvas* plot_compare_stack(TString dirname, TString histname, TString histdirectory="Hists/"){
  TCanvas* c=new TCanvas;
  TSystemDirectory dir(dirname,dirname);
  THStack* stack=new THStack(histname,histname);
  vector<TH1D*> hists;
  TH1D* hdata;

  //loop for files in directory
  for(int i=0;i<dir.GetListOfFiles()->GetSize();i++){
    TSystemFile *file=(TSystemFile*)dir.GetListOfFiles()->At(i);
    if(file->IsDirectory()) continue;
    bool isdata=0;
    if(strstr(file->GetName(),"data")!=NULL) isdata=1;
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
TH1D* get_hist_data_bgsub(TString dirname, TString histname){
  TSystemDirectory dir(dirname,dirname);
  vector<TH1D*> hists;
  TH1D* hdata;
  for(int i=0;i<dir.GetListOfFiles()->GetSize();i++){
    TSystemFile *file=(TSystemFile*)dir.GetListOfFiles()->At(i);
    if(file->IsDirectory()) continue;
    bool isdata=0;
    if(strstr(file->GetName(),"data")!=NULL) isdata=1;
    TFile tfile(TString(file->GetTitle())+TString(file->GetName()));
    TH1D* hist=(TH1D*)tfile.Get("Hists/"+histname);
    if(!hist) continue;
    hist->GetXaxis()->SetTitle(hist->GetName());
    TString title=file->GetName();
    TString subtitle=title(title.Index("_")+1,title.Index("cat")-title.Index("_")-2);
    hist->SetTitle(subtitle);
    hist->SetDirectory(0);
    if(isdata) hdata=hist;
    else if(strstr(file->GetName(),"DY")!=NULL){
      TH1D* tauhist=(TH1D*)tfile.Get("Hists/tau_"+histname);
      if(!tauhist) continue;
      tauhist->GetXaxis()->SetTitle(hist->GetName());
      tauhist->SetTitle("tautau");
      tauhist->SetDirectory(0);
      hists.push_back(tauhist);
    }
    else hists.push_back(hist);
  }
  for(int i=0;i<hists.size();i++){
    TH1D* hist=hists.at(i);
    hdata->Add(hist,-1.);
  }
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(0.6);
  return hdata;
}
TH1D* get_hist_mc_dy(TString dirname, TString histname){
  TSystemDirectory dir(dirname,dirname);
  TH1D* hdy;
  for(int i=0;i<dir.GetListOfFiles()->GetSize();i++){
    TSystemFile *file=(TSystemFile*)dir.GetListOfFiles()->At(i);
    if(file->IsDirectory()) continue;
    if(strstr(file->GetName(),"DY")!=NULL){
      TFile tfile(TString(file->GetTitle())+TString(file->GetName()));
      TH1D* hist=(TH1D*)tfile.Get("Hists/"+histname);
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
TGraphErrors* get_graph(TString dirname,void *histfunction,bool getgenlevel=0){
  const int immax=5;
  TGraphErrors* graph=new TGraphErrors;
  TString prefix="";
  if(getgenlevel) prefix="gen";
  for(int im=0;im<immax;im++){
    TH1D* pthist=(*histfunction)(dirname,Form("%sdipt_m%d",prefix.Data(),im));
    TH1D* mhist=(*histfunction)(dirname,Form("%sdimass_m%d",prefix.Data(),im));
    graph->SetPoint(im,mhist->GetMean()*mhist->GetMean(),pthist->GetMean());
    graph->SetPointError(im,2*mhist->GetMeanError(),pthist->GetMeanError());
  }
  return graph;
}
TGraphErrors* get_graph_data_corrected(TString dirname){
  TGraphErrors *mc_reco=get_graph(dirname,get_hist_mc_dy);
  TGraphErrors *mc_gen=get_graph(dirname,get_hist_mc_dy,1);
  TGraphErrors *data_reco=get_graph(dirname,get_hist_data_bgsub);
  TGraphErrors *graph=new TGraphErrors;
  
  for(int i=0;i<mc_reco->GetN();i++){
    double x=data_reco->GetX()[i]*mc_gen->GetX()[i]/mc_reco->GetX()[i];
    double y=data_reco->GetY()[i]*mc_gen->GetY()[i]/mc_reco->GetY()[i];
    graph->SetPoint(i,x,y);
    graph->SetPointError(i,sqrt(pow(data_reco->GetErrorX(i)*mc_gen->GetX()[i]/mc_reco->GetX()[i],2)+pow(data_reco->GetX()[i]*mc_gen->GetErrorX(i)/mc_reco->GetX()[i],2)+pow(data_reco->GetX()[i]*mc_gen->GetX()[i]/mc_reco->GetX()[i]/mc_reco->GetX()[i]*mc_reco->GetErrorX(i),2)),data_reco->GetErrorY(i));
  }
  return graph;
}

TCanvas* plot_final(TString dirname){
  TCanvas* c=new TCanvas;
  TGraphErrors *mc_reco=get_graph(dirname,get_hist_mc_dy);
  TGraphErrors *mc_gen=get_graph(dirname,get_hist_mc_dy,1);
  TGraphErrors *data_reco=get_graph(dirname,get_hist_data_bgsub);
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

void SaveAll(TString dirname,TString format="pdf",TString option=""){
  TString histnames[]={"dimass","dipt","met","nvtx","l1pt","l2pt","l1eta","l2eta"};
  for(int i=0;i<sizeof(histnames)/sizeof(TString);i++){
    TCanvas *c=plot_compare_stack(dirname,histnames[i]);
    c->SaveAs(histnames[i]+"."+format);
  }
  if(option.Contains("MORE")||option.Contains("more")){
    for(int im=0;im<5;im++){
      for(int i=0;i<sizeof(histnames)/sizeof(TString);i++){
	TCanvas *c=plot_compare_stack(dirname,Form("%s_m%d",histnames[i].Data(),im));
	c->SaveAs(Form("%s_m%d.%s",histnames[i].Data(),im,format.Data()));
      }
    }
  }
  TCanvas* c=plot_final(dirname);
  c->SaveAs("graph."+format);
}
