# ISR2016MuonAnalyzer

###########################
A. SETUP
###########################

@@ It works only at tamsa1.snu.ac.kr(147.47.242.42). 

1. Get CatAnalyzer. If you already have CatAnalyzer, skip this.

   git clone https://github.com/jalmond/LQanalyzer CatAnalyzer  
   cd CatAnalyzer  
   git checkout -b CatAnalyzer_13TeV_v8-0-7.41 v8-0-7.41  



2. Get ISR2016MuonAnalyzer

   git pull --no-edit https://github.com/CMSSNU/ISR2016MuonAnalyzer CatAnalyzer_13TeV_v8-0-7.41  

You may nead to set user.name. It is OK with any string.  
    git config user.name ANYNAME



3. Add ISR2016MuonAnalyzer class at root dictionary.

   sed -i 's/#endif/#pragma link C++ class ISR2016MuonAnalyzer+;\n#endif/' LQAnalysis/Analyzers/include/LQAnalysis_LinkDef.h  

Or, you may add below line before "#endif" in LQAnalysis/Analyzers/include/LQAnalysis_LinkDef.h  
    #pragma link C++ class ISR2016MuonAnalyzer+;



4. setup CanAnalyzer

   source setup.sh

You should regist your vaid email when it is first time for setup.


###########################
B. Making ROOT files
###########################

1. Add below lines in $LQANALYZER_DIR/LQRun/txt/list_user_mc.sh

   declare -a ISRlist_DY=( 'DYJets' 'DYJets_10to50' )  
   declare -a ISRlist_background=( 'TT_powheg' 'WJets' 'WW' 'WZ' 'ZZ' )  



2. Run below commands

   sktree -a ISR2016MuonAnalyzer -S SingleMuon -list ISRlist_background -s SKTree_DiLepSkim //For 2016 full data and background MC  
   sktree -a ISR2016MuonAnalyzer -list ISRlist_DY -s FLATCAT //For Drell-Yan MC  



###########################
C. Plotting
###########################
1. execute below lines

root <<EOF  
.L scripts/plot.cc  
SaveAll("/data2/CAT_SKTreeOutput/JobOutPut/$USER/LQanalyzer/data/output/CAT/ISR2016MuonAnalyzer/periodBtoH/")  
.q  
EOF  
