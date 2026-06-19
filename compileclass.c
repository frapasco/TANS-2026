void Compileclass (TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){opt="kfg";} else {opt ="kg"; }
  
  gSystem->SetBuildDir("./build",true);
  gSystem->AddIncludePath("-I./include");
  gSystem->CompileMacro("./src/Vertex.cxx", opt.Data());
  gSystem->CompileMacro("./src/Track.cxx", opt.Data()); 
  gSystem->CompileMacro("./src/Intpoint.cxx", opt.Data()); 
  gSystem->CompileMacro("./macros/simulazione.cxx", opt.Data());
  gSystem->CompileMacro ("./macros/ricostruzione.cxx", opt.Data()); 
  gSystem->CompileMacro ("./macros/analysis.cxx", opt.Data());
}

void Clean(TString myopt="build"){
  gSystem->Exec("rm -rf ./build");
  if(myopt.Contains("all")){
    gSystem->Exec("rm -rf ./data.txt");
    gSystem->Exec("rm -rf ./simulazione.root");
    gSystem->Exec("rm -rf ./ricostruzione.root");
    gSystem->Exec("rm -rf ./istogrammi.root");
  }
}
