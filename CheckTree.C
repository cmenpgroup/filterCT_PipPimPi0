void CheckTree(char *inFile){

  TFile *f = new TFile(inFile);
  TTree *myTree = (TTree*)f->Get("Data");

  myTree->Print();
}
