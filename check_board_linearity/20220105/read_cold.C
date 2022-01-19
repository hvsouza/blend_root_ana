

void read_cold(){
  TFile *f = new TFile("cold_data.root","RECREATE");
  TTree *t1 = new TTree("t1","t1");
  
  Double_t vnom,vin,vout,vinstd,voutstd,areain,areainstd,areaout,areaoutstd,offset,offsetstd,dataset,basein,baseinstd,baseout,baseoutstd,baselineoffset;
  
  t1->Branch("vnom",&vnom,"vnon/D");
  t1->Branch("vin",&vin,"vin/D");
  t1->Branch("vout",&vout,"vout/D");
  t1->Branch("vinstd",&vinstd,"vinstd/D");
  t1->Branch("voutstd",&voutstd,"voutstd/D");
  t1->Branch("areain",&areain,"areain/D");
  t1->Branch("areainstd",&areainstd,"areainstd/D");
  t1->Branch("areaout",&areaout,"areaout/D");
  t1->Branch("areaoutstd",&areaoutstd,"areaoutstd/D");
  t1->Branch("offset",&offset,"offset/D");
  t1->Branch("offsetstd",&offsetstd,"offsetstd/D");
  t1->Branch("dataset",&dataset,"dataset/D");
  t1->Branch("basein",&basein,"basein/D");
  t1->Branch("baseinstd",&baseinstd,"baseinstd/D");
  t1->Branch("baseout",&baseout,"baseout/D");
  t1->Branch("baseoutstd",&baseoutstd,"baseoutstd/D");
  t1->Branch("baselineoffset",&baselineoffset,"baselineoffset/D");
  
  ifstream fin; 
  fin.open("20220105_cold_measurements.txt");
  
  if(fin.good() && fin.is_open()){ // Ok
  }
  else{
    cout << "ERROR reading file" << endl;
    return;
  }
  string headers;
  getline(fin,headers);
  while(!fin.fail()){
    vnom = 0; 
    vin = 0; 
    vout = 0; 
    vinstd = 0; 
    voutstd = 0; 
    areain = 0; 
    areainstd = 0; 
    areaout = 0; 
    areaoutstd = 0; 
    offset = 0; 
    offsetstd = 0;
    dataset = 0;
    baseout = 0;
    baseoutstd = 0;
    baselineoffset = 0;
       
    fin >> vnom >> vin >> vinstd >> areain >> areainstd >> vout >> voutstd >> areaout >> areaoutstd >> offset >> offsetstd >> dataset >> basein >> baseinstd >> baseout >> baseoutstd >> baselineoffset;
//     cout << vnom << " " <<  vin << " " <<  vout << " " <<  vinstd << " " <<  voutstd << " " <<  areain << " " <<  areainstd << " " <<  areaout << " " <<  areaoutstd << " " <<  offset << " " <<  offsetstd << " " <<  dataset << endl;
    if(fin.bad() || fin.fail()){
      break;
    }
    
    t1->Fill();
  
  }
  
  f->Write();
  
  t1->Draw("(areaout-baseout):(areain-basein)","","goff");
  TGraph *g = new TGraph(t1->GetEntries(),t1->GetV2(),t1->GetV1());
  g->Draw("AP");
  g->SetMarkerStyle(21);
  g->SetMarkerSize(0.7);
  
  
  
  
  
  
}
