#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom.h"

const Int_t NDIM = 3;
const Int_t NPROP = 3;
const Char_t* PROP_NAME[] = {"E","Njets", "H"};

struct particle{
    

};

/*
 * structure:
 * tree -> PROP_NAME[0]
        -> PROP_NAME[1]
        ...
 * 
 * initially coords will be based on properties (E, phi,  # jets...)
 * but use evolution to pair signal and noise to similar areas
 * then look at how those coordinates are based on properties
*/
void Initialize(TTree*);
void evolve();

int main(){
    evolve();
    return 0;
}

void
evolve(){
    //TFile;
    TTree* tree = new TTree("tree", "evolution tree");
    Initialize(tree);
    tree->BuildIndex("id");
    Display(tree);
}

Float_t
Distance(TTree* tree, Int_t id1, Int_t id2){
    Float_t d2=0;
    Float_t* p1 = new Float_t[NDIM];
    Float_t* p2 = new Float_t[NDIM];
    
    tree->SetBranchAddress("x", p1);
    tree->GetEntryWithIndex(id1);
    tree->SetBranchAddress("x", p2);
    tree->GetEntryWithIndex(id2);
    
    for(Int_t i=0; i<NDIM; i++){
        d2 += pow(p1[i] - p2[i], 2);
    }

    delete p1;
    delete p2;
    
    return sqrt( d2 );
}

Float_t
Calculate_rating(TTree* tree){
    Int_t Nentries = tree->GetEntries();

    Float_t d;
    Float_t ratodist, oneodist;
    Int_t rating;
    
    //tree->SetBranchAddress("U", &U);
    tree->SetBranchAddress("rating", &rating);
    
    for(Int_t i=0; i<Nentries; i++){
        for(Int_t j=0; j<Nentries; j++){
            if(i != j){
                tree->GetEntryWithIndex(i);
                d = Distance(tree, id, i);
                ratodist += rating/d;
                oneodist += 1.0/d;
            }
        }
    }
    return rating =  ratodist/oneodist;
    //This is not what I want

}
/*
Float_t
Potential(TTree* tree, Int_t id){
    Float_t U;
    tree->SetBranchAddress("U", &U);
    tree->GetEntryWithIndex(id);
    return U;
}
*/
void
Add_point(TTree* tree, Int_t id, Float_t* prop, Float_t* x){
    //set branches
    //x0,x1,x2...
    Int_t rating;
    Float_t U;    
    tree->SetBranchAddress("id",     &id    );
    tree->SetBranchAddress("rating", &rating);
    tree->SetBranchAddress("U",      &U     );
    for(Int_t i=0; i<NPROP; i++){
        tree->SetBranchAddress(PROP_NAME[i], &prop[i]); //error here
    }
    tree->SetBranchAddress("x", x);

    //U =  
    
    tree->Fill();
    //Iterate();
}

void
Iterate(){
}

void
Step(){
    
}

void
Display(TTree* tree){
    Char_t* field = "";
    for(Int_t i=0; i<NDIM; i++){
        if(i==0) field  = Form("x[%i]", i);
        else     field  = Form("%s:x[%i]",field,i);
    }
    printf("%s\n",field);
    tree->Draw(field);
}

void
Initialize(TTree* tree){
    //give random position to everyone
    //or fill by assigning signal and bkgrd status
    Int_t id=-1, rating=-1;
    Float_t* x = new Float_t[NDIM];
    Float_t* prop = new Float_t[NPROP];
    Float_t U=0;
    
    tree->Branch("id",     &id,     "id/I");
    tree->Branch("rating", &rating, "rating/I");
    tree->Branch("U",      &U,      "U/F");
    for(Int_t i=0; i<NPROP; i++){
        tree->Branch(PROP_NAME[i], &prop[i], Form("%s/F",PROP_NAME[i]));
    }
    tree->Branch("x", x, Form("x[%i]/F",NDIM));
    
    //get value from somewhere
    for(Int_t i=0; i<25; i++){
        id = i;
        for(Int_t j=0; j<NPROP; j++){
            prop[j] = gRandom->Gaus(0,.02);
        }
        for(Int_t j=0; j<NDIM; j++){
            x[j]=prop[j];
        }
        //gRandom->Rndm() < 0.5) charge = -1;
        Add_point(tree, id, prop, x);
        //tree->Fill();
    }
  
    //Iterate();
    delete x;
    delete prop;
}


