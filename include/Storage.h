#ifndef STORAGE_H
#define STORAGE_H

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Storage class                                                         //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#include "TObject.h"
#include "TBrowser.h"
#include "TMath.h"
#include "TH1I.h"
class Storage : public TObject {

private:

public:
    Int_t fnmovies;
    Int_t fnusers;
    Int_t fnentries;

    Int_t fndim;
    Int_t fstep_mode;
    Int_t fniterations;
    
    Float_t fminx;
    Float_t fmaxx;
    
    Float_t fstep_size;
    Float_t fmaxdist;
    Float_t frscale;

    TH1I* huserlist;
    TH1I* hnrec;
    
   Storage();// : fnmovies(0), fnusers(0), fnentries(0) { }
   Storage(Int_t, Int_t, Int_t, Int_t);//: fnmovies(nm), fnusers(nu) { }
   ~Storage();
   virtual void Print(Option_t *option = "") const;
   void Set_Particlenums(Int_t, Int_t);
   void Add_Userlist(TH1I*);
   Bool_t IsFolder() const { return kTRUE; }
   void Browse(TBrowser* b) {
       //b->Add(&fndim);
       //b->Add(&fniterations);
       //if (!fDoubleWrap) fDoubleWrap = new TParameter<double>("fDouble", fDouble);
       //if (!fIntWrap)    fIntWrap    = new TParameter<int>("fInt", fInt);
       //b->Add(fDoubleWrap);
       //b->Add(fIntWrap);
   }
   ClassDef(Storage,1)  //Storage class
};

#endif

