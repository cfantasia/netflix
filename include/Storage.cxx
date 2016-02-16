/////////////////////////////////////////////////////////////////////////////////
// Storage.cxx                                                          //
/////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "Storage.h"

ClassImp(Storage)

Storage::Storage(){
    fnmovies=0;
    fnusers=0;
    fnentries=fnmovies + fnusers;

    fndim=2;
    fstep_mode=2;
    fniterations=0;
    
    fmaxx= 10.0;
    fminx=-10.0;
    
    fstep_size=0.5;
    fmaxdist=sqrt(fndim)*(fmaxx-fminx)/2;
    frscale=1.0;

    huserlist = NULL;
}
Storage::Storage(Int_t num_m, Int_t num_u, Int_t sm, Int_t nd){
    fnmovies = num_m;
    fnusers = num_u;
    fnentries = fnmovies+fnusers;

    fndim=nd;
    fstep_mode=sm;
    fniterations=0;
    
    fmaxx= 10.0;
    fminx=-10.0;
    
    fstep_size=0.5;
    fmaxdist=sqrt(fndim)*(fmaxx-fminx)/2;
    frscale=1.0;

    huserlist = NULL;
}

Storage::~Storage() {
// Destructor
    if(huserlist) delete huserlist;
}

void
Storage::Print(Option_t *option) const{
    printf("fnmovies: %i\n", fnmovies);
    printf("fnusers: %i\n", fnusers);
    printf("fnentries: %i\n", fnentries); 

    printf("fndim: %i\n", fndim);;
    printf("fstep_mode: %i\n", fstep_mode);;
    printf("fniterations: %i\n", fniterations);
    
    printf("fminx: %.5f\n", fminx);
    printf("fmaxx: %.5f\n", fmaxx);
    
    printf("fstep_size: %.8f\n", fstep_size);
    printf("fmaxdist: %.5f\n", fmaxdist);
    printf("frscale: %.2f\n", frscale);

} 

void
Storage::Set_Particlenums(Int_t nm, Int_t nu){
    fnmovies = nm;
    fnusers = nu;
    fnentries = fnmovies + fnusers;
}

void
Storage::Add_Userlist(TH1I* hu){
    huserlist = hu;
}
