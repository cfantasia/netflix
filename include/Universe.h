#include "Particle.h"
class Universe{
  public:
    Particle* users;
    Particle* movies;
    Info* recommend;
    
    Int_t nmovies;
    Int_t nusers;
    Int_t nentries;

    Int_t ndim;
    Int_t step_mode;
    //Int_t NPROP;
    //Char_t* PROP_NAME[NPROP] = {"E","Njets", "H"};
    
    Float_t MINX;
    Float_t MAXX;
    
    Float_t STEP_SIZE;
    Float_t RTHRESH;
    Float_t MAXDIST;
    Float_t RSCALE;
    Float_t ALPHA;
    Float_t BETA;
    
    
    Char_t output_base[256];    
    Char_t pos_filename[256];
    Char_t rmse_filename[256];
    Char_t ufile_filename[256];
    Char_t loopstats_filename[256];
    Char_t movefile_filename[256];
    Char_t ssizefile_filename[256];
    Char_t movettl_filename[256];
    Char_t recdist_filename[256]; 
    Char_t reg_filename[256];
    Char_t universe_filename[256];
    Char_t uplot_filename[256];
    Char_t compare_filename[256];

    
    Char_t* user_filename;
    Char_t* predict_filename;
    Char_t* rating_filename;
    Char_t* numpredict_filename;
    Char_t* stats_filename;
    Char_t* realdist_filename;
    
    Char_t* probe_filename;
    Char_t* qual_filename;
   
    //Universe();
    Universe(Int_t, Int_t, Int_t, Int_t);
    //Universe(Int_t, Int_t, Float_t, Float_t Int_t, Float_t, Bool_t);
    ~Universe();
    
    void Update_Filenames(); 
    void Set_Prop(Int_t, Int_t, Int_t, Int_t);
    void Set_Stepsize(Float_t);
  private:
    static const Int_t NMOVIES_GIVEN=17770;
    static const Int_t NUSERS_GIVEN =480189;
   };

Universe::Universe(Int_t num_m=0, Int_t num_u=0, Int_t mode=2, Int_t num_dim=3){
    nmovies = num_m;
    nusers = num_u;
    nentries = nmovies+nusers;
     
    step_mode = mode;
    ndim = num_dim;
    
    //NPROP = 3;
    //PROP_NAME[NPROP] = {"E","Njets", "H"};
    
    MINX = -10.0;
    MAXX = 10.0;
    STEP_SIZE = 0.5;
    MAXDIST = sqrt(ndim)*(MAXX-MINX)/2;

    RTHRESH=3.604;//3.0
    RSCALE = 1.0;

    if(nusers > 0){
        users = new Particle[nusers];
        for(Int_t i=0; i<nusers; i++) users[i].Set_Ndim(ndim);
    }else{
        users = NULL;
    }
    
    if(nmovies > 0){
        movies = new Particle[nmovies];
        for(Int_t i=0; i<nmovies; i++) movies[i].Set_Ndim(ndim);
        recommend = new Info[nmovies];
    }else{
        movies = NULL;
        recommend = NULL;
    }
     
    Update_Filenames();

    user_filename =("/home/fantasia/work/netflix/output/user_list.txt");
    predict_filename=("/home/fantasia/work/netflix/output/my_predictions.txt");
    rating_filename =("/home/fantasia/work/netflix/output/real_ratings.txt");
    numpredict_filename=("/home/fantasia/work/netflix/output/numpredict.txt");
    stats_filename =("/home/fantasia/work/netflix/output/stats.txt");
    realdist_filename =("/home/fantasia/work/netflix/output/all_real_dist.ps");
    
    probe_filename=("/home/fantasia/work/netflix/data/probe.txt");
    qual_filename =("/home/fantasia/work/netflix/data/qualifying.txt");
}

Universe::~Universe(){
    if(users) delete[] users;
    if(movies) delete[] movies;
    if(recommend) delete[] recommend;
}

void
Universe::Update_Filenames(){
    sprintf(output_base,       "/home/fantasia/work/netflix/output/mode_%i_%02i", step_mode, ndim);
    sprintf(pos_filename,      "%s/positions.txt",     output_base);
    sprintf(rmse_filename,     "%s/rmse.txt",          output_base);
    sprintf(ufile_filename,    "%s/ufile.txt",         output_base);
    sprintf(loopstats_filename,"%s/loopstats.txt",     output_base);
    sprintf(movefile_filename, "%s/movefile.txt",      output_base);
    sprintf(ssizefile_filename,"%s/ssizefile.txt",     output_base); 
    sprintf(movettl_filename,  "%s/movettl.ps",        output_base);
    sprintf(recdist_filename,  "%s/recommend_dist.ps", output_base);
    sprintf(reg_filename,      "%s/regression.ps",     output_base);
    sprintf(universe_filename, "%s/universe.ps",       output_base);
    sprintf(uplot_filename,    "%s/uplot.ps",          output_base);
    sprintf(compare_filename,    "%s/compare.ps",          output_base);
}

void
Universe::Set_Prop(Int_t num_m, Int_t num_u, Int_t mode, Int_t num_dim){
    nmovies = num_m;
    nusers = num_u;
    nentries = nmovies+nusers;

    step_mode = mode;
    ndim = num_dim;

    //NPROP = 3;
    //PROP_NAME[NPROP] = {"E","Njets", "H"};
    
    MAXDIST = sqrt(ndim)*(MAXX-MINX)/2;
    
    if(users) delete[] users;
    if(nusers > 0){
        users = new Particle[nusers];
        for(Int_t i=0; i<nusers; i++) users[i].Set_Ndim(ndim);
    }else{
        users = NULL;
    }
    
    if(movies) delete[] movies;
    if(nmovies > 0){
        movies = new Particle[nmovies];
        for(Int_t i=0; i<nmovies; i++) movies[i].Set_Ndim(ndim);
        recommend = new Info[nmovies];
    }else{
        movies = NULL;
    }
    
    if(recommend) delete[] recommend;
    if(nmovies > 0){
        recommend = new Info[nmovies];
    }
    
    Update_Filenames();
}

void
Universe::Set_Stepsize(Float_t value){
   STEP_SIZE = value;
   printf("New STEP_SIZE = %.4f\n", STEP_SIZE);
}
