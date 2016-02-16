/* Author:Cory Fantsia
 * Netflix Project
 *
 *
 * To Do:
 * Work on Step Fn
 * Movement size depends on rating; check this
 * Run for a bit
 * Identify movie vs person
 * speed up
 * make a class
 * create mega arrays rather than disk reads 1st priority
 * dynamic sizing of universe
 * add points fn
 *
 * need user index array
 * npredicts array
 * recommends array
 * NUSERS
 * if(positions) 
 * 
 */

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom.h"
#include "TPaveStats.h"
#include <fstream>

#include "consts.h"
#include "Info.h"
#include "Piece.h"
#include "Particle.h"

using namespace std;

Int_t* users = NULL;
Particle* universe=NULL;
Piece* head = NULL;
Info* recommend = NULL;
//Int_t num_predict[NMOVIES];
/* structure:
 * tree -> PROP_NAME[0]
        -> PROP_NAME[1]
        ...
 * 
 * initially coords will be based on properties (E, phi,  # jets...)
 * but use evolution to pair signal and noise to similar areas
 * then look at how those coordinates are based on properties
 */

Int_t Initialize();
void Display();
void Display(TTree*);
void netflix(Int_t, Int_t, Int_t, Bool_t, Int_t);
void Fill_tree(Particle*, TTree*);
void Step(Int_t, Int_t, Int_t);
void Iterate();

void Load_Stats();
void Write_Stats();
Int_t Load_Positions();
Int_t Write_Positions();

Int_t Index_Lookup(Int_t);
Int_t Index_Lookup(Piece*, Int_t*);
//void Make_Userlist();
Int_t Add_Piece(Int_t);
Int_t Make_Userlist();

void Write_Userlist();
Int_t Load_Userlist();
void Save_List(FILE*, Piece*);

Int_t Size_List(Piece*);
void Fill_Array(Piece*, Int_t*);

void Write_Predictions(Bool_t, Bool_t);
Float_t Compare_Predictions();
Float_t Make_Prediction(Int_t, Int_t);
Int_t Real_Rating(Int_t, Int_t);
void Optimize_RSCALE();
Float_t Calc_U();

void Test();

Int_t Recommend_Tables();
//Int_t Load_numpredict();

int main(int argc, char ** argv){
    if(argc < 2 || argc > 6 || argv[1][1]=='?'){
		fprintf(stderr,"%s usage: %s nloops [step_mode] [ndim] [predict=0] [reset]\n",argv[0],argv[0]);
        fprintf(stderr, "  step_mode=-1: Loop over 0<STEP_MODE<4\n");
        fprintf(stderr, "  ndim     =-1: Loop over 0<NDIM<15\n");
        exit( 1 );
	}

    Int_t input_default[5]={0,2,3,0,0};
    
    if(atoi(argv[1]) == -1) Test();
    else{
        if      (argc == 2) netflix(atoi(argv[1]), input_default[1], input_default[2], input_default[3], input_default[4]);
        else if (argc == 3) netflix(atoi(argv[1]), atoi(argv[2])   , input_default[2], input_default[3], input_default[4]);
        else if (argc == 4) netflix(atoi(argv[1]), atoi(argv[2])   , atoi(argv[3])   , input_default[3], input_default[4]);
        else if (argc == 5) netflix(atoi(argv[1]), atoi(argv[2])   , atoi(argv[3])   , atoi(argv[4])   , input_default[4]);
        else if (argc == 6) netflix(atoi(argv[1]), atoi(argv[2])   , atoi(argv[3])   , atoi(argv[4])   , atoi(argv[5]));
    }
    return 0;
}

void
netflix(Int_t nloops, Int_t step_mode, Int_t ndim, Bool_t predict, Int_t reset){
    if(ndim == -1 && step_mode == -1){
        for(Int_t dim=1; dim<15; dim++){
            for(Int_t mode=0; mode<3; mode++){
                netflix(nloops, mode, dim, predict, reset);
            }
        }
    }else if(ndim == -1){
        for(Int_t dim=1; dim<15; dim++){
            netflix(nloops, step_mode, dim, predict, reset);
        }
    }else if(step_mode == -1){
        for(Int_t mode=0; mode<3; mode++){
            netflix(nloops, mode, ndim, predict, reset);
        }
    }

    STEP_MODE = step_mode;
    NDIM = ndim;
    
    Int_t list_size, universe_size;
    time_t start, end;
    time_t start_iterate, end_iterate;
    time_t start_fn, end_fn;
    Double_t dif, dif_iterate;

    time(&start);
    printf("Going to iterate %i times in mode %i in %i dim.\n", nloops, STEP_MODE, NDIM);

    ifstream statsfile (stats_filename);
    if(reset==2 || !statsfile){
        statsfile.close();
        NMOVIES = 17770;
        NUSERS = Make_Userlist();
        NENTRIES = NMOVIES + NUSERS;
        Write_Stats();
        Write_Userlist();
    }else{
        statsfile.close();
        Load_Stats();
        ifstream userfile (user_filename);
        if(!userfile){ 
            userfile.close();
            printf("%s does not exist.\n", user_filename);
            list_size = Make_Userlist();
            Write_Userlist();
        }else{
            userfile.close();
            printf("%s exists.\n", user_filename);
            list_size = Load_Userlist();
        }  
        if(list_size != NUSERS){
            printf("File lengths not equal!!!!\n");
            return;
        }
    }
    //users = new Int_t[NUSERS];
    recommend = new Info[NMOVIES];
    universe = new Particle[NENTRIES];
        
    Recommend_Tables();

    sprintf(pos_filename, "%s_%i_%i.txt", posbase_filename, STEP_MODE, NDIM);                                
    ifstream posfile (pos_filename);
    if(reset || !posfile){ 
        printf("%s does not exist.\n", pos_filename);
        posfile.close();
        universe_size = Initialize(); 
    }else{
        printf("%s exists.\n", pos_filename);
        posfile.close();
        universe_size = Load_Positions();
    }
    
    for(Int_t loop=0; loop<nloops; loop++){
        printf("\tBeginning Loop %i of %i\n", loop+1, nloops);
        time(&start_iterate);
        STEP_COUNT = 0;
        MOVE_TTL = 0;
        Iterate();
        time(&end_iterate);
        dif_iterate = difftime(end_iterate, start_iterate);
        printf("movie 1469 pos is ");
        for(Int_t dim=0; dim<NDIM; dim++){
            printf(" %.4lf ",universe[NUSERS+1469-1].x[dim]);
        }printf("\n");
        //if(loop%10 == 0) printf("U is %.lf\n", Calc_U());
        printf("STEP_COUNT = %i\n", STEP_COUNT);
        printf("MOVE_TTL = %.6lf\n", MOVE_TTL);
        printf("Iterate Loop %i of %i took %.2lf seconds or %.2lf hours to run.\n", loop+1, nloops, dif_iterate, dif_iterate/3600.0);
    	
        if(myu->STEP_SIZE < 0.000001){
            printf("Move size=%.7f is below threshold.  Breaking\n", myu->STEP_SIZE); 
            break;
        }
    }
    
    if(NDIM <= 3) Display();
    Write_Positions();
       
    //These functions aren't needed everytime
    if(predict){
        time(&start_fn);
        Bool_t make_real=false, make_mine=true;//defaults
        if(reset == 2){
            make_real=true;
            make_mine=true;
        }else{
            ifstream realfile (rating_filename);
            if( !realfile) make_real=true;
            realfile.close();

            ifstream minefile (predict_filename);
            if( !minefile) make_mine=true;
            minefile.close();
        }

        Write_Predictions(make_real,make_mine); 
        time(&end_fn);
        printf("Write_Predictions took %li seconds to run.\n", end_fn-start_fn);
        printf("RMSE is %.2lf BEFORE\n", Compare_Predictions());
        
        //time(&start_fn);
        Optimize_RSCALE();
        //time(&end_fn);
        //printf("Optimize_RSCALE took %li seconds.\n", end_fn-start_fn);

        printf("Going to redo Write_Predictions...\n");
        Write_Predictions(false, true); //Rewrite predictions using new RSCALE 
        printf("done.\n");
        
        //time(&start_fn);
        Float_t error = Compare_Predictions();
        //time(&end_fn);
        //printf("Compare_Predictions took %li seconds to run.\n", end_fn-start_fn);

        printf("RMSE is %.2lf\n", error);
    }
    //TTree* etree = new TTree("etree", "evolve tree");
    //Fill_tree(universe, etree);
    //etree->BuildIndex("index");
    //Display(etree);
    delete[] users;
    delete[] universe;
    delete head; //Uses recursive delete in class destructor
    delete[] recommend;
    
    time(&end);
    dif = difftime(end, start);
    printf("Program took %.2lf seconds or %.2lf hours to run.\n", dif, dif/3600.0);
}


Float_t
Distance2(Int_t index1, Int_t index2, Float_t* dist_x2=NULL){
    if(index1 != index2){
        Float_t sum2=0, d2;
        Float_t x1, x2, dx;
             
        for(Int_t dim=0; dim<NDIM; dim++){
            x1 = universe[index1].x[dim];
            x2 = universe[index2].x[dim];
            dx = TMath::Abs(x1-x2);   
            if( dx < (MAXX-MINX)/2) d2 = pow(dx,2);
            else                    d2 = pow(MAXX-MINX-dx,2);
            
            if(dist_x2 != NULL) dist_x2[dim] = d2;
            sum2 += d2;
        }
        return  sum2;
    }
    printf("Distance between same point %i is 0!\n", index1);
    return 0;
}


Float_t
Distance(Int_t index1, Int_t index2, Float_t* delta_x2=NULL){
    return sqrt( Distance2(index1, index2, delta_x2));
}

Float_t
Calc_U(){
    Info* temp;
    Int_t len;
    Float_t U=0;
    Int_t movie_index;
    for(Int_t movie=1; movie<NMOVIES+1; movie++){
        temp = recommend + (movie-1);
        len = temp->length;
        movie_index = NUSERS + (movie-1);
        for(Int_t j=0; j<len; j++){
            //if rating is < 3 equiv to same sign
            U += -1.0*(temp->rating[j]-RTHRESH)/Distance(movie_index, temp->user_index[j]);
        }
    }    
    return U;
}

void
Calc_rating(){
    /*
    Float_t d;
    Float_t ratodist, oneodist;
 
    for(Int_t i=0; i<NENTRIES; i++){
        ratodist = 0;
        oneodist = 0;
        for(Int_t j=0; j<NENTRIES; j++){
            if(i != j){
                //tree->GetEntryWithIndex(i);
                d = Distance(i, j);
                ratodist += (universe[j].rating)/d;
                oneodist += 1.0/d;
            }
        }
        universe[i].rating =  ratodist/oneodist;  
    }
    */
}

/*
Float_t
Potential(Particle* universe, Int_t index){
    Float_t U;
    tree->SetBranchAddress("U", &U);
    tree->GetEntryWithIndex(index);
    return U;
}
*/

void
Add_Point(Int_t id, Float_t* x, Bool_t is_user){
    //set branches
    //x0,x1,x2...
    /*
    Int_t rating;
    Float_t U;    
    tree->SetBranchAddress("index",     &index    );
    tree->SetBranchAddress("rating", &rating);
    tree->SetBranchAddress("U",      &U     );
    for(Int_t i=0; i<NPROP; i++){
        tree->SetBranchAddress(PROP_NAME[i], &prop[i]); //error here
    }
    tree->SetBranchAddress("x", x);

    //U =  
    
    tree->Fill();
    //Iterate();
    NENTRIES++; //update this
    */

    //users, universe, maybe recommend NUSERS, NENTRIES
    if( is_user ){
        Int_t* users2 = new Int_t[NUSERS+1];
        Int_t offset = 0;
        for(Int_t i=0; i<NUSERS+1; i++){
            if(users[i-offset] != id){
                if(i != NUSERS+1){
                    if(offset || users[i-offset] < id){
                        users2[i] = users[i-offset];
                    }else{
                        users2[i] = id;
                        offset = 1;
                    }
                }else{
                    if(offset) users2[i] = users[i-offset];
                    else       users2[i] = id;
                }
            }else{
                printf("Tried to add an existing member %i\n!", id);
                delete[] users2;
                return;
            }
        }
        delete[] users;
        users = users2;
        NUSERS++;
    }else{  // It's a movie!
        //Add entry to recommend
        Info* recommend2 = new Info[NMOVIES+1];
        for(Int_t i=0; i<NMOVIES; i++){
            
        }
        delete[] recommend;
        recommend = recommend2;
        NMOVIES++;
    }  
    NENTRIES = NUSERS + NMOVIES;
    Particle* universe2 = new Particle[NENTRIES];
    for(Int_t i=0; i<NENTRIES; i++){
        //insert new entry
    }
    delete[] universe;
    universe = universe2;
}

void
Parse(ifstream* infile, Int_t* user_id, Int_t* rating){
    //, Int_t* year, Int_t* month, Int_t* day){ 
    //1488844,3,2005-09-06
    Char_t line[15];
    
    infile->getline(line, 10, ','); //get customer name 
    *user_id = atoi(line);
    
    infile->getline(line, 10, ','); //get rating
    *rating = atoi(line);

    infile->getline(line, 15);
    /*
    infile>>*user_id;
    infile.ignore(1);
    infile>>*rating;
    infile.ignore(1);
    infile.ignore(256, '\n');
    */
    /*
    infile->getline(line, 10, '-'); //get year
    *year = atoi(line);
    
    infile->getline(line, 10, '-'); //get month
    *month = atoi(line);
    
    infile->getline(line, 10); //get day
    *day = atoi(line);
    */
}

Char_t*
Filename(Int_t movie_id){
    //There are 7 numbers in ID
    Char_t* zeros;
    if     (movie_id<10) zeros = Form("000000"); //how many 0's do I need
    else if(movie_id<100) zeros = Form("00000");
    else if(movie_id<1000) zeros = Form("0000");
    else if(movie_id<10000) zeros = Form("000");
    else if(movie_id<100000) zeros = Form("00");
    else if(movie_id<1000000) zeros = Form("0");
    else if(movie_id<10000000) zeros = Form("");
    else                        zeros = Form("!!!!!!!");

    //if(movie_id%10000==1) printf("movie_id:%i zeros:%s \n", movie_id, zeros);
    
    return Form("/home/fantasia/work/netflix/data/training_set/mv_%s%i.txt", zeros, movie_id); 
}

void
Iterate(){
    //Char_t* in_filename;
    //Char_t line[10];
    //Int_t user_id, rating;//, year, month, day;
    //Int_t user_index;
    //Int_t count=0;
    Info* temp=NULL;
    /*for(Int_t movie=NUSERS; movie<NENTRIES;movie++){
        //open file 'movie'
        in_filename = Filename(movie-NUSERS+1);
        ifstream infile (in_filename);
        
        if (infile.is_open()){
            //Cory: try ignore fn
            //infile.ignore(100, '\n');
                          
            infile.getline(line, 10); //get movie name
            if(atoi(line) != movie-NUSERS+1) printf("movie:%i file:%i\n", movie, atoi(line));
            while( !infile.eof()){
                //Parse(&infile, &user_id, &rating);//, &year, &month, &day);
                infile>>user_id;
                infile.ignore(1);
                infile>>rating;
                infile.ignore(256, '\n');
                //printf("movie:%i user:%i rating:%i\n", movie, user_id, rating);
                if(user_id != 0){
                    if(count%10000000==0) printf("file:%s count:%i\n", in_filename, count);
                    count++;
                    //user_index = Index_Lookup(head, &user_id); //Cory:Change this to array fn
                    user_index = Index_Lookup(user_id);
                    if(user_index == -1) printf("file:%s user_id: %i rating:%i\n", in_filename, user_id, rating);

                    Step(movie, user_index, rating);
                    STEP_COUNT++;
                }
            }
            
            infile.close();
        }else printf("Unable to open file %s\n", in_filename);
    }
    */
    Int_t length;
    Int_t movie_index;
    for(Int_t movie=1; movie<NMOVIES+1; movie++){
        temp = recommend + movie - 1;
        length = temp->length;
        movie_index = NUSERS + (movie - 1);
        //if( movie < 10) printf("length of movie %i is %i\n", movie, length);
        for(Int_t i=0; i<length; i++){
            Step(movie_index, temp->user_index[i], temp->rating[i]);
            STEP_COUNT++;
        }
    }
}

Int_t
Sign(Int_t index1, Int_t index2, Int_t rating, Int_t dim){
//returns the sign of the move
/*
  1 if away && index1.x > index2.x  //what if index1.x[0] > but index1.x[1] <
  1 if towards && index2.x > index2.x
  -1 else
*/
    Int_t dir=0;
    Float_t dx = universe[index1].x[dim]-universe[index2].x[dim];  
          
    if( dx > 0 ) dir =  1;
    else         dir = -1;

    if( TMath::Abs(dx) > (MAXX-MINX)/2) dir*=-1;
    else                                dir*= 1;
    
    if(rating < RTHRESH) dir *=  1;
    else                 dir *= -1;   
    
    return dir;
}

Float_t
Rating_Mag(Int_t rating){
    //Should this fn do anything with rating = 3?
    return TMath::Abs((rating-RTHRESH)/2.0);
}

Float_t
OutofBounds(Float_t pos){
    while(pos>MAXX || pos <= MINX){
        //printf("pos before:%.4lf ", pos);
        if      (pos >  MAXX){
            pos = pos - MAXX + MINX;
        }else{
            pos = pos - MINX + MAXX;
        }
        //printf("pos after:%.2lf\n", pos);
    }
    
    //if(pos > MAXX || pos <=MINX) printf("Still out of bounds pos=%.2lf\n", pos);
    return pos;
}

void
Step(Int_t index1, Int_t index2, Int_t rating){
/*should two objects move closer or further based on rating
 * move = STEP_SIZE*(sign*mag)*(dist_x2/dist2)
 * should this be devided by 1/dist2 or 1/distx2
 * 
 * STEP_MODE: 0=1/r
 *            1=random
 *            2=Hooke
 */
    Float_t sign, move, mag_move, rat_mag;
    //Float_t dist  = Distance (index1, index2);
    Float_t* dist_x2 = new Float_t[NDIM];
    Float_t dist2 = Distance2(index1, index2, dist_x2);
    Float_t scale=0;
    if      (STEP_MODE == 0) scale = (STEP_SIZE/dist2)/dist2; //
    else if (STEP_MODE == 1) scale = (STEP_SIZE/dist2); 
    else                     scale = (STEP_SIZE/dist2)*sqrt(dist2);
    
    for(Int_t dim=0; dim<NDIM; dim++){
        sign    = Sign(index1, index2, rating, dim);
        rat_mag = Rating_Mag(rating);
      
        mag_move = rat_mag*scale*dist_x2[dim];
      
        if(mag_move > MAXDIST){
            mag_move = MAXDIST;
            printf("BIGGER\n");
        }
        //mag_move = TMath::Min(mag_move, (Float_t)sqrt(dist_x2[dim])/2);

        move = sign*mag_move;
        universe[index1].x[dim] += move;
        universe[index2].x[dim] -= move;
        
        universe[index1].x[dim] = OutofBounds( universe[index1].x[dim] );
        universe[index2].x[dim] = OutofBounds( universe[index2].x[dim] );

        if(STEP_COUNT %20000000 == 0) printf("d2=%.4lf dx2=%.4lf |move| is %.4lf\n", dist2, dist_x2[dim], TMath::Abs(move)); 
        MOVE_TTL += TMath::Abs(move);
    }
    delete[] dist_x2;
}

void
Display(){
    TCanvas* canvas1;
    TH3F* husers=NULL;
    TH3F* hmovies=NULL;
    Float_t* pos;
    
    canvas1 = (TCanvas*)gROOT->FindObject("canvas1");
    if(canvas1) delete canvas1;

    //if     (NDIM == 2) grid = (TH2F*)gROOT->FindObject("grid");
    if(NDIM == 3){
        husers  = (TH3F*)gROOT->FindObject("husers");
        hmovies = (TH3F*)gROOT->FindObject("hmovies");
    }
    if(husers)  delete husers;
    if(hmovies) delete hmovies;
       
    canvas1 = new TCanvas("canvas1" "The Universe");
    //if     (NDIM == 2) grid = new TH2F("grid", "The Universe", 40, MINX, MAXX, 40, MINX, MAXX);
    if(NDIM == 3){
        husers  = new TH3F("husers", "The Universe", 40, MINX, MAXX,
                           40, MINX, MAXX, 40, MINX, MAXX);
        hmovies = new TH3F("hmovies", "The Universe", 40, MINX, MAXX,
                           40, MINX, MAXX, 40, MINX, MAXX);
    }
    husers->SetMarkerColor(kGreen);
    hmovies->SetMarkerColor(kRed);
    for(Int_t i=0; i<NENTRIES; i++){
        pos=universe[i].x;
        if(pos[0] > MAXX || pos[0] <= MINX) printf("Out of Bounds Index:%i pos[0]:%.lf\n", i, pos[0]);
        if(pos[1] > MAXX || pos[1] <= MINX) printf("Out of Bounds Index:%i pos[1]:%.lf\n", i, pos[1]);
        if(NDIM == 3) if(pos[2] > MAXX || pos[2] <= MINX) printf("Out of Bounds Index:%i pos[2]:%.lf\n", i, pos[2]);          
        if(universe[i].is_user){
            if      (NDIM == 2) husers->Fill(pos[0], pos[1]);
            else if (NDIM == 3) husers->Fill(pos[0], pos[1], pos[2]);    
        }else{
            if      (NDIM == 2) hmovies->Fill(pos[0], pos[1]);
            else if (NDIM == 3) hmovies->Fill(pos[0], pos[1], pos[2]); 
        }
    }
    husers->Draw();
    hmovies->Draw("sames");
    canvas1->SaveAs(pic_filename);

    delete husers;
    delete hmovies;
    delete canvas1;
}

            
void
Display(TTree* tree){
    Char_t* field = Form("");;
    for(Int_t i=0; i<NDIM; i++){
        if(i==0) field  = Form("x[%i]", i);
        else     field  = Form("%s:x[%i]",field,i);
    }
    printf("%s\n",field);
/*
    TH2F* htemp;
    //tree->SetMarkerSize(100.0);
    tree->SetMarkerColor(kBlue);
    tree->Draw(field, "rating>0", "");
    htemp = (TH2F*) gROOT->FindObject("htemp");
    Float_t min[NDIM];
    Float_t max[NDIM];
    for(Int_t i=0; i<NDIM; i++){
        TAxis* axis;
        if(i == 0) axis = htemp->GetXaxis();
        if(i == 1) axis = htemp->GetYaxis();
        if(i == 2) axis = htemp->GetZaxis();

        min[NDIM] = axis->GetXmin();
        max[NDIM] = axis->GetXmax();
    }
    tree->SetMarkerColor(kRed);
    tree->Draw(field, "rating<0", "same");
    
    Float_t max = tree->GetMaximum("x");
    gPad->SetX1(min,max);
    gPad->SetYaxis()->SetAxisRange(min,max);
    //gPad->Update();
*/
    tree->Draw(field, "", "");
}

Int_t
Initialize(){
    //give random position to everyone
    //or fill by assigning signal and bkgrd status

    printf("\tEntering Initialize\n");

    //Char_t line[10];
    //Int_t users[NUSERS];
    /*
    ifstream myfile (user_filename);
    if (myfile.is_open()){
        for(Int_t i=0; !myfile.eof(); i++){
            myfile.getline(line, 10); //get user id
            users[i] = atoi(line);
            //if(i!=0 && users[i-1]>=users[i]) printf("i:%d %d>=%d\n", i, users[i-1], users[i]);
            //if(i%1000==0) printf(" i=%i\n", i);
        }
        if(!myfile.eof()) printf("NUSERS:%i \n", NUSERS); //Cory: check why this is bugging out
        
       myfile.close();
0    }else printf("Unable to open file %s\n", user_filename);
    */
    
    //tree->ReadFile("cory.dat", "evt/i:t/i:c0/f:c1/f:c2/f:c3/f:c4/f:c5/f:c6/f:c7/f:c8/f:c9/f");
    Int_t obj_id; //get from file
    gRandom->SetSeed(0); //Argument = 0 means use computer clock
    printf("Current Seed is %i\n", gRandom->GetSeed());
    //Piece* temp = head;
    /*
      for(Int_t i=0; i<20; i++){
        if(temp != NULL){
            printf("head[%i]=%i\n", i, temp->id);
            temp = temp->less;
        }
    }
    */
    for(Int_t i=0; i<NENTRIES; i++){
        //don't repeat numbers
        universe[i].index = i; //give INDEX # to every1
        
        if(i<NUSERS){
            obj_id = users[i];
            universe[i].is_user=true;
        }else{
            obj_id = i+1 - NUSERS;
            universe[i].is_user=false;
        }
        universe[i].id = obj_id;
        
        for(Int_t j=0; j<NDIM; j++){
            //should be using linear
            universe[i].x[j] = gRandom->Uniform(MINX, MAXX);
            //universe[i].x[j] = gRandom->Gaus(50,10);//universe[i].prop[j];
        }
        
        //gRandom->Rndm() < 0.5) charge = -1;
        //Add_point(universe, index, prop, x);
        
    }
    //Calc_U(universe);
    //Calc_rating(universe);

    printf("\tLeaving Initialize\n");
    return NENTRIES; //return # of users assigned pos
}

void
Fill_tree(Particle* universe, TTree* tree){
/*
    Int_t index=-1;
    Float_t rating=-1;
    Float_t x[NDIM];
    Float_t prop[NPROP];
    Float_t U=0;

    tree->Branch("index",     &index,      "index/I");
    tree->Branch("rating", &rating, "rating/F");
    tree->Branch("U",      &U,      "U/F");
    
    tree->Branch("x",       x,       Form("x/F",NDIM));//Cory: started breaking it here

    for(Int_t i=0; i<NPROP; i++){
        tree->Branch(PROP_NAME[i], &prop[i], Form("%s/F",PROP_NAME[i]));
    }
    
    for(Int_t i=0; i<NENTRIES; i++){
        index = universe[i].index;
        rating = universe[i].rating;
        for(Int_t j=0; j<NPROP; j++){
            prop[j] = universe[i].prop[j];
        }
        for(Int_t j=0; j<NDIM; j++){
            x[j] = universe[i].x[j];
        }
        U = universe[i].U;
        tree->Fill();
    }
*/
}

/*
void
Make_Userlist(){
    //verify num of unique customers
    //use linked list (thanks ben!)

    
    //get value from somewhere
    //file: probe.txt
    //IMPORT FORMAT:
    //movieid:
    //customerid, rating, year-month-day
    //customerid, rating, year-month-day
    //......
    
    Int_t users_entered=0;
    Char_t line[25];
    Char_t* filename;
    //Char_t* zeros;
    Int_t id;
    Bool_t repeat;
    for(Int_t i=1; i<NUSERS+1 && users_entered<NUSERS; i++){
        filename = Filename(i);
        if(i%100==1) printf("users_entered: %i\n", users_entered);

        ifstream myfile (filename);
        if (myfile.is_open()){
            myfile.getline(line, 25); //get movie name
            //printf("1st Line: %s\n", line);//take off the colon and put into variable
            //printf("before while\n");
            while (!myfile.eof() ){
                myfile.getline (line, 25, ',');
                //printf("good part: %s\n", line);
                id=atoi(line);//truncate at "," and put in variable
                myfile.getline (line, 25);
                //printf("rest if: %s\n", line);
                repeat = false;
                for(Int_t j=0; j<users_entered;j++){
                    if(id==users[j]){
                        repeat = true;
                        break;
                    }
                }
                if(!repeat) users[users_entered++]=id;
            }
            //printf("end of while\n");
            myfile.close();
        }
        else printf("Unable to open file %s\n", filename); 
    }
    printf("Begin Sort...");
    sort(users, users+NUSERS); //sort array using standard call
    printf("Done.\n");
    //write list to a file
    FILE* outfile;
    outfile = fopen (user_filename, "w");
    for (Int_t i=0 ; i<NUSERS ; i++){
        fprintf (outfile, "%i\n",users[i]);
    }
    fclose (outfile);

}
*/
void
Write_Stats(){
    FILE* outfile;
    outfile = fopen (stats_filename, "w");
    fprintf(outfile, "NENTRIES: %i\n", NENTRIES);
    fprintf(outfile, "NMOVIES: %i\n", NMOVIES);
    fprintf(outfile, "NUSERS: %i\n", NUSERS);
    fclose (outfile);
}

void
Load_Stats(){
    ifstream statsfile(stats_filename);
    if(statsfile.is_open()){
        statsfile.ignore(256, ':');
        statsfile>>NENTRIES;
        statsfile.ignore(256, '\n');

        statsfile.ignore(256, ':');
        statsfile>>NMOVIES;
        statsfile.ignore(256, '\n');

        statsfile.ignore(256, ':');
        statsfile>>NUSERS;
        printf("NENTRIES:%i NMOVIES:%i NUSERS:%i\n", NENTRIES, NMOVIES, NUSERS);
    }else printf("Countn't open file %s\n", stats_filename);
}

Int_t
Write_Positions(){
   //Format:
    //index, id, x[0], ..., x[NDIM-1], is_user 
    Particle* temp;
    
    FILE* outfile;
    printf("About to open %s for writing with base:%s\n", pos_filename, posbase_filename);
    outfile = fopen (pos_filename, "w");
    for (Int_t i=0; i<NENTRIES; i++){
        temp = universe + i; 
        
        fprintf (outfile, "%i, %i, ",temp->index, temp->id);
        for(Int_t dim=0; dim<NDIM; dim++){
            fprintf (outfile, "%f, ",temp->x[dim]);
        }
        fprintf (outfile, "%i\n",temp->is_user);
    }
    fclose (outfile);
    return NENTRIES;
}

Int_t
Load_Positions(){
    printf("\tEntering Load_Positions.\n");
    Int_t obj_entered=0;
    ifstream infile (pos_filename);
    if (infile.is_open()){
        Particle* temp;
        while (!infile.eof() && obj_entered < NENTRIES){
            //printf("Entering while count:%i\n", count);
            //READ IN PARAMETERS
            //Format:
            //index, id, x[0], ..., x[NDIM-1], is_user 
            if(obj_entered%100000 == 0) printf("obj_entered:%i\n", obj_entered);
            temp = universe + obj_entered; 
            
            infile>>temp->index;//.getline(line, 20, ','); //get index 
            infile.ignore(256, ',');//temp->index = atoi(line);
            
            infile>>temp->id;//infile.getline(line, 20, ','); //get ID
            infile.ignore(256, ',');//temp->id = atoi(line);
            
            for(Int_t dim=0; dim<NDIM; dim++){
                infile>>temp->x[dim];//infile.getline(line, 20, ','); //get x[dim]
                infile.ignore(256, ',');//temp->x[dim] = atof(line); //needs to be a float
            }
            infile>>temp->is_user;
            infile.ignore(256, '\n');
            obj_entered++;
            //should verify that users have positions
        }
        //printf("At the end %i\n", obj_entered);
        /*
        while(!infile.eof()){
            infile>>id;
            if( id != 0 ){
                printf("user %i was not added\n", id);
            }
        }
        */
        infile.close();
    }
    else printf("Unable to open file %s\n", pos_filename);
        
    printf("\tLeaving Load Positions.\n");
    return obj_entered;
}

Int_t
Index_Lookup(Int_t user_id){
    //Given the user_id, find the index number;
    Int_t first=0;
    Int_t last= NUSERS-1;
    
    while (first <= last) {
        Int_t mid = (first + last) / 2;  // compute mid point.
        if (user_id > users[mid]) 
            first = mid + 1;  // repeat search in top half.
        else if (user_id < users[mid]) 
            last = mid - 1; // repeat search in bottom half.
        else
            return mid;     // found it. return position /////
    }
    printf("ALERT!!!! Could not find user_id=%i\n", user_id);
    return -1;    // failed to find key
}

Int_t
Index_Lookup(Piece* temp, Int_t* user_id){
/* fn is really slow.  don't use*/
    
    if(temp!= NULL){
        if(*user_id >= temp->id){
            return 1 + Index_Lookup(temp->more, user_id) + Index_Lookup(temp->less, user_id); 
        }else{
            return Index_Lookup(temp->less, user_id); 
        }
    }
    return 0; //what if it really is the first 1??
}

Int_t
Add_Piece(Int_t user_id){
    Piece* temp=head;
    //printf("temp: %x\n", temp);
    if(head != NULL){
        Piece* backup= temp;
        while(temp!=NULL){
            if      (user_id < temp->id){
                backup=temp;
                temp = temp->less;
            }else if(user_id > temp->id){
                backup=temp;
                temp = temp->more;
            }else{
                return -1; 
            }
        }
        //Fell off list
        if(user_id < backup->id){
            backup->less = new Piece(user_id);
            temp = backup->less;
        }else{
            backup->more = new Piece(user_id);
            temp = backup->more;
        }
    }else{
        temp = new Piece(user_id);
        head = temp;
        //printf("head:%x\n", head);
    }
    return 0;
}

Int_t
Make_Userlist(){
    //verify num of unique customers
    //use linked list (thanks ben!)
    
    //get value from somewhere
    //file: 
    //IMPORT FORMAT:
    //movieid:
    //customerid, rating, year-month-day
    //customerid, rating, year-month-day
    //......
    printf("\tEntering Make_Userlist.\n");
    Int_t users_entered=0, count;
    Char_t* filename;
    Int_t user_id, rating, movie;

    FILE* outfile;
    outfile = fopen (numpredict_filename, "w");
    
    for(Int_t i=1; i<NMOVIES+1; i++){
        count = 0;
        filename = Filename(i);
        if(i%1000==1) printf("file:%s users_entered: %i\n", filename, users_entered);
        
        ifstream infile (filename);
        if (infile.is_open()){
            infile>>movie;//get movie name
            infile.ignore(256, '\n');
            //printf("Read in movie:%i\n", movie);
            if(i != movie) printf("i:%i movie:%i\n", i, movie);  
            while (!infile.eof()){
                infile>>user_id;
                infile.ignore(1);
                infile>>rating;
                infile.ignore(256, '\n');
                //printf("movie:%i user:%i rating:%i\n", movie, user_id, rating);
                if(user_id != 0){
                    count++;
                    if( !Add_Piece(user_id) ){
                        //users[users_entered] = user_id;
                        users_entered++;
                    }
                    //else printf("Couldn't add %i\n", id);
                }
            }
            fprintf(outfile, "%i: %i\n", movie, count);
            //if(i < 10) printf("movie:%i length:%i\n", i, count); 
            infile.close();
        }else printf("Unable to open file %s\n", filename); 
    }
    printf("Number of unique users entered:%i\n", users_entered);
    fclose(outfile);

    users = new Int_t[users_entered];
    
    count=0;
    Fill_Array(head, &count);
    printf("Filled arrray\n");
  
    printf("\tLeaving Make_Userlist.\n");
    return users_entered;
}    

Int_t
Load_Userlist(){
    printf("\tEntering Load_Userlist\n");
    Int_t id, users_entered=0;
    users = new Int_t[NUSERS];
    ifstream infile (user_filename);
    if (infile.is_open()){
        while ( infile>>id ){
            //if(users_entered%100000==0) printf("Load List users entered:%i\n", users_entered);
            users[users_entered] = id;
            users_entered++;
        }
        infile.close();
    }else printf("Unable to open file %s\n", user_filename); 
    if(NUSERS != users_entered) printf("NUSERS:%i NOT EQUAL users_entered:%i\n", NUSERS, users_entered);
    //printf("Number of unique users entered:%i\n", users_entered);
    printf("\tLeaving Load_Userlist\n");
    return users_entered;
}

void
Write_Userlist(){
    FILE* outfile = fopen (user_filename, "w");
    for(Int_t i=0; i<NUSERS; i++){
        fprintf(outfile, "%i\n", users[i]);
    }
    fclose (outfile);
}

void
Save_List(FILE* outfile, Piece* temp){
    if(temp != NULL){
        Save_List(outfile, temp->less);
        fprintf(outfile, "%i\n", temp->id);
        Save_List(outfile, temp->more);
    }
}

Int_t
Size_List(Piece* temp){
    if(temp != NULL)
        return 1 + Size_List(temp->less) + Size_List(temp->more);
    return 0;
}

void
Fill_Array(Piece* temp, Int_t* count){
    if(temp != NULL){
        Fill_Array(temp->less, count);
        users[*count]=temp->id;
        (*count)++;
        Fill_Array(temp->more, count);
    }

}

void
Write_Predictions(Bool_t real, Bool_t mine){
    printf("\tEntering Write_Predictions\n");
//Watch out for 0's being read in
    if(real || mine){
        printf("Make Real:%i\tMake Mine:%i\n", real, mine);
        Int_t col_count = 0;
        Int_t* col_pos = new Int_t[NMOVIES]; 
        Int_t npredict=0;
        ifstream infile (probe_filename);
        if (infile.is_open()){
            Int_t movie, user_id, user_index;
            Int_t real_rating;
            Float_t guess;
            Char_t* line;
            Bool_t new_movie;
            FILE* predictfile=NULL;
            FILE* ratingfile=NULL;
            if(mine) predictfile = fopen (predict_filename, "w");
            if(real) ratingfile  = fopen (rating_filename,  "w");
            
            while(!infile.eof() ){
                new_movie = false;
                line = Form(EMPTYLINE);
                infile.getline(line, 25);
                for(Int_t i=0; i<25; i++){//if it contains a colon
                    if(line[i] == ':'){
                        new_movie = true;
                        break;
                    }
                }
                if(new_movie){
                    movie = atoi(line);
                    if(mine) fprintf(predictfile, "%i:\n", movie);
                    if(real) fprintf(ratingfile,  "%i:\n", movie);
                }else{
                    user_id = atoi(line);//get customer name
                    if(user_id != 0){
                        if(mine){
                            user_index = Index_Lookup(user_id);
                            //printf("user_id:%i user_index:%i\n", user_id, user_index);
                            guess = Make_Prediction(movie, user_index);
                            fprintf (predictfile, "%.2f\n", guess);
                            //printf("Made Prediction movie:%i user_id:%i user_index:%i pred:%.2lf\n",movie,user_id, user_index, guess);
                        }
                        if(real){
                            real_rating = Real_Rating(movie, user_id);
                            fprintf (ratingfile,  "%i\n", real_rating);
                        }
                        
                        //if(npredict%100000 == 0) printf("Made %i predictions\n", npredict);
                        npredict++;
                    }else{
                        //printf("Found user 0 in movie %i: %s\n", movie, line);
                    }
                }
            }
            
            if(mine) fclose (predictfile);
            if(real) fclose (ratingfile);
            infile.close();
        }else printf("Unable to open file %s\n", probe_filename);
        printf("Made %i Predictions\n", npredict);
        delete[] col_pos;
    }else printf("Didn't Want Real Ratings or My Predictions\n");
    printf("\tLeaving Write_Predictions\n");
}

Float_t
Make_Prediction(Int_t movie_id, Int_t user_id){
    Int_t movie_index = NUSERS + (movie_id - 1);
    Float_t dist = sqrt(Distance2(movie_index, user_id));
    Float_t prediction;
    if       (STEP_MODE == 0){
        prediction = 5.0*MAXDIST/(4.0*dist+MAXDIST); 
    }else if (STEP_MODE == 1){
        prediction = 4.0*(1-dist/MAXDIST)+1.0;
    }else{
        prediction = -4.0*pow(dist/MAXDIST,2) + 5.0;
    }    
    //printf("ANS = %.2lf\n", prediction);
    //Float_t prediction = 5.0*RSCALE;
    return RSCALE*prediction;
}

Float_t
Compare_Predictions(){
//Calculate RMSE for probe
    printf("\tEntering Compare_Predictions\n");
    Float_t rmse = -1;
    ifstream predictfile (predict_filename);
    if (predictfile.is_open()){
        ifstream ratingfile (rating_filename);
        if (ratingfile.is_open()){
            Int_t movie, count=0;
            Float_t sum2=0;
            Float_t real_rating, my_rating;
            Char_t* predictline;
            Char_t* ratingline;
            Bool_t new_movie;
       
            TH1F hreal("hreal", "Real Ratings", 100, 0, 6);
            TH1F hmine("hmine", "Mine Ratings", 100, 0, 6);
            
            TH2F hreg("hreg", "Real vs Mine", 100, 0, 6, 100, 0, 6);
            
            while(!ratingfile.eof() && !predictfile.eof()){
                new_movie = false;
                ratingline = Form(EMPTYLINE);
                predictline = Form(EMPTYLINE);

                ratingfile.getline(ratingline, 25);
                predictfile.getline(predictline,25);
                for(Int_t i=0; i<10; i++){//if it contains a colon
                    if(ratingline[i] == ':'){
                        new_movie = true;
                        break;
                    }
                }
                if(new_movie){  
                    movie = atoi(ratingline);
                    //printf("New Movie: %s\n", probeline);
                }else{
                    real_rating = atoi(ratingline);
                    my_rating = atof(predictline);
                    if(real_rating && my_rating){
                        if(count%10000 == 0) printf("Count:%i real:%.2lf my:%.2lf\n", count, real_rating, my_rating);
                        hreal.Fill(real_rating);
                        hmine.Fill(my_rating);
                        hreg.Fill(my_rating, real_rating); //(x,y)
                        count++;
                        sum2 += pow(real_rating-my_rating,2); 
                    }
                }
            }
            //printf("sum2:%.2lf count:%i\n",sum2, count);
            rmse = pow(sum2/count, 0.5);

            //TPaveStats *st;
            TCanvas c_recommend("c_recommend", "Recommend Plot");
            hreal.SetLineColor(kGreen);
            hmine.SetLineColor(kRed);
            hreal.Draw();
          
            //st = (TPaveStats*)hmine.FindObject("stats");
            //st->SetY1NDC(0.7); //new x start position
            //st->SetY2NDC(0.8); //new x end position
            hmine.Draw("sames");
            c_recommend.SaveAs(recdist_filename);

           
            TCanvas c_reg("c_reg", "Regression Plot");
            hreg.Draw();
            c_reg.SaveAs(reg_filename);
            
            FILE* outfile;
            outfile = fopen (rmse_filename, "w");
            fprintf(outfile, "RMSE: %.2lf\n", rmse);
            fclose (outfile);
            
            ratingfile.close();
        }else printf("Unable to open file %s\n", probe_filename);
        predictfile.close();
    }else printf("Unable to open file %s\n", predict_filename);
    printf("\tLeaving Compare_Predictions\n");
    return rmse;
}

Int_t
Real_Rating(Int_t movie, Int_t user_id){
    Char_t line[25];
    Int_t id, rating=0;
    Char_t* in_filename = Filename(movie);//what should be the call?
    ifstream infile (in_filename);
    if (infile.is_open()){
        infile.getline(line, 10); //get movie name
        while( !infile.eof() ){
            Parse(&infile, &id, &rating);//, &year, &month, &day);
            if(user_id == id) break;
        }
        
        infile.close();
    }else printf("Unable to open file %s\n", in_filename);
    
    return rating;
}

void
Optimize_RSCALE(){
    //RSCALE = ( SUM(RT*RM) ) / ( SUM(RM*RM) )
    printf("\tEntering Optimize_RSCALE\n");
    printf("RSCALE is %.2lf before\n", RSCALE);
    Int_t real;
    Float_t mine;
    Float_t sum_rt_rm = 0;
    Float_t sum_rm_rm = 0;
    //open predict file
    //open rating file
    ifstream predictfile (predict_filename);
    if (predictfile.is_open()){
        ifstream ratingfile (rating_filename);
        if (ratingfile.is_open()){
            Int_t count=0;
            Int_t movie;
            Char_t* predictline;
            Char_t* ratingline;
            Bool_t new_movie;
            
            while(!ratingfile.eof() && !predictfile.eof()){
                new_movie = false;
                ratingline = Form(EMPTYLINE);
                predictline = Form(EMPTYLINE);

                ratingfile.getline(ratingline, 25);
                predictfile.getline(predictline,25);
                for(Int_t i=0; i<10; i++){//if it contains a colon
                    if(ratingline[i] == ':'){
                        new_movie = true;
                        break;
                    }
                }
                if(new_movie){  
                    movie = atoi(ratingline);
                    //printf("New Movie: %s\n", probeline);
                }else{
                    count++;
                    real = atoi(ratingline);
                    mine = atof(predictline);
                    if( real && mine){
                        sum_rt_rm += real*mine;
                        sum_rm_rm += mine*mine;
                    }
                }
            }
            //printf("count:%i sum_rt_rm:%.2lf sum_rm_rm:%.2lf\n", count, sum_rt_rm, sum_rm_rm);
            RSCALE = sum_rt_rm / sum_rm_rm;
            printf("RSCALE is now %.2lf.\n", RSCALE);
            
            ratingfile.close();
        }else printf("Unable to open file %s\n", rating_filename);
        predictfile.close();
    }else printf("Unable to open file %s\n", predict_filename);
    printf("\tLeaving Optimize_RSCALE\n");
}

void Test(){
    /*time_t start_iterate, end_iterate;
    time_t start_step, end_step;
    Double_t dif_iterate, dif_step;
    time(&start_iterate);
    Iterate();
    time(&end_iterate);
    dif_iterate = difftime(end_iterate, start_iterate);
    printf("Iterate Loop took %.2lf seconds or %.2lf hours to run.\n", dif_iterate, dif_iterate/3600);

    time(&start_step);
    Step(NUSERS, 1488844, 3);
    time(&end_step);
    dif_step = difftime(end_step, start_step);
    printf("Iterate Loop took %.2lf seconds or %.2lf hours to run.\n", dif_step, dif_step/3600);
    */
    Int_t count = 0;
    Int_t* col_pos = new Int_t[NMOVIES]; 
    Int_t* mov_pos = new Int_t[NMOVIES];
    //Int_t npredict=0;
    ifstream infile (probe_filename);
    if (infile.is_open()){
        //Char_t c;
        count=0;
        Int_t limit;
        while( !infile.eof() && count < 17000){
            infile.ignore(999999, ':');
            //infile.ignore(1);
            //printf("AA\n");
            Int_t pos = 1 + infile.tellg();
            col_pos[count] = pos;

            if(count > 0){
                infile.unget();
                while( infile.peek() != '\n'){
                    infile.unget();
                }
                infile.ignore(1);
                Int_t pos2 = infile.tellg();
                mov_pos[count] = pos2;
                infile.seekg( col_pos[count] );
            }else mov_pos[count] = 0;
            //c = infile.peek();
            //printf("pos is %i next char is %c\n", pos, c);
            
            count++;            
            limit = pos;
        }
        
        //printf("limit is %i\n", limit);
        infile.clear();
        infile.seekg (0, ios::beg);
        //Int_t pos = infile.tellg();
        //printf("pos is %i\n", pos );

        count=0;
        Int_t input, movie;
        while( count < 1700){
            //printf("count:%i mov_pos:%i col_pos:%i\n", count, mov_pos[count], col_pos[count]);
            infile.seekg( mov_pos[count] );
            infile>>movie;
            printf("movie %i ", movie);
            infile.seekg( col_pos[count] );
            //Int_t pos = infile.tellg();
            //printf("pos %i col_pos:%i\n", pos, col_pos[col_count]);
            infile>>input;
            printf("input %i\n", input);
            infile.ignore(256, '\n');
            count++;
        }
    }else printf("file not open\n");
    delete[] col_pos;
    delete[] mov_pos;
    return;
}


Int_t
Recommend_Tables(){
    printf("\tEntering Recommend_Tables\n");

    Int_t user_id, movie_index, user_index, rating, count, input_movie;
    Char_t* in_filename;
    Int_t* num_predict=NULL;
    ifstream recfile(numpredict_filename);
    if(recfile.is_open()){
        num_predict = new Int_t[NMOVIES];
   
        for(Int_t movie=1; movie<NMOVIES+1; movie++){
            recfile.ignore(256, ':');
            recfile>>num_predict[movie-1];
            recfile.ignore(256, '\n');
        }
        recfile.close();
    }else printf("Unable to open file %s\n", numpredict_filename);

    for(Int_t movie=1; movie<NMOVIES+1;movie++){
        //open file 'movie'
        movie_index = NUSERS + (movie-1);
        in_filename = Filename(movie);
        ifstream infile (in_filename);
        if (infile.is_open()){
            count = 0;
            if(movie%2000==0) printf("movie: %i\n", movie);
            //printf("movie: %i length is %i\n", movie_index,num_predict[movie_index]);
            recommend[movie-1].length     = num_predict[movie-1];
            recommend[movie-1].user_id    = new Int_t[recommend[movie-1].length];
            recommend[movie-1].rating     = new Int_t[recommend[movie-1].length];
            recommend[movie-1].user_index = new Int_t[recommend[movie-1].length];
            
            infile>>input_movie; //get movie name
            infile.ignore(256, '\n');
            if(input_movie != movie) printf("movie:%i input_val:%i\n", movie, input_movie);
            while( !infile.eof()){
                //Parse(&infile, &user_id, &rating);//, &year, &month, &day);
                infile>>user_id;
                infile.ignore(1);
                infile>>rating;
                infile.ignore(256, '\n');
                //printf("movie:%i user:%i rating:%i\n", movie, user_id, rating);
                if(user_id != 0){
                    //if(count%1000==0) printf("file:%s count:%i\n", in_filename, count);
                    user_index = Index_Lookup(user_id);
                    recommend[movie-1].user_id   [count] = user_id;
                    recommend[movie-1].rating    [count] = rating;
                    recommend[movie-1].user_index[count] = user_index;
                    count++;
                }
            }
            
            infile.close();
        }else printf("Unable to open file %s\n", in_filename);
    }
    delete[] num_predict;
    printf("\tLeaving Recommend_Tables\n");
    return 0;
}
