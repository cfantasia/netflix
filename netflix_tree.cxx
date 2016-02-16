/* Author:Cory Fantsia
 * Netflix Project
 *
 * To Do:
 * Work on Step Fn ->done
 * Movement size depends on rating -> 3 modes of movement
 * Run for a bit -> 10 iterations
 * Make a class -> Universe, Particle, Info
 * create mega arrays rather than disk reads -> done
 * dynamic sizing of myu -> done
 * add points fn ->IP!
 * need user index array
 * npredicts array -> done
 * recommends array -> done
 * make universe class equal to array of particle class objects->done
 * getline vs >> vs ignore
 * should prediction be rounded?
 * create new data member of Info called day for date of recommend
 * make new special dim depend on day
 * Should move direction me random?
 * Add in new movies a few at a time until equalibrium is reached!
 * 5 rating should not give 0 dist!
 */

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom.h"
#include "TPaveStats.h"
#include <fstream>

//#include "consts.h"
#include "Info.h"
#include "Piece.h"
#include "Universe.h"
//#include "Storage_dict.h"
#include "Storage.h"
//#include "Storage.cxx"

using namespace std;

Int_t* users_array = NULL;
Universe* myu=NULL;

//Info* recommend = NULL;
//Int_t num_predict[myu->nmovies];
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
void Fill_Rectree(TTree*, Storage*);
void Fill_Postree(TTree*, Storage*, Bool_t);
void Step(Particle*, Particle*, Short_t);
void Iterate();

void Load_Stats(Int_t*, Int_t*, Int_t*, Float_t*);
void Write_Stats(Int_t, Int_t, Int_t, Float_t);
Int_t Load_Positions();
Int_t Write_Positions();

Int_t Index_Lookup(Int_t*, Int_t, Int_t);
Int_t Index_Lookup(Piece*, Int_t*);
//void Make_Userlist();
Int_t Add_Piece(Piece* head, Int_t);
Int_t Make_Userlist(Int_t);
Int_t Count_Movies();
Int_t Count_Users(Int_t);
void Fill_Userlist(Storage*);

void Write_Userlist();
Int_t Load_Userlist(Int_t);
void Save_List(FILE*, Piece*);

Int_t Size_List(Piece*);
void Fill_Array(Piece*, Int_t*);

void Write_Predictions(Bool_t, Bool_t);
Float_t Compare_Predictions();
Float_t Make_Prediction(Int_t, Int_t);
Short_t Real_Rating(Int_t, Int_t);
void Optimize_RSCALE();
Float_t Calc_U();

void Test();

Int_t Recommend_Tables();
Int_t Find_Nearby(Particle*, Particle*, Int_t, Float_t);
//Int_t Load_numpredict();

Int_t STEP_COUNT;
Float_t MOVE_TTL;

Int_t special=1;//6734;
Int_t numspecial = TMath::Min(200, 17770);

Float_t gstyle_orig = gStyle->GetStatY();

int main(int argc, char ** argv){
    if(argc < 2 || argc > 6 || argv[1][1]=='?'){
		fprintf(stderr,"%s usage: %s nloops [step_mode] [ndim] [predict=0] [reset]\n",argv[0],argv[0]);
        fprintf(stderr, "  step_mode=-1: Loop over 0<STEP_MODE<4\n");
        fprintf(stderr, "  ndim     =-1: Loop over 0<NDIM<15\n");
        exit( 1 );
	}
    
    Int_t input_default[6]={0,1,2,3,0,0};
    Int_t input        [6];
    for(Int_t i=1; i<6; i++){
        if(i>=argc || argv[i][0] =='!') input[i] = input_default[i];
        else                            input[i] = atoi(argv[i]);
    }
    if(atoi(argv[1]) == -1) Test();
    else netflix(input[1], input[2], input[3], input[4], input[5]);
    
    return 0;
}

void
netflix(Int_t nloops, Int_t step_mode, Int_t ndim, Bool_t predict, Int_t reset){
    TROOT myroot("myroot","We're testing TROOT here");

    TFile* frec = new TFile("recommendations.root","UPDATE");
    TTree* trec;
    
    
      
    TFile* fpos= new TFile(Form("simulation_%i_%02i.root", step_mode, ndim),"UPDATE");
    Storage* params = new Storage();
    TTree* tpos;

    trec = (TTree*)frec->Get("recommendations");
    if(!trec) Fill_Rectree(trec, params);
    
    tpos = (TTree*)fpos->Get("universe");
    if(!tpos) Fill_Postree(tpos, params, true);

    //TTree* trec = (TTree*);
    //t3->Write("",TObject::kOverwrite);

    trec->Write("",TObject::kOverwrite);
    tpos->Write("",TObject::kOverwrite);
    frec->Close();
    fpos->Close();
    
    Int_t list_size, myu_size;
    time_t start, end;
    time_t start_iterate, end_iterate;
    time_t start_fn, end_fn;
    Double_t dif, dif_iterate;

    Int_t nusers, nmovies, nentries;
    Float_t step_size;
    time(&start);
    printf("Going to iterate %i time(s) in mode %i in %i dim(s).\n", nloops, step_mode, ndim);
    myu = new Universe();
    ifstream statsfile (myu->stats_filename);
    if(reset==2 || !statsfile){
        printf("Making file %s\n", myu->stats_filename);
        statsfile.close();
        nmovies = Count_Movies();
        nusers = Make_Userlist(nmovies);
        nentries = nmovies + nusers;
        step_size = 0.5;

        myu->Set_Prop(nmovies, nusers, step_mode, ndim);
        myu->Set_Stepsize(step_size);
        Write_Stats(nmovies, nusers, nentries, step_size);
        Write_Userlist();
    }else{
        printf("Loading %s\n", myu->stats_filename);
        statsfile.close();
        Load_Stats(&nmovies, &nusers, &nentries, &step_size);
        myu->Set_Prop(nmovies, nusers, step_mode, ndim);
        myu->Set_Stepsize(step_size);
        ifstream userfile (myu->user_filename);
        if(!userfile){ 
            userfile.close();
            printf("%s does not exist.\n", myu->user_filename);
            list_size = Make_Userlist(nmovies);
            Write_Userlist();
        }else{
            userfile.close();
            printf("%s exists.\n", myu->user_filename);
            list_size = Load_Userlist(nusers);
        }  
        if(list_size != nusers){
            printf("File lengths not equal!!!!\n");
            return;
        }
    }
    
    //users_array = new Int_t[myu->nusers];
    //recommend = new Info[myu->nmovies];
    //myu = new Universe[myu->nentries];
    if(nloops > 0) Recommend_Tables();
        
    ifstream posfile (myu->pos_filename);
    if(reset || !posfile){ 
        printf("%s does not exist.\n", myu->pos_filename);
        posfile.close();
        myu_size = Initialize(); 
    }else{
        printf("%s exists.\n", myu->pos_filename);
        posfile.close();
        myu_size = Load_Positions();
    }
    
    //delete[] users_array;
    // 0  1  9  10 11 20 21 NLOOPS

    // 0  1  1  1  2  2  3  WANTED

    //-1  0  8  9  10 19 20 i-1
    // 0  0  0  0  1  1  2  i-1)/10
    // 1  1  1  1  2  2  3  i-1)/10 + 1

    
    // 0  1  9  10 11 20 21 NLOOPS
    
    // 0  1  1  2  2  3  3  WANTED

    // 0  0  0  1  1  2  2  i/10
    // 1  1  1  2  2  3  3  i/10 + 1

    //tpos->GetUserInfo()->Add();

    Float_t prev_ttl=0;
    Int_t hour_count = 1;
    Float_t U=0;
    TH1F hmove("hmove", "MOVE_TTL vs loop #", nloops-1, 0, nloops-1);
    TH1F huplot("huplot", "U vs loop#", 1+(nloops-1)/10, 0, nloops-1);
    TH1F hssize("hssize", "Step_SIZE vs loop #", nloops-1, 0, nloops-1);
    for(Int_t loop=0; loop<nloops; loop++){
        printf("\tBeginning Loop %i of %i\n", loop+1, nloops);
        time(&start_iterate);
        STEP_COUNT = 0;
        MOVE_TTL = 0;
        Iterate();
        hmove.Fill(loop, MOVE_TTL);
        hssize.Fill(loop, myu->STEP_SIZE);
        
        printf("movie %i pos is ", special);
        for(Int_t dim=0; dim<myu->ndim; dim++){
            printf(" %.4lf ",myu->movies[special-1].x[dim]);
        }printf("\n");
        
        if(loop%10 == 0){
            U = Calc_U();
            printf("U is %.2f\n", U);
            huplot.Fill(loop, U);
        }
        printf("STEP_COUNT = %i\n", STEP_COUNT);
        printf("MOVE_TTL = %.6lf STEP_SIZE = %.6lf\n", MOVE_TTL, myu->STEP_SIZE);
        if(loop>0 && MOVE_TTL<prev_ttl && MOVE_TTL>0.999999*prev_ttl) myu->Set_Stepsize(myu->STEP_SIZE/2);
        prev_ttl = MOVE_TTL;
        time(&end_iterate);
        dif_iterate = difftime(end_iterate, start_iterate);
        printf("Iterate Loop %i of %i took %.2lf seconds or %.2lf hours to run.\n", loop+1, nloops, dif_iterate, dif_iterate/3600.0);

        if(myu->STEP_SIZE < 0.000000001){
            printf("Move size=%.9f is below threshold.  Breaking\n", myu->STEP_SIZE); 
            break;
        }
        if(difftime(end_iterate, start) > hour_count*3600){
            printf("Writing Position Data after %i hour(s)\n", hour_count);
            hour_count++;
            Write_Positions();
            Write_Stats(myu->nmovies, myu->nusers, myu->nentries, myu->STEP_SIZE);
            Fill_Postree(tpos, params, false);
            tpos->Write("",TObject::kOverwrite);
        }
    }
    if(nloops > 0){
        TCanvas cmove("cmove", "");
        cmove.Divide(1,3);
        cmove.cd(1);
        //gPad->SetLogy(1);
        hmove.Draw();
        hmove.Write("",TObject::kOverwrite);

        cmove.cd(2);
        huplot.Draw("L");
        huplot.Write("",TObject::kOverwrite);
  
        cmove.cd(3);
        //gPad->SetLogy(1);
        hssize.Draw("L");
        hssize.Write("",TObject::kOverwrite);

        cmove.SaveAs(myu->movettl_filename);
    }
    /*
    //printf("Should be 1468 %i\n", (myu->movies + 1469-1)->index);
    Int_t num;
    for(Int_t i=0; i<myu->nmovies; i++){
        num = Find_Nearby(myu->movies + i, myu->movies, myu->nmovies, 0.1);
        if(num > 100) printf("movie %i count:%i\n", i+1, num);
        }
    */
    if(myu->ndim >= 2) Display();
    
    Write_Positions();
    Write_Stats(myu->nmovies, myu->nusers, myu->nentries, myu->STEP_SIZE);
    Fill_Postree(tpos, params, false);
    tpos->Write("",TObject::kOverwrite);
    
    if(predict){
        time(&start_fn);
        Bool_t make_real=false, make_mine=true;//defaults
        if(reset == 2){
            make_real=true;
            make_mine=true;
        }else{
            ifstream realfile (myu->rating_filename);
            if( !realfile) make_real=true;
            realfile.close();

            ifstream minefile (myu->predict_filename);
            if( !minefile) make_mine=true;
            minefile.close();
        }

        Write_Predictions(make_real,make_mine); 
        time(&end_fn);
        printf("Write_Predictions took %li seconds to run.\n", end_fn-start_fn);

        printf("RMSE is %.2lf BEFORE\n", Compare_Predictions());
        
        Optimize_RSCALE();
        
        printf("Going to redo Write_Predictions...\n");
        Write_Predictions(false, true); //Rewrite predictions using new RSCALE 
                
        printf("RMSE is %.2lf AFTER\n", Compare_Predictions());
    }
    //TTree* etree = new TTree("etree", "evolve tree");
    //Fill_Tree(myu, etree);
    //etree->BuildIndex("index");
    //Display(etree);

    Write_Stats(myu->nmovies, myu->nusers, myu->nentries, myu->STEP_SIZE);
    
    if(users_array) delete[] users_array;
    if( myu)  delete myu;

    fpos->Write("",TObject::kOverwrite);
    fpos->Close();   
    //frec->Write("",TObject::kOverwrite);
    frec->Close();
    
    time(&end);
    dif = difftime(end, start);
    printf("Program took %.2lf seconds or %.2lf hours to run.\n", dif, dif/3600.0);
}


Float_t
Distance2(Particle* part1, Particle* part2, Float_t* dist_x2=NULL){
    if(part1 != part2){
        Float_t sum2=0, d2;
        Float_t x1, x2, dx;
             
        for(Int_t dim=0; dim<myu->ndim; dim++){
            x1 = part1->x[dim];
            x2 = part2->x[dim];
            dx = TMath::Abs(x1-x2);   
            if( dx < (myu->MAXX-myu->MINX)/2) d2 = pow(dx,2);
            else d2 = pow(myu->MAXX-myu->MINX-dx,2);
            
            if(dist_x2 != NULL) dist_x2[dim] = d2;
            sum2 += d2;
        }
        return  sum2;
    }
    printf("Distance between same point %i is 0!\n", part1->index);
    return 0;
}


Float_t
Distance(Particle* part1, Particle* part2, Float_t* delta_x2=NULL){
    return sqrt( Distance2(part1, part2, delta_x2));
}

Float_t
Calc_U(){
    Info* temp;
    Int_t len;
    Float_t U=0, dist2;
    Int_t movie_index;
    for(Int_t movie=1; movie<myu->nmovies+1; movie++){
        temp = myu->recommend + (movie-1);
        len = temp->length;
        movie_index = movie-1;
        for(Int_t j=0; j<len; j++){
            //rating is < RTHRESH equiv to same sign particles
            //Cory: What should I do here?
            dist2 = Distance2(myu->movies + movie_index, myu->users + temp->user_index[j]);//Cory: check 2nd
            if       (myu->step_mode == 0){ //Cory: Formula
                //prediction = 5.0*myu->MAXDIST/(4.0*dist+myu->MAXDIST); 
                U += -1.0*(temp->rating[j]-myu->RTHRESH)/sqrt(dist2);
            }else if (myu->step_mode == 1){
                //prediction = 4.0*(1-dist/myu->MAXDIST)+1.0;
                U +=  1.0*(temp->rating[j]-myu->RTHRESH)*sqrt(dist2);
            }else{
                //prediction = -4.0*pow(dist/myu->MAXDIST,2) + 5.0;
                U +=  1.0*(temp->rating[j]-myu->RTHRESH)*dist2;
            }
        }
    }    
    return U;
}

void
Calc_rating(){
    /*
    Float_t d;
    Float_t ratodist, oneodist;
 
    for(Int_t i=0; i<myu->nentries; i++){
        ratodist = 0;
        oneodist = 0;
        for(Int_t j=0; j<myu->nentries; j++){
            if(i != j){
                //tree->GetEntryWithIndex(i);
                d = Distance(i, j);
                ratodist += (myu->part[j].rating)/d;
                oneodist += 1.0/d;
            }
        }
        myu->part[i].rating =  ratodist/oneodist;  
    }
    */
}

/*
Float_t
Potential(Particle* myu, Int_t index){
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
    myu->nentries++; //update this
    */

    //users_array, myu, maybe recommend myu->nusers, myu->nentries
    if( is_user ){
        Int_t* users_array2 = new Int_t[myu->nusers+1];
        Int_t offset = 0;
        for(Int_t i=0; i<myu->nusers+1; i++){
            if(users_array[i-offset] != id){
                if(i != myu->nusers+1){
                    if(offset || users_array[i-offset] < id){
                        users_array2[i] = users_array[i-offset];
                    }else{
                        users_array2[i] = id;
                        offset = 1;
                    }
                }else{
                    if(offset) users_array2[i] = users_array[i-offset];
                    else       users_array2[i] = id;
                }
            }else{
                printf("Tried to add an existing member %i\n!", id);
                delete[] users_array2;
                return;
            }
        }
        delete[] users_array;
        users_array = users_array2;
        myu->nusers++;
    }else{  // It's a movie!
        //Add entry to recommend
        Info* recommend2 = new Info[myu->nmovies+1];
        for(Int_t i=0; i<myu->nmovies; i++){
            
        }
        delete[] myu->recommend;
        myu->recommend = recommend2;
        myu->nmovies++;
    }  
    myu->nentries = myu->nusers + myu->nmovies;
    Universe* myu2 = new Universe[myu->nentries];
    for(Int_t i=0; i<myu->nentries; i++){
        //insert new entry
    }
    delete[] myu;
    myu = myu2;
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
    Info* temp=NULL;
    Int_t len, movie_index;
    for(Int_t movie=/*1*/special; movie<special+numspecial/*myu->nmovies+1*/; movie++){//Cory: Make upper limit special
        //printf("Iterating with movie:%i\n", movie);
        temp = myu->recommend + movie-1;
        len = temp->length;
        movie_index = movie - 1;
        for(Int_t i=0; i<len; i++){
            //if(i>len-2) printf("movie:%i id:%i index[%i] = %i\n", movie, temp->user_id[i], i, temp->user_index[i]);
            Step(myu->movies + movie_index, myu->users + temp->user_index[i], temp->rating[i]);
            STEP_COUNT++;
        }
    }
}
/*
void
Iterate2(){
    Info* temp=NULL;
    Int_t len, movie_index, rating;
    trec->SetBranchAddress("movie", &movie);
    trec->SetBranch????
    for(Int_t movie=special; movie<special+numspecial; movie++){//Cory: Make upper limit special
        //printf("Iterating with movie:%i\n", movie);
        temp = myu->recommend + movie-1;
        len = hnrec->GetBinContent(i);//temp->length;
        movie_index = movie - 1;
        for(Int_t i=0; i<len; i++){
            trec->GetEntry(count);
            //if(i>len-2) printf("movie:%i id:%i index[%i] = %i\n", movie, temp->user_id[i], i, temp->user_index[i]);
            Step(myu->movies + movie_index, myu->users + temp->user_index[i], temp->rating[i]);
            STEP_COUNT++;
        }
    }
}
*/
Int_t
Sign(Particle* part1, Particle* part2, Short_t rating, Int_t dim){
//returns the sign of the move
/*
   1 if away && index1.x > index2.x  //what if index1.x[0] > but index1.x[1] <
   1 if towards && index2.x > index2.x
  -1 else
*/
    Int_t dir = -1;//NOTICE THIS!
    Float_t dx = part1->x[dim] - part2->x[dim];  
          
    if( dx > 0 ) dir *=  1;
    else         dir *= -1;

    if( TMath::Abs(dx) > (myu->MAXX-myu->MINX)/2) dir*=-1;
    else                                          dir*= 1;
    /*
    if(rating < myu->RTHRESH) dir *=  1;
    else                 dir *= -1;   
    */
    return dir;
}

Float_t
Rating_Mag(Short_t rating){
    //Should this fn do anything with rating = RTHRESH?
    return TMath::Abs((rating-myu->RTHRESH)/2.0);
}

Float_t
OutofBounds(Float_t pos){
    while(pos>myu->MAXX || pos <= myu->MINX){
        //printf("pos before:%.4lf ", pos);
        if      (pos >  myu->MAXX){
            pos = pos - myu->MAXX + myu->MINX;
        }else{
            pos = pos - myu->MINX + myu->MAXX;
        }
        //printf("pos after:%.2lf\n", pos);
    }
    
    //if(pos > myu->MAXX || pos <=myu->MINX) printf("Still out of bounds pos=%.2lf\n", pos);
    return pos;
}

void
Step(Particle* part1, Particle* part2, Short_t rating){
/*should two objects move closer or further based on rating
 * move = STEP_SIZE*(sign*mag)*(dist_x2/dist2)
 * should this be devided by 1/dist2 or 1/distx2
 * 
 * step_mode: 0=1/r
 *            1=random
 *            2=Hooke
 */
    Float_t sign, move, mag_move, rat_mag;
    Float_t* dist_x2 = new Float_t[myu->ndim];
    Float_t dist2 = Distance2(part1, part2, dist_x2);
    Float_t scale=0; //4*(MAX-1)/4 + 1 = MAX
    Float_t d0 = (myu->MAXDIST-2.0)*(5.0-rating)/4.0 + 1.0; //Formula
                                                //Cory:should this be max_dist??? or depend on it squared?
    if      (myu->step_mode == 0) scale = (myu->STEP_SIZE/dist2)/(sqrt(dist2)-d0); //Cory: don't know about this!
    else if (myu->step_mode == 1) scale = (myu->STEP_SIZE/dist2); 
    else                          scale = (myu->STEP_SIZE/dist2)* (sqrt(dist2) - d0);
    //rat_mag = Rating_Mag(rating);

    for(Int_t dim=0; dim<myu->ndim; dim++){
        sign    = Sign(part1, part2, rating, dim);
        mag_move = scale*dist_x2[dim]; //Cory: took out rat_mag
      
        //mag_move = TMath::Min(mag_move, myu->MAXDIST/2);
        //mag_move = TMath::Max(mag_move, (Float_t)sqrt(dist_x2[dim])/2);

        move = sign*mag_move;
        part1->x[dim] += 0.5*move;  //Cory:!!!!!!!!!!!!!!!!!!!!!!!
        part2->x[dim] -= 0.5*move;
        
        part1->x[dim] = OutofBounds( part1->x[dim] );
        part2->x[dim] = OutofBounds( part2->x[dim] );

        if(STEP_COUNT %20000000 == 0) printf("d2=%.4lf dx2=%.4lf |move| is %.4lf\n", dist2, dist_x2[dim], TMath::Abs(move)); 
        MOVE_TTL += TMath::Abs(move);
    }
    delete[] dist_x2;
}

void
Display(){
    TCanvas* canvas1;
    canvas1 = (TCanvas*)gROOT->FindObject("canvas1");
    
    if(canvas1) delete canvas1;
    canvas1 = new TCanvas("canvas1" "The Universe");

    Float_t* pos;
    if(myu->ndim == 2){
        TH2F* husers=NULL;
        TH2F* hmovies=NULL;
        TH2F* hspecial=NULL;
        
        husers  = (TH2F*)gROOT->FindObject("husers");
        hmovies = (TH2F*)gROOT->FindObject("hmovies");
    
        if(husers)  delete husers;
        if(hmovies) delete hmovies;
      
        husers  = new TH2F("husers", "The Universe", 40, myu->MINX, myu->MAXX, 40, myu->MINX, myu->MAXX);
        hmovies = new TH2F("hmovies", "The Universe", 40, myu->MINX, myu->MAXX, 40, myu->MINX, myu->MAXX);
        hspecial = new TH2F("hspecial", "The Universe", 40, myu->MINX, myu->MAXX, 40, myu->MINX, myu->MAXX);
    
        hspecial->SetMarkerColor(kBlue);
        husers->SetMarkerColor(kGreen);
        hmovies->SetMarkerColor(kRed);
        printf("!!!!!!!!!!!!!\n");
        //printf("6 is %i\n", Index_Lookup(myu->recommend[special-1].user_id, myu->recommend[special-1].length-1, 6)); 
        //for(Int_t i=0; i< 10; i++) printf("%i\n", Index_Lookup(myu->recommend[special-1].user_id, myu->recommend[special-1].length-1, myu->users[i].id));
        /*
        for(Int_t i=0; i<10 ; i++){
            Int_t marker=-1;
            Bool_t kg = true;
            for(Int_t j=0; kg && j<myu->recommend[special-1].length-1; j++){
                if(myu->users[i].id == myu->recommend[special-1].user_id[j]){
                    marker = j;
                    kg = false;
                }
            }
            if(marker != -1) printf("i:%i id:%i index:%i\n", i, myu->recommend[special-1].user_id[marker], marker);
            else             printf("Couldn't find %i\n", myu->users[i].id );
        }
        */
        Float_t sum2=0;
        Int_t countrmse = 0;
        for(Int_t mov=special; mov<special+numspecial; mov++){
            hmovies->Fill(myu->movies[mov-1].x[0], myu->movies[mov-1].x[1]);
            for(Int_t i=0; i<myu->nusers/*myu->nentries*/; i++){ //Cory: Make upper limit nusers + 1
                if(mov == special){
                    Int_t marker=-1;
                    Bool_t kg = true;
                    for(Int_t j=0; kg && j<myu->recommend[mov-1].length; j++){
                        if(myu->users[i].id == myu->recommend[mov-1].user_id[j]){
                            marker = j;
                            kg = false;
                        }
                    }
                    if( marker != -1 ){
                        Float_t dist = Distance(myu->movies+mov-1, myu->users+i);
                        Short_t real = myu->recommend[mov-1].rating[marker];
                        Float_t mine = 5.0-4.0*(dist-1.0)/(myu->MAXDIST-1);//Formula
                        sum2 += pow(real-mine, 2);
                        countrmse++;
                        if(i > myu->nusers-10) printf("found i:%i id:%i marker:%i dist:%.2f rat:%i pred:%.2f\n", i, myu->users[i].id, marker, dist, real, mine);
                        if(i < myu->nusers) pos=myu->users[i].x;
                        else                pos=myu->movies[i-myu->nusers].x;

                        hspecial->Fill(pos[0], pos[1]);
                    }
                }
                pos = myu->users[i].x;
                if(pos[0] > myu->MAXX || pos[0] <= myu->MINX) printf("Out of Bounds Index:%i pos[0]:%.lf\n", i, pos[0]);
                if(pos[1] > myu->MAXX || pos[1] <= myu->MINX) printf("Out of Bounds Index:%i pos[1]:%.lf\n", i, pos[1]);
                
                if(i<myu->nusers){
                    //husers->Fill(pos[0], pos[1]);
                }else{
                    //hmovies->Fill(pos[0], pos[1]);
                }
            }
        }
        printf("RMSE for me is %.6f\n", sqrt(sum2/countrmse));

        gStyle->SetStatY(gstyle_orig);
        hmovies->Draw();
        hmovies->Write("",TObject::kOverwrite);
        canvas1->Update();

        gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
        husers->Draw("sames");
        husers->Write("",TObject::kOverwrite);
        canvas1->Update();
        
        gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
        hspecial->Draw("sames");
        hspecial->Write("",TObject::kOverwrite);
        
        canvas1->SaveAs(myu->universe_filename);
        
        delete husers;
        delete hmovies;
        delete hspecial;
    }else if(myu->ndim >= 3){
        TH3F* husers=NULL;
        TH3F* hmovies=NULL;
        TH3F* hspecial=NULL;
        
        husers  = (TH3F*)gROOT->FindObject("husers");
        hmovies = (TH3F*)gROOT->FindObject("hmovies");
        hspecial = (TH3F*)gROOT->FindObject("hspecial");
    
        if(husers)  delete husers;
        if(hmovies) delete hmovies;
        if(hspecial) delete hspecial;

        husers  = new TH3F("husers", "The Universe", 40, myu->MINX, myu->MAXX, 40, myu->MINX, myu->MAXX,40, myu->MINX, myu->MAXX);
        hmovies = new TH3F("hmovies", "The Universe", 40, myu->MINX, myu->MAXX, 40, myu->MINX, myu->MAXX,40, myu->MINX, myu->MAXX);
        hspecial = new TH3F("hspecial", "The Universe", 40, myu->MINX, myu->MAXX, 40, myu->MINX, myu->MAXX,40, myu->MINX, myu->MAXX);
                
        hspecial->SetMarkerColor(kBlue);
        husers->SetMarkerColor(kGreen);
        hmovies->SetMarkerColor(kRed);
        printf("!!!!!!!!!!!!!\n");
     
        Float_t sum2=0;
        Int_t countrmse = 0;
        for(Int_t mov=special; mov<special+numspecial; mov++){
            hmovies->Fill(myu->movies[mov-1].x[0], myu->movies[mov-1].x[1], myu->movies[mov-1].x[2]);
            for(Int_t i=0; i<myu->nusers/*myu->nentries*/; i++){ //Cory: Make upper limit nusers + 1
                if(mov == special){
                    Int_t marker=-1;
                    Bool_t kg = true;
                    for(Int_t j=0; kg && j<myu->recommend[mov-1].length; j++){
                        if(myu->users[i].id == myu->recommend[mov-1].user_id[j]){
                            marker = j;
                            kg = false;
                        }
                    }
                    if( marker != -1 ){
                        Float_t dist = Distance(myu->movies+mov-1, myu->users+i);
                        Short_t real = myu->recommend[mov-1].rating[marker];
                        Float_t mine = 5.0-4.0*(dist-1.0)/(myu->MAXDIST-2.0);//Formula
                        sum2 += pow(real-mine, 2);
                        countrmse++;
                        if(i > myu->nusers-10) printf("found i:%i id:%i marker:%i dist:%.2f rat:%i pred:%.2f\n", i, myu->users[i].id, marker, dist, real, mine);
                        if(i < myu->nusers) pos=myu->users[i].x;
                        else                pos=myu->movies[i-myu->nusers].x;

                        hspecial->Fill(pos[0], pos[1], pos[2]);
                    }
                }
                pos = myu->users[i].x;
                if(pos[0] > myu->MAXX || pos[0] <= myu->MINX) printf("Out of Bounds Index:%i pos[0]:%.lf\n", i, pos[0]);
                if(pos[1] > myu->MAXX || pos[1] <= myu->MINX) printf("Out of Bounds Index:%i pos[1]:%.lf\n", i, pos[1]);
                if(pos[2] > myu->MAXX || pos[2] <= myu->MINX) printf("Out of Bounds Index:%i pos[2]:%.lf\n", i, pos[2]);
                
                if(i<myu->nusers){
                    //husers->Fill(pos[0], pos[1]);
                }else{
                    //hmovies->Fill(pos[0], pos[1]);
                }
            }
        }
        printf("RMSE for me is %.6f\n", sqrt(sum2/countrmse));

        gStyle->SetStatY(gstyle_orig);
        hmovies->Draw();
        hmovies->Write("",TObject::kOverwrite);
        canvas1->Update();

        gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
        husers->Draw("sames");
        husers->Write("",TObject::kOverwrite);
        canvas1->Update();
        
        gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
        hspecial->Draw("sames");
        hspecial->Write("",TObject::kOverwrite);
        canvas1->SaveAs(myu->universe_filename);
      
        delete husers;
        delete hmovies;
        delete hspecial;
    }

    delete canvas1;
}

            
void
Display(TTree* tree){
    Char_t* field = Form("");;
    for(Int_t i=0; i<myu->ndim; i++){
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
    Float_t min[myu->ndim];
    Float_t max[myu->ndim];
    for(Int_t i=0; i<myu->ndim; i++){
        TAxis* axis;
        if(i == 0) axis = htemp->GetXaxis();
        if(i == 1) axis = htemp->GetYaxis();
        if(i == 2) axis = htemp->GetZaxis();

        min[myu->ndim] = axis->GetXmin();
        max[myu->ndim] = axis->GetXmax();
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
    
    //tree->ReadFile("cory.dat", "evt/i:t/i:c0/f:c1/f:c2/f:c3/f:c4/f:c5/f:c6/f:c7/f:c8/f:c9/f");
    Particle* temp;
    gRandom->SetSeed(0); //Argument = 0 means use computer clock
    printf("Current Seed is %i\n", gRandom->GetSeed());
  
    for(Int_t i=0; i<myu->nentries; i++){
        if(i<myu->nusers){
            temp = myu->users + i;
            temp->index = i;
            temp->id = users_array[i];
            temp->is_user=true;
        }else{
            temp = myu->movies + i-myu->nusers;
            temp->index = i-myu->nusers;
            temp->id = i+1 - myu->nusers;
            temp->is_user=false;
        }
        
        for(Int_t j=0; j<myu->ndim; j++){
            temp->x[j] = gRandom->Uniform(myu->MINX, myu->MAXX);
            //myu->part[i].x[j] = gRandom->Gaus(50,10);//myu->part[i].prop[j];
        }
        
        //gRandom->Rndm() < 0.5) charge = -1;
        //Add_point(myu, index, prop, x);
        
    }
    //Calc_U(myu);
    //Calc_rating(myu);

    printf("\tLeaving Initialize\n");
    return myu->nentries; //return # of users assigned pos
}

void
Fill_Postree(TTree* ptree, Storage* params, Bool_t init=false){
    printf("\tEntering Fill_Postree\n");
    Int_t ndim = params->fndim;  
    Particle* ptemp;// = new Particle(ndim);
    printf("D\n");
        
    printf("K\n");
    //TBranch* bparticle = ptree->Branch("particle", ptemp, Form("index/I:id/I:x[%i]/F:is_user/O", ndim));
    
    TBranch* bindex = ptree->Branch("index",  &ptemp->index, "index/I");
    TBranch* bid = ptree->Branch("id",     &ptemp->id, "id/I");
    TBranch* bx = ptree->Branch("x", ptemp->x, Form("x[%i]/F", ndim)); 
    TBranch* bis_user = ptree->Branch("is_user", &ptemp->is_user, "is_user/O");
    
    if(init){
        ptemp = new Particle(ndim);
        //bparticle->SetAddress(ptemp);
        
        bindex->SetAddress(&ptemp->index);
        bid->SetAddress(&ptemp->id);
        bx->SetAddress(ptemp->x);
        bis_user->SetAddress(&ptemp->is_user);
        
        gRandom->SetSeed(0); //Argument = 0 means use computer clock
        printf("Current Seed is %i\n", gRandom->GetSeed());
        TH1I* hu = params->huserlist;
        printf("M\n");
        for(Int_t i=0; i<params->fnmovies; i++){
            ptemp->index = i;
            ptemp->id = i+1;
            ptemp->is_user=false;
            
            for(Int_t j=0; j<ndim; j++){
                ptemp->x[j] = gRandom->Uniform(params->fminx, params->fmaxx);
                //xtemp[j] = ptemp->x[j];
            }
            ptree->Fill();
            //if(i<2) ptemp->Print(ndim);
        }
      
        printf("N\n");
        for(Int_t i=0; i<params->fnusers; i++){
            //printf("i:%i\n",i);
            ptemp->index = i;
            ptemp->id = (Int_t) hu->GetBinContent(i+1);
            ptemp->is_user=true;
            
            for(Int_t j=0; j<ndim; j++){
                ptemp->x[j] = gRandom->Uniform(params->fminx, params->fmaxx);
                //xtemp[j] = ptemp->x[j];
            }
            ptree->Fill();
        }
        printf("T\n");
        delete ptemp;
        
    }else{
        //delete ptemp;
        //branch->SetAddress(Void *address)
           
        for(Int_t i=0; i<params->fnmovies; i++){
            ptemp = myu->movies + i;
            //bparticle->SetAddress(ptemp);
            
            bindex->SetAddress(&ptemp->index);
            bid->SetAddress(&ptemp->id);
            bx->SetAddress(ptemp->x);
            bis_user->SetAddress(&ptemp->is_user);
            
            if(i<2) ptemp->Print(ndim);
            ptree->Fill();
        }
        printf("# entries is %li\n", (long Int_t)ptree->GetEntries() );
        
        for(Int_t i=0; i<params->fnusers; i++){
            ptemp = myu->users + i;
            //bparticle->SetAddress(ptemp);
            
            bindex->SetAddress(&ptemp->index);
            bid->SetAddress(&ptemp->id);
            bx->SetAddress(ptemp->x);
            bis_user->SetAddress(&ptemp->is_user);
            
            ptree->Fill();
        }
        printf("# entries is %li\n", (long Int_t)ptree->GetEntries() );
    }

    ptree->Show(0);
    ptree->Show(1);
    
    printf("\tExiting Fill_Postree\n");
}

void
Fill_Rectree(TTree* rtree, Storage* params){
    printf("\tEntering Fill_Rectree\n");

    Int_t movie;
    Int_t year;//[5];
    Int_t month;//[5];
    Int_t day;//[5];
    Int_t user_id;//[5];
    //Int_t user_index;//[5];
    Int_t rating;//[5];
           
    rtree->Branch("movie",  &movie, "movie/I");
    rtree->Branch("user_id", &user_id, "user_id/I");
    rtree->Branch("rating", &rating, "rating/I");
    rtree->Branch("year", &year, "year/I");
    rtree->Branch("month", &month, "month/I");
    rtree->Branch("day", &day, "day/I");

    printf("Y\n");
    //printf("params->fnmovies: %i\n", params->fnmovies);
    params->hnrec = new TH1I("hnrec", "", params->fnmovies, 1,params->fnmovies +1);
    printf("Z\n");
    //rtree->GetUserInfo()->Add(hnrec);
    
    Char_t* filename;
    Int_t count;
    for(Int_t i=1; i<params->fnmovies; i++){ 
        if(i%5000==1) printf("i:%i\n",i);
        filename = Form("/home/fantasia/work/netflix/data/training_set/mv_%07i.txt",i);
        FILE* infile = fopen (filename, "r");
        if( !infile ){
            printf("Couldn't open file %s\n", filename);
            break;
        }
        
        fscanf(infile, "%i: ", &movie);
        //printf("movie: %i\n", movie);
        count = 0;
        
        while( ! feof (infile) ){
        
            if (fscanf(infile, "%d,%d,%d-%d-%d ", &user_id, &rating, &year, &month, &day) ) count++;
            if (count<3 && i==1) printf("%i %i %i %i %i\n", user_id, rating, year, month, day);
            //if( fgets(line, 256, infile) ) count++;

            rtree->Fill();
        }
        
        params->hnrec->Fill(movie, count);
        if (i<3) printf("count is %i\n", count);
                      
        fclose(infile);
    }
    
    rtree->Show(0);
    rtree->Show(1);
    printf("\tLeaving Fill_Rectree\n");
}

Int_t
Count_Movies(){
    Char_t* filename;
    Int_t count = 0;
    while(1){
        filename = Form("/home/fantasia/work/netflix/data/training_set/mv_%07i.txt", count+1);
        //if(count < 10) printf("%s\n", filename);
        //Filename(count+1); //movies start at 1
        ifstream infile (filename);
        if(!infile) break;
        count++;
    }
    printf("%i unique movies found\n", count);
    return count;
}

void
Fill_Hist(Piece* temp, Int_t* count, TH1I* husers){
    if(temp != NULL){
        Fill_Hist(temp->less, count, husers);
        husers->Fill(*count,temp->id);
        (*count)++;
        Fill_Hist(temp->more, count, husers);
    }
}

Int_t
Count_Users(Int_t nmovies){
    Char_t* filename;
    Int_t users_entered = 0;
    Int_t movie;
    Int_t user_id;
    Piece* head=NULL;
    
    for(Int_t i=1; i<nmovies+1; i++){
        filename = Form("/home/fantasia/work/netflix/data/training_set/mv_%07i.txt", i);
        if(i%5000==1) printf("file:%s users_entered: %i\n", filename, users_entered);
        ifstream infile (filename);
        if (infile.is_open()){
            infile>>movie;//get movie name
            infile.ignore(256, '\n');
            //printf("Read in movie:%i\n", movie);
            if(i != movie) printf("i:%i movie:%i\n", i, movie);  
            while (infile>>user_id){
                infile.ignore(256, '\n');
                //printf("movie:%i user:%i\n", movie, user_id);
                //printf("head before:%p\n", head);
                if(!head){
                    head = new Piece(user_id);
                    users_entered++;
                }else
                    if( !Add_Piece(head, user_id) ) users_entered++;
                //else printf("Couldn't add %i\n", id);
                //printf("head after:%p\n", head);
            }
            infile.close();
        }else printf("Unable to open file %s\n", filename); 
    }
    printf("Number of unique users entered:%i\n", users_entered);
    
    if(head) delete head;
    return users_entered;
}

void
Fill_Userlist(Storage* params){
    Char_t* filename;
    Int_t users_entered = 0;
    Int_t movie;
    Int_t user_id;
    Piece* head=NULL;
    
    for(Int_t i=1; i<(params->fnmovies+1); i++){
        filename = Form("/home/fantasia/work/netflix/data/training_set/mv_%07i.txt", i);
        if(i%5000==1) printf("file:%s users_entered: %i\n", filename, users_entered);
        ifstream infile (filename);
        if (infile.is_open()){
            infile>>movie;//get movie name
            infile.ignore(256, '\n');
            //printf("Read in movie:%i\n", movie);
            if(i != movie) printf("i:%i movie:%i\n", i, movie);  
            while (infile>>user_id){
                infile.ignore(256, '\n');
                //printf("movie:%i user:%i\n", movie, user_id);
                //printf("head before:%p\n", head);
                if(!head){
                    head = new Piece(user_id);
                    users_entered++;
                }else{
                    if( !Add_Piece(head, user_id) ) users_entered++;
                    //else printf("Couldn't add %i\n", id);
                    //printf("head after:%p\n", head);
                }
            }
            infile.close();
        }else printf("Unable to open file %s\n", filename); 
    }
    printf("Number of unique users entered:%i\n", users_entered);
    params->huserlist = new TH1I("huserlist", "user ids", users_entered, 0, users_entered);

    TCanvas ctest("ctest", "ctest");

    Int_t count = 0;
    Fill_Hist(head, &count, params->huserlist);
    
    params->huserlist->Draw();
    //params->Add_Userlist(huserlist);
    //printf("inside address: %p\n", params->huserlist);
    //ctest.SaveAs("AAA.ps");
    
    if(head) delete head;
}


void
Write_Stats(Int_t nmovies, Int_t nusers, Int_t nentries, Float_t step_size){
    FILE* outfile;
    outfile = fopen (myu->stats_filename, "w");
    fprintf(outfile, "nentries: %i\n", nentries);
    fprintf(outfile, "nmovies: %i\n", nmovies);
    fprintf(outfile, "nusers: %i\n", nusers);
    fprintf(outfile, "step_size: %f\n", step_size);
    fclose (outfile);
}

void
Load_Stats(Int_t* nmovies, Int_t* nusers, Int_t* nentries, Float_t* step_size){
    ifstream statsfile(myu->stats_filename);
    if(statsfile.is_open()){
        statsfile.ignore(256, ':');
        statsfile>>*nentries;
        statsfile.ignore(256, '\n');

        statsfile.ignore(256, ':');
        statsfile>>*nmovies;
        statsfile.ignore(256, '\n');

        statsfile.ignore(256, ':');
        statsfile>>*nusers;
        statsfile.ignore(256, '\n');

        statsfile.ignore(256, ':');
        statsfile>>*step_size;

        printf("nentries:%i nmovies:%i nusers:%i\n", *nentries, *nmovies, *nusers);
    }else printf("Countn't open file %s\n", myu->stats_filename);
}

Int_t
Write_Positions(){
   //Format:
    //index, id, x[0], ..., x[myu->ndim-1], is_user 
    Particle* temp;
    
    FILE* outfile;
    //printf("About to open %s for writing with base:%s\n", myu->pos_filename, myu->posbase_filename);
    outfile = fopen (myu->pos_filename, "w");
    for (Int_t i=0; i<myu->nentries; i++){
        //printf("i:%i\n", i);
        if(i<myu->nusers) temp = myu->users  + i;
        else              temp = myu->movies + i - myu->nusers;

        fprintf (outfile, "%i, %i, ",temp->index, temp->id);
        for(Int_t dim=0; dim<myu->ndim; dim++){
            fprintf (outfile, "%f, ",temp->x[dim]);
        }
        fprintf (outfile, "%i\n",temp->is_user);
    }
    fclose (outfile);
    return myu->nentries;
}

Int_t
Load_Positions(){
    printf("\tEntering Load_Positions.\n");
    Char_t c;
    Int_t obj_entered=0;
    ifstream infile (myu->pos_filename);
    if (infile.is_open()){
        Particle* temp = myu->users;
        while ( infile>>temp->index >> c >> temp->id >> c){
            //index, id, x[0], ..., x[ndim-1], is_user 
            //if(obj_entered%100000 == 0) printf("obj_entered:%i\n", obj_entered);
            

            for(Int_t dim=0; dim<myu->ndim; dim++) infile>>temp->x[dim]>>c;
            infile>>temp->is_user;
            infile.ignore(256, '\n');
            obj_entered++;
            if(obj_entered < myu->nusers) temp = myu->users + obj_entered;
            else                          temp = myu->movies + obj_entered-myu->nusers;
            /*
            infile>>temp->index;//.getline(line, 20, ','); //get index 
            infile.ignore(256, ',');//temp->index = atoi(line);
            
            infile>>temp->id;//infile.getline(line, 20, ','); //get ID
            infile.ignore(256, ',');//temp->id = atoi(line);
            
            for(Int_t dim=0; dim<myu->ndim; dim++){
                infile>>temp->x[dim];//infile.getline(line, 20, ','); //get x[dim]
                infile.ignore(256, ',');//temp->x[dim] = atof(line); //needs to be a float
            }
            infile>>temp->is_user;
            */
           
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
    else printf("Unable to open file %s\n", myu->pos_filename);
    printf("Total # obj_entered:%i\n", obj_entered);
    printf("\tLeaving Load Positions.\n");
    return obj_entered;
}

Int_t
Index_Lookup(Int_t* array, Int_t last, Int_t value){
    //Given the user_id, find the index number;
    Int_t first=0;
        
    while (first <= last) {
        Int_t mid = (first + last) / 2;  // compute mid point.
        if (value > array[mid]) 
            first = mid + 1;  // repeat search in top half.
        else if (value < array[mid]) 
            last = mid - 1; // repeat search in bottom half.
        else
            return mid;     // found it. return position /////
    }
    //printf("ALERT!!!! Could not find user_id=%i\n", user_id);
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
Add_Piece(Piece* head, Int_t user_id){
    Piece* temp=head;
    //printf("temp: %p\n", temp);
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
                return -1; //already in list
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
        printf("head:%p\n", head);
    }
    return 0;
}

Int_t
Make_Userlist(Int_t nmovies){
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
    Int_t user_id, movie;
    Short_t rating;

Piece* head=NULL;

    FILE* outfile;
    outfile = fopen (myu->numpredict_filename, "w");
    
    for(Int_t i=1; i<nmovies+1; i++){
        count = 0;
        filename = Filename(i);
        if(i%5000==1) printf("file:%s users_entered: %i\n", filename, users_entered);
        ifstream infile (filename);
        if (infile.is_open()){
            infile>>movie;//get movie name
            infile.ignore(256, '\n');
            //printf("Read in movie:%i\n", movie);
            if(i != movie) printf("i:%i movie:%i\n", i, movie);  
            while (infile>>user_id){
                infile.ignore(1);
                infile>>rating;
                infile.ignore(256, '\n');
                //printf("movie:%i user:%i rating:%i\n", movie, user_id, rating);
                count++;
                if( !Add_Piece(head, user_id) ) users_entered++;
                //else printf("Couldn't add %i\n", id);
            }
            fprintf(outfile, "%i: %i\n", movie, count);
            if(i >= special && i < special+numspecial) printf("movie:%i length:%i\n", i, count); 
            infile.close();
        }else printf("Unable to open file %s\n", filename); 
    }
    printf("Number of unique users entered:%i\n", users_entered);
    fclose(outfile);

    users_array = new Int_t[users_entered];
    
    count=0;
    Fill_Array(head, &count);
    printf("Filled arrray with %i members\n", count);
    for(Int_t i=0; i<10; i++) printf("%i: %i\n", i, users_array[i]);  
    if(head) delete head;
    printf("\tLeaving Make_Userlist.\n");
    return users_entered;
}    

Int_t
Load_Userlist(Int_t nusers){
    printf("\tEntering Load_Userlist\n");
    Int_t id, users_entered=0;
    users_array = new Int_t[nusers];
    ifstream infile (myu->user_filename);
    if (infile.is_open()){
        while ( infile>>id ){
            infile.ignore(256, '\n');
            //if(users_entered%100000==0) printf("Load List users entered:%i\n", users_entered);
            users_array[users_entered] = id;
            users_entered++;
        }
        infile.close();
    }else printf("Unable to open file %s\n", myu->user_filename); 
    if(nusers != users_entered) printf("nusers:%i NOT EQUAL users_entered:%i\n", nusers, users_entered);
    //printf("Number of unique users entered:%i\n", users_entered);
    printf("\tLeaving Load_Userlist\n");
    return users_entered;
}

void
Write_Userlist(){
    FILE* outfile = fopen (myu->user_filename, "w");
    for(Int_t i=0; i<myu->nusers; i++){
        fprintf(outfile, "%i\n", users_array[i]);
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
        users_array[*count]=temp->id;
        (*count)++;
        Fill_Array(temp->more, count);
    }
}

void
Write_Predictions(Bool_t real, Bool_t mine){
    printf("\tEntering Write_Predictions\n");
    if(real || mine){
        printf("Make Real:%i\tMake Mine:%i\n", real, mine);
        //Int_t col_count = 0;
        //Int_t* col_pos = new Int_t[myu->nmovies]; 
        Int_t npredict=0;
        ifstream infile (myu->probe_filename);
        if (infile.is_open()){
            Int_t movie=-1, user_id, user_index;
            Short_t real_rating;
            Int_t value;
            Float_t guess;
            //Char_t* line;
            //Bool_t new_movie;
            FILE* predictfile=NULL;
            FILE* ratingfile=NULL;
            if(mine) predictfile = fopen (myu->predict_filename, "w");
            if(real) ratingfile  = fopen (myu->rating_filename,  "w");
            
            while( !infile.eof() ){
                while( infile>>value && infile.peek() != ':'){
                    user_id = value;
                    if(mine){
                        user_index = Index_Lookup(users_array, myu->nusers-1, user_id);
                        //printf("user_id:%i user_index:%i\n", user_id, user_index);
                        guess = Make_Prediction(movie, user_index);
                        fprintf (predictfile, "%.2f\n", guess);
                        //printf("Made Prediction movie:%i user_id:%i user_index:%i pred:%.2lf\n",movie,user_id, user_index, guess);
                    }
                    if(real){
                        real_rating = Real_Rating(movie, user_id);
                        fprintf (ratingfile,  "%i\n", real_rating);
                        if(npredict%100000 == 0) printf("Made %i predictions\n", npredict);
                    }
                                        
                    npredict++;
                }
                if( infile.peek() == ':' ){
                    movie = value;
                    if(mine) fprintf(predictfile, "%i:\n", movie);
                    if(real) fprintf(ratingfile,  "%i:\n", movie);
                    infile.ignore(256, '\n');
                }
            }
            
            if(mine) fclose (predictfile);
            if(real) fclose (ratingfile);
            infile.close();
        }else printf("Unable to open file %s\n", myu->probe_filename);
        printf("Made %i Predictions\n", npredict);
        //delete[] col_pos;
    }else printf("Didn't Want Real Ratings or My Predictions\n");
    printf("\tLeaving Write_Predictions\n");
}

Float_t
Make_Prediction(Int_t movie, Int_t user_index){
    Int_t movie_index = movie-1;
    Float_t dist = sqrt(Distance2(myu->movies + movie_index, myu->users + user_index));
    Float_t prediction;
    if       (myu->step_mode == 0){
        prediction = 5.0*myu->MAXDIST/(4.0*dist+myu->MAXDIST); //Formula 
    }else if (myu->step_mode == 1){
        prediction = 4.0*(1-dist/myu->MAXDIST)+1.0;
    }else{//HERE
        prediction = 5.0 - 4.0*(dist-1.0)/(myu->MAXDIST-2.0);//-5.0*pow(dist/myu->MAXDIST,2) + 5.5;
    }    
    //printf("ANS = %.2lf\n", prediction);
    //Float_t prediction = 5.0*myu->RSCALE;
    //if(myu->RSCALE == 1.0) return round(myu->RSCALE*prediction);
    Float_t m=0.1;
    return myu->RSCALE*(m*prediction + (1-m)*round(prediction));
}

Float_t
Compare_Predictions(){
//Calculate RMSE for probe
    printf("\tEntering Compare_Predictions\n");
    Float_t rmse = -1;
    ifstream predictfile (myu->predict_filename);
    if (predictfile.is_open()){
        ifstream ratingfile (myu->rating_filename);
        if (ratingfile.is_open()){
            Int_t count=0;
            Float_t sum2=0;
            Short_t real_rating;
            Float_t my_rating;
            //Char_t* predictline;
            //Char_t* ratingline;
            //Bool_t new_movie;
       
            TH1F hreal("hreal", "Real Ratings", 100, 0, 6);
            TH1F hmine("hmine", "Mine Ratings", 100, 0, 6);
            TH2F hreg("hreg", "Real vs Mine", 100, 0, 6, 100, 0, 6);
            while(!ratingfile.eof() && !predictfile.eof()){
                while ( ratingfile>>real_rating  && ratingfile.peek() != ':'){
                    predictfile>>my_rating;// && predictfile.peek() != ':'){
                    hreal.Fill(real_rating);
                    hmine.Fill(my_rating);
                    hreg.Fill(my_rating, real_rating); //(x,y)
                    count++;
                    sum2 += pow(real_rating-my_rating,2); 
                    //if(count>1400000 == 0) printf(" real:%i mine:%f\n", real_rating, my_rating); 
                }
                predictfile>>my_rating;
                if (ratingfile.peek() == ':' && predictfile.peek() == ':'){
                    //a movie line
                    //if(count > 1400000 == 0) printf("MOVIES real:%i mine:%f\n", real_rating, my_rating); 
                    ratingfile.ignore(256, '\n');
                    predictfile.ignore(256, '\n');
                }
            }
            
            //printf("sum2:%.2lf count:%i\n",sum2, count);
            rmse = pow(sum2/count, 0.5);

            //TPaveStats *st;
            TCanvas c_recommend("c_recommend", "Recommend Plot");
            hreal.SetLineColor(kGreen);
            hmine.SetLineColor(kRed);

            gStyle->SetStatY(gstyle_orig);
            hreal.Draw();
            hreal.Write("",TObject::kOverwrite);
            c_recommend.Update();
            gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
            hmine.Draw("sames");
            hmine.Write("",TObject::kOverwrite);
            c_recommend.SaveAs(myu->recdist_filename);

            gStyle->SetStatY(gstyle_orig);
            TCanvas c_reg("c_reg", "Regression Plot");
            hreg.Draw("lego2");
            hreg.Write("",TObject::kOverwrite);
            c_reg.SaveAs(myu->reg_filename);
            
            FILE* outfile;
            outfile = fopen (myu->rmse_filename, "w");
            fprintf(outfile, "RMSE: %.2lf\n", rmse);
            fclose (outfile);
            
            ratingfile.close();
        }else printf("Unable to open file %s\n", myu->probe_filename);
        predictfile.close();
    }else printf("Unable to open file %s\n", myu->predict_filename);
    printf("\tLeaving Compare_Predictions\n");
    return rmse;
}

Short_t
Real_Rating(Int_t movie, Int_t user_id){
    Char_t line[25];
    Int_t id;
    Short_t rating=0;
    Char_t* in_filename = Filename(movie);//what should be the call?
    ifstream infile (in_filename);
    if (infile.is_open()){
        infile.getline(line, 10); //get movie name
        while( !infile.eof() ){
            while(infile>>id){
                if(user_id == id) break;
                infile.ignore(256, '\n');
            }
        }
        if( infile.eof() ) printf("Couldn't find id:%i in movie:%i\n", user_id, movie);
        infile.close();
    }else printf("Unable to open file %s\n", in_filename);
    
    return rating;
}

void
Optimize_RSCALE(){
    //RSCALE = ( SUM(RT*RM) ) / ( SUM(RM*RM) )
    printf("\tEntering Optimize_RSCALE\n");
    printf("RSCALE is %.2lf before\n", myu->RSCALE);
    Int_t real;
    Float_t mine;
    Float_t sum_rt_rm = 0;
    Float_t sum_rm_rm = 0;
    //open predict file
    //open rating file
    ifstream predictfile (myu->predict_filename);
    if (predictfile.is_open()){
        ifstream ratingfile (myu->rating_filename);
        if (ratingfile.is_open()){
            Int_t count=0;
            //Int_t movie;
            //Char_t* predictline;
            //Char_t* ratingline;
            //Bool_t new_movie;
            
            while(!ratingfile.eof() && !predictfile.eof()){
                while( ratingfile>>real && ratingfile.peek() != ':'){
                    //Its not a movie!
                    predictfile>>mine;
                    sum_rt_rm += real*mine;
                    sum_rm_rm += mine*mine;
                    count++;
                }
                predictfile>>mine;
                if(ratingfile.peek() == ':' && predictfile.peek() == ':'){
                    //its a movie!
                    //printf("real:%i mine:%f\n", real, mine); 
                    ratingfile.ignore(256, '\n');
                    predictfile.ignore(256, '\n');
                }
                /*    
                new_movie = false;
                ratingline = Form(myu->EMPTYLINE);
                predictline = Form(myu->EMPTYLINE);

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
                    */
            }
            //printf("count:%i sum_rt_rm:%.2lf sum_rm_rm:%.2lf\n", count, sum_rt_rm, sum_rm_rm);
            myu->RSCALE = sum_rt_rm / sum_rm_rm;
            printf("RSCALE is %.2lf after\n", myu->RSCALE);
            
            ratingfile.close();
        }else printf("Unable to open file %s\n", myu->rating_filename);
        predictfile.close();
    }else printf("Unable to open file %s\n", myu->predict_filename);
    printf("\tLeaving Optimize_RSCALE\n");
}


Int_t
Recommend_Tables(){
    printf("\tEntering Recommend_Tables\n");

    Int_t user_id, user_index, count, input_movie;
    Short_t rating;
    Char_t* in_filename;
    Char_t c;
    Int_t* num_predict=NULL;
    ifstream nrecfile(myu->numpredict_filename);
    if(nrecfile.is_open()){
        //printf("myu->nmovies is %i\n", myu->nmovies);
        num_predict = new Int_t[myu->nmovies];
   
        for(Int_t movie=1; movie<myu->nmovies+1; movie++){
            nrecfile>>input_movie>>c;//ignore(256, ':');
            nrecfile>>num_predict[input_movie-1];
            nrecfile.ignore(256, '\n');
        }
        nrecfile.close();
    }else printf("Unable to open file %s\n", myu->numpredict_filename);

    for(Int_t movie=special/*1*/; movie<special+numspecial/*myu->nmovies+1*/;movie++){ //Cory: Changed this too
        //open file 'movie'
        //movie_index = (movie-1);
        in_filename = Filename(movie);
        ifstream infile (in_filename);
        if (infile.is_open()){
            count = 0;
            if(movie%5000==0) printf("movie: %i\n", movie);
            if(movie >= special && movie < special+numspecial) printf("movie: %i length is %i\n", movie, num_predict[movie-1]);
           
            infile>>input_movie; //get movie name
            infile.ignore(256, '\n');
            if(input_movie != movie) printf("movie:%i input_val:%i\n", movie, input_movie);

            myu->recommend[movie-1].Prep_Info(num_predict[movie-1]);
            /*
            //Parse(&infile, &user_id, &rating);//, &year, &month, &day);
            infile>>user_id;
            infile.ignore(1);
            infile>>rating;
            infile.ignore(256, '\n');
            */
            while(infile >> user_id >> c >> rating){
                user_index = Index_Lookup(users_array, myu->nusers-1, user_id);
                myu->recommend[movie-1].user_id   [count] = user_id;
                myu->recommend[movie-1].rating    [count] = rating;
                myu->recommend[movie-1].user_index[count] = user_index;
                count++;
                infile.ignore(256, '\n');
            }
                    
            infile.close();
        }else printf("Unable to open file %s\n", in_filename);
    }
    delete[] num_predict;
    printf("\tLeaving Recommend_Tables\n");
    return 0;
}

void
Find_Allrecmean(){
     TH1F real_dist("real_dist", "All Real Ratings Dist", 100, 0, 6);
    Int_t count;
    Char_t* filename;
    Char_t c;
    Int_t nmovies = 17770;
    Int_t id;
    Short_t rating;
    for(Int_t i=1; i<nmovies+1; i++){
        count = 0;
        filename = Filename(i);
        if(i%1000==1) printf("file:%s count: %i\n", filename, count);
        
        ifstream infile (filename);
        if (infile.is_open()){
            infile.ignore(256, '\n'); //get movie name
            while( !infile.eof() ){    
                while(infile>>id>>c>>rating){
                    real_dist.Fill(rating);
                    infile.ignore(256, '\n');
                    count++;
                }
            }
        }else printf("Counldn't open\n");
    }

    TCanvas c_real_dist("c_real_dist", "All Real Ratings Dist");
    real_dist.Draw();
    c_real_dist.SaveAs(myu->realdist_filename);
}

Int_t
Find_Nearby(Particle* target, Particle* array, Int_t array_len, Float_t dist){
    Int_t count=0;
    Particle* temp;
    Float_t d;    
    for(Int_t i=0; i<array_len; i++){
        //printf("i: %i\n", i);
        temp = array + i;
        if(temp != target){
            d = Distance(target, temp);
            if( d < dist ){
                //printf("Index %i is %.2f away from target index %i\n", temp->index, d, target->index);
                count++;
            }
        }
    }
    //printf("Leaving Find_Nearby\n");
    return count;
}

void Test(){
    /*
    Int_t count = 0;
    Int_t value1, value2, value3;
    ifstream infile ("./output/test.txt");
    if (infile.is_open()){
        Char_t c;
        count=0;
        //while( !infile.eof() ){
            infile>>value1;
            infile.ignore(256, '\n');
            infile>>value2;
            printf("movie:%i value2:%i\n", value1, value2);
      
            while( infile>>value1 && infile.peek() != ':' && infile>>c>>value2>>c>>value3){
                //infile>>c>>value2>>c>>value3;
                printf("value1: %i value2:%i value3:%i\n", value1, value2, value3);
                count++;            
            }
            if( infile.peek() == ':'){
                printf("movie: %i\n", value1);
                infile.ignore(256, '\n');
            }
    
            //}
    }else printf("file not open\n");
    */

    /*
    time_t start1, start2, start3, start4, end1, end2, end3, end4;
    Double_t dif;
    Char_t c;
    Int_t ivalue;
    Float_t fvalue;
    //Char_t* pos_filename = ("./output/positions/positions_2_10.txt");
    Int_t obj_entered;
    Char_t line[25];
    Char_t* filename =("/home/fantasia/work/netflix/output/mode_2_10/positions.txt");
    Char_t* filename2 =("/home/fantasia/work/netflix/output/mode_2_10/test.txt");
    //printf("Entered OK\n");
    time(&start1);
    obj_entered = 0;
    for(Int_t i=1; i<10+1; i++){
        //printf("Loop %i\n", i);
        ifstream infile1 (filename);
        if (infile1.is_open()){
            while ( !infile1.eof()){
                if( infile1.eof() ) break;
                
                infile1.getline(line, 20, ','); //get index 
                ivalue = atoi(line);

                //if(i==1 && obj_entered < 10) printf("obj_entered:%i index:%i\n", obj_entered, ivalue);
                
                infile1.getline(line, 20, ','); //get ID
                ivalue = atoi(line);
                
                for(Int_t dim=0; dim<10; dim++){
                    infile1.getline(line, 20, ','); //get x[dim]
                    fvalue = atof(line); //needs to be a float
                }
                infile1.getline(line, 20, '\n'); //get ID
                ivalue = atoi(line);
                
                obj_entered++;
            }
            infile1.close();
        }else printf("Unable to open file %s\n", filename);
    }
    time(&end1);
    dif = difftime(end1, start1);
    printf("Program 1 took %.2lf seconds or %.2lf hours to run.\n", dif, dif/3600.0);


    time(&start2);
    obj_entered = 0;
    for(Int_t i=1; i<10 + 1; i++){
        ifstream infile (filename);
        if (infile.is_open()){
            while ( infile>>ivalue>>c>>ivalue>>c){
                //infile.ignore(1);
                //infile>>ivalue;
                //infile.ignore(1);
                //index, id, x[0], ..., x[ndim-1], is_user 
                //if(i==1 && obj_entered < 10) printf("obj_entered:%i index:%i\n", obj_entered, ivalue);
                for(Int_t dim=0; dim<10; dim++){
                    infile>>fvalue>>c;
                    //infile.ignore(1);
                }
                infile>>ivalue;
                infile.ignore(256, '\n');
                obj_entered++;
            }
            infile.close();
        }
        else printf("Unable to open file %s\n", filename);
    }
    time(&end2);
    dif = difftime(end2, start2);
    printf("Program took %.2lf seconds or %.2lf hours to run.\n", dif, dif/3600.0);
        
    time(&start3);
    obj_entered = 0;
    Int_t index, id, is_user;
    Float_t x;
    for(Int_t i=1; i<10 + 1; i++){
        FILE* infile = fopen (filename, "r");
        //fscanf ( FILE * stream, const char * format, ... );
        //fscanf(infile, "%i%*s ", &movie);
        //printf("movie: %i\n", movie);
        while( ! feof (infile) ){
            //if(obj_entered < 10) printf("obj_entered:%i index:%i\n", obj_entered, ivalue);
            fscanf(infile, "%i%*c %i%*c ", &index, &id);
            for(Int_t dim=0; dim<10; dim++){
                fscanf(infile, "%f%*c ", &x);
            }
            fscanf(infile, "%i ", &is_user);
            //if(i==1 && obj_entered < 10) printf("index:%i id:%i x[10]:%f is_user:%i\n", index, id, x, is_user);
            obj_entered++;
        }
        fclose(infile);
    }
    time(&end3);
    dif = difftime(end3, start3);
    printf("Program took %.2lf seconds or %.2lf hours to run.\n", dif, dif/3600.0);

    time(&start4);
    obj_entered = 0;
    Int_t num;
    TTree* tree = new TTree("tree", "");
    for(Int_t i=1; i<10 + 1; i++){
        
        num = tree->ReadFile(filename2, "index/I:id/I:x[10]/F:is_user/O");//"evt/i:t/i:c0/f:c1/f:c2/f:c3/f:c4/f:c5/f:c6/f:c7/f:c8/f:c9/f");
        
        //printf("num lines read in is %i\n", num);
        //tree->Show(0);
        //tree->Show(1);
        //tree->Show(2);
        //tree->Show(3);
        
        //delete tree;
    }
    delete tree;
    
    time(&end4);
    dif = difftime(end4, start4);
    printf("Program took %.2lf seconds or %.2lf hours to run.\n", dif, dif/3600.0);
    */
    /*
    TFile* ftest = new TFile("test.root", "UPDATE");
    TTree* tree = new TTree("tree", "");
    Storage* stest = new Storage(5,11, 2, 2);
    Storage* stest2 = new Storage(100,2, 2, 2);
    stest->Print();
    stest2->Print();
    tree->GetUserInfo()->Add(stest);
    tree->GetUserInfo()->Add(stest2);
    ftest->Write();
    ftest->Close();

    
    TFile* ktest = new TFile("test.root", "NEW");
    if(ktest) printf("THIS IS \n");
    if(ktest->IsOpen()) printf("AA\n");
    printf("BB\n");
    TTree* b = (TTree*) ktest->Get("tree");
    const Char_t* name= b->GetUserInfo()->At(0)->GetName();
    printf("%s\n", name);
    Storage* me = (Storage*) b->GetUserInfo()->FindObject(name);
    me->Print();

    const Char_t* name2= b->GetUserInfo()->At(1)->GetName();
    printf("%s\n", name2);
    Storage* me2 = (Storage*) b->GetUserInfo()->FindObject(name2);
    me2->Print();
    
    ktest->Close();
    
    Universe* temp_u = new Universe();
    printf("%.2f %s\n", temp_u->RTHRESH, temp_u->pos_filename);
    delete temp_u;
    */
    
    Int_t step_mode = 2, ndim = 7;  
    Storage* params=NULL;
    TFile* fpos= new TFile(Form("simulation_%i_%02i.root", step_mode, ndim),"UPDATE");
         
    TTree* tpos=NULL;
    tpos = (TTree*)fpos->Get("tpos");
    printf("A\n");
    if(!tpos){
        tpos = new TTree("tpos", "Position Tree");
        //Count # of movies
        Int_t nmovies = Count_Movies();
        //Count # of users
        Int_t nusers = 480189;//Cory: Count_Users(nmovies);
            
        params = new Storage(nmovies, nusers, step_mode, ndim);
        Fill_Userlist(params);
        //tpos->GetUserInfo()->Add(params);
        params->Print();

        Fill_Postree(tpos, params, true);
    }else{
        //const Char_t* name= tpos->GetUserInfo()->At(0)->GetName();
        //printf("%s\n", name);
        params = (Storage*) fpos->Get("params");
        //params = (Storage*)tpos->GetUserInfo("params");//new Storage();
    }
    //Fill myu position data
    myu = new Universe();
    myu->Set_Prop(params->fnmovies, params->fnusers,
                  params->fstep_mode, params->fndim);
    //Particle* temp = new Particle(ndim);
    /*
    Int_t index, id;
    Float_t* x = new Float_t[ndim];
    Bool_t is_user;
    */
    Particle* ptemp;
    //tpos->SetBranchAddress("particle", ptemp);
    
   
    
    for(Int_t i=0; i<params->fnmovies; i++){
        ptemp = myu->movies + i;

        tpos->SetBranchAddress("index", &ptemp->index);
        tpos->SetBranchAddress("id", &ptemp->id);
        tpos->SetBranchAddress("x", ptemp->x);
        tpos->SetBranchAddress("is_user", &ptemp->is_user);
        
        tpos->GetEntry(i);
        /*       
        ptemp->index = index;
        ptemp->id = id;
        for(Int_t j=0; j<ndim; j++) temp->x[j] = x[j];
        ptemp->is_user = is_user;
        */
        if(i<2) tpos->Show();
        if(i<2) ptemp->Print(ndim);
    }
    for(Int_t i=0; i<params->fnusers; i++){
        ptemp = myu->users + i;

        tpos->SetBranchAddress("index", &ptemp->index);
        tpos->SetBranchAddress("id", &ptemp->id);
        tpos->SetBranchAddress("x", ptemp->x);
        tpos->SetBranchAddress("is_user", &ptemp->is_user);
        
        tpos->GetEntry(i+params->fnmovies);
        /*
        temp->index = index;
        temp->id = id;
        for(Int_t j=0; j<ndim; j++) temp->x[j] = x[j];
        temp->is_user = is_user;
        */
    }
    //delete[] x;
    
    printf("Done with Storage and TPos Stuff\n");
    
    //TFile* frec=new TFile("recommendations.root", "UPDATE");
    TTree* trec=NULL;
    trec = (TTree*)fpos->Get("trec");
    if(!trec){
        trec = new TTree("trec", "Recommendation Table");
        Fill_Rectree(trec, params);
    }
    printf("C\n");
   
    tpos->Write("", TObject::kOverwrite);
    params->Write("params", TObject::kOverwrite);
    trec->Write("", TObject::kOverwrite);

    printf("D\n");

    //Cory: Call Step
    Int_t nloops = 10;
    TH1F* hu = new TH1F("hu", "Potential", nloops, 0, nloops);
    //Iterate2();

    //printf("About to delete tpos\n");
    //delete tpos;
    //printf("Delete tpos\n");

    TTree* tpos2 = new TTree("tpos2", "Position Tree");

    Fill_Postree(tpos2, params, false);
    
    tpos = tpos2;

    tpos->SetObject("tpos", "Position Tree"); 
    //frec->Close();
    //fpos->Close();

    printf("E\n");    

    tpos2->Write("", TObject::kOverwrite);
    params->Write("params", TObject::kOverwrite);
    trec->Write("", TObject::kOverwrite);

    delete myu;
    //delete params;
    return;
}   

