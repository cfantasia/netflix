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
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom.h"
#include "TPaveStats.h"
#include "TThread.h"
#include <fstream>

//#include "consts.h"
#include "Info.h"
#include "Piece.h"
#include "Universe.h"

using namespace std;

Int_t* users_array = NULL;
Universe* myu=NULL;
Piece* head = NULL;
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
 *
 */

struct Worth{
    Particle* part1;
    Particle* part2;
    Short_t rating;
    Float_t* move;
};

void check();

Int_t Initialize();
void Display();
void Display(TTree*);
void netflix(Int_t, Int_t, Int_t, Bool_t, Int_t);
void Fill_Rectree(TTree*);
void Fill_Postree(TTree*);

Float_t Ideal(Short_t);
void Step(Particle*, Particle*, Short_t, Float_t*);
void Iterate();
Float_t OutofBounds(Float_t);
Float_t Distance2(Particle*, Particle*, Float_t*);

void Load_Stats(Int_t*, Int_t*, Int_t*);
void Write_Stats(Int_t, Int_t, Int_t);
Int_t Load_Positions();
Int_t Write_Positions();

Int_t Index_Lookup(Int_t*, Int_t, Int_t);
Int_t Index_Lookup(Piece*, Int_t*);
//void Make_Userlist();
Int_t Add_Piece(Int_t);
Int_t Make_Userlist(Int_t);
Int_t Count_Movies();

void Write_Userlist();
Int_t Load_Userlist(Int_t);
void Save_List(FILE*, Piece*);

Int_t Size_List(Piece*);
void Fill_Array(Piece*, Int_t*);

void Write_Predictions(Bool_t, Bool_t, Int_t);
Float_t Compare_Predictions();
Float_t newCompare_Predictions();
Float_t Make_Prediction(Int_t, Int_t, Int_t);
Short_t Real_Rating(Int_t, Int_t);
void Optimize_RSCALE();
void Calc_Fit_Param();
Bool_t Is_In_Probe(Int_t, Int_t);

Float_t Calc_U();

void Test();

Int_t Recommend_Tables();
Int_t Find_Nearby(Particle*, Particle*, Int_t, Float_t);
//Int_t Load_numpredict();

Int_t STEP_COUNT;
Float_t MOVE_TTL;

Int_t special=1;//6734;
Int_t numspecial = TMath::Min(17770, 17770);

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
    Int_t list_size, myu_size;
    time_t start, end;
    time_t start_iterate, end_iterate;
    time_t start_fn, end_fn;
    Double_t dif, dif_iterate;

    Int_t nusers, nmovies, nentries;
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

        myu->Set_Prop(nmovies, nusers, step_mode, ndim);
        Write_Stats(nmovies, nusers, nentries);
        Write_Userlist();
    }else{
        printf("Loading %s\n", myu->stats_filename);
        statsfile.close();
        Load_Stats(&nmovies, &nusers, &nentries);
        myu->Set_Prop(nmovies, nusers, step_mode, ndim);
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
    if(nloops > 0){
        Int_t tot_nloops;
        Int_t prev_nloops;
        Int_t nupoints;
        Char_t* flag;
        ifstream loopstatsfile (myu->loopstats_filename);
        if( !loopstatsfile){
            printf("%s does not exist\n", myu->loopstats_filename);
            prev_nloops = 0;
            tot_nloops = nloops;
            nupoints = 1+(tot_nloops-1)/10;
            flag = "w";
        }else{
            loopstatsfile.ignore(256, ':');
            loopstatsfile>>prev_nloops;
            tot_nloops = prev_nloops + nloops;
            //fscanf(infile, "NumUpoints: %i ", &nupoints);
            nupoints = 1+(tot_nloops-1)/10;
            flag = "a+";
        }loopstatsfile.close();
        
        printf("prev_nloops: %i\n", prev_nloops);
        printf("tot_nloops: %i\n", tot_nloops);
        printf("nupoints: %i\n", nupoints);
        printf("flag: %s\n", flag);
        
        Float_t prev_ttl=0;
        Int_t hour_count = 1;
        Float_t U=0, ss=0;
        Int_t i, lasti=-1;
        
        TH1F hmove("hmove", "MOVE_TTL vs loop #", tot_nloops-1, 0, tot_nloops-1);
        FILE* mfile = fopen (myu->movefile_filename, flag);//Cory: should this be "a"
        if( flag == "w" ){
            
        }else{
            rewind(mfile);
            while( fscanf(mfile, "%i, %f", &i, &MOVE_TTL) && !feof(mfile) ){
                //printf("i:%i MOVE_TTL:%f eof:%i\n", i, MOVE_TTL, -1); //feof(mfile));
                hmove.Fill(i, MOVE_TTL); 
                lasti = i;
            }
            fflush (mfile);    // flushing or repositioning required
        }
        printf("Loaded %i points in mfile\n", lasti);
        
        TH1F huplot("huplot", "U vs loop#", nupoints, 0, tot_nloops-1);
        FILE* ufile = fopen (myu->ufile_filename, flag);//Cory: should this be "a"
        if( flag == "w"){

        }else{
            rewind(ufile);
            while( fscanf(ufile, "%i, %f", &i, &U) && !feof(ufile)){
                //printf("i:%i U:%f\n", i, U);
                huplot.Fill(i, U);
                lasti = i;
            }
            fflush (ufile);
        }
        printf("Loaded %i points in ufile\n", lasti);
        
        TH1F hssize("hssize", "Step_SIZE vs loop #", tot_nloops-1, 0, tot_nloops-1);
        FILE* ssfile = fopen (myu->ssizefile_filename, flag);//Cory: should this be "a"
        if( flag == "w" ){
            
        }else{
            rewind(ssfile);
            while( fscanf(ssfile, "%i, %f", &i, &ss) && !feof(ssfile)){
                //printf("i:%i ss:%f\n", i, ss);
                hssize.Fill(i, ss);
                lasti = i;
            }
            myu->Set_Stepsize(ss);
            fflush(ssfile);
        }
        printf("Loaded %i point in ssfile\n", lasti);
        
        for(Int_t loop=prev_nloops; loop<tot_nloops; loop++){
            printf("\tBeginning Loop %i of %i\n", loop+1, tot_nloops);
            time(&start_iterate);
            STEP_COUNT = 0;
            MOVE_TTL = 0;
            U = 0;
            Iterate();
            hmove.Fill(loop, MOVE_TTL);
            fprintf(mfile, "%i, %f\n", loop, MOVE_TTL);
            hssize.Fill(loop, myu->STEP_SIZE);
            fprintf(ssfile, "%i, %.9f\n", loop, myu->STEP_SIZE);
            
            printf("movie %i pos is ", special);
            for(Int_t dim=0; dim<myu->ndim; dim++){
                printf(" %.4lf ",myu->movies[special-1].x[dim]);
            }printf("\n");
            
            if(loop%10 == 0){
                U = Calc_U();
                printf("U is %.2f\n", U);
                huplot.Fill(loop, U);
                fprintf(ufile, "%i, %f\n", loop, U);
            }
            printf("STEP_COUNT = %i\n", STEP_COUNT);
            printf("MOVE_TTL = %.6lf STEP_SIZE = %.9lf\n", MOVE_TTL, myu->STEP_SIZE);
            //if(loop>0 && MOVE_TTL<prev_ttl && MOVE_TTL>0.999999*prev_ttl) myu->Set_Stepsize(myu->STEP_SIZE/2);
            //if(loop%1000 == 0) myu->Set_Stepsize(myu->STEP_SIZE/2);
            prev_ttl = MOVE_TTL;
            time(&end_iterate);
            dif_iterate = difftime(end_iterate, start_iterate);
            printf("Iterate Loop %i of %i took %.2lf seconds or %.2lf hours to run.\n", loop+1, tot_nloops, dif_iterate, dif_iterate/3600.0);
            /*
            if(myu->STEP_SIZE < 0.0000000000000001){
                printf("Move size=%.9f is below threshold.  Breaking\n", myu->STEP_SIZE); 
                tot_nloops = loop;
                break;
            }    
            */
            
            if(difftime(end_iterate, start) > hour_count*3600){
                printf("Writing Position Data after %i hour(s)\n", hour_count);
                hour_count++;
                Write_Positions();
                Write_Stats(myu->nmovies, myu->nusers, myu->nentries);
                
                FILE* loopstats = fopen( myu->loopstats_filename, "w");
                fprintf(loopstats, "Prev_Nloops: %i ", loop);
                fclose( loopstats );

                fclose( mfile ); //make sure it writes to disk
                mfile = fopen( myu->movefile_filename, "a");
                               
                fclose( ufile );
                ufile = fopen( myu->ufile_filename, "a");
                
                fclose( ssfile );
                ssfile = fopen( myu->ssizefile_filename, "a");
            }
        }
        fclose( mfile );
        fclose( ufile );
        fclose( ssfile );
        
        TCanvas cmove("cmove", "");
        cmove.Divide(1,3);
        cmove.cd(1);
        gPad->SetLogy(1);
        hmove.Draw();
        
        cmove.cd(2);
        //gPad->SetLogy(1);
        huplot.Draw("L");
        
        cmove.cd(3);
        gPad->SetLogy(1);
        hssize.Draw("L");
        
        cmove.SaveAs(myu->movettl_filename);
        
        FILE* loopstats = fopen( myu->loopstats_filename, "w");
        if( !loopstats ){
            printf("Couldn't open !!\n");
        }else{
            fprintf(loopstats, "Prev_Nloops: %i ", tot_nloops);
        }
        fclose( loopstats );

        
        
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
        Write_Stats(myu->nmovies, myu->nusers, myu->nentries);
    }
    
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

        Calc_Fit_Param();
        Write_Predictions(make_real, make_mine, 7);
        Float_t ans = Compare_Predictions();
        Float_t ans2 = newCompare_Predictions();
        printf(" old answer for rmse is %.4f\n", ans);
        printf(" new answer for rmse is %.4f\n", ans2);
        /*
        Float_t rmse_array[3];
        for(Int_t predict_mode=0; predict_mode<3; predict_mode++){
            if(predict_mode > 0) make_real = false;
            
            Write_Predictions(make_real, make_mine, predict_mode); 

            rmse_array[predict_mode] = Compare_Predictions();
            printf("mode: %i RMSE is %.4lf BEFORE\n",
                   predict_mode, rmse_array[predict_mode]);
            Optimize_RSCALE();
            printf("Going to redo Write_Predictions...\n");
            Write_Predictions(false, true, predict_mode); //Rewrite predictions using new RSCALE 
            rmse_array[predict_mode] = Compare_Predictions();
            printf("mode: %i RMSE is %.4lf AFTER\n",
                   predict_mode, rmse_array[predict_mode]);
        }
        for(Int_t predict_mode=0; predict_mode<3; predict_mode++){
            printf("mode: %i RMSE is %.4lf\n",
                   predict_mode, rmse_array[predict_mode]);
        }
        */
        time(&end_fn);
        printf("Predictions took %li seconds to run.\n", end_fn-start_fn);
    }
    //TTree* etree = new TTree("etree", "evolve tree");
    //Fill_Tree(myu, etree);
    //etree->BuildIndex("index");
    //Display(etree);

    Write_Stats(myu->nmovies, myu->nusers, myu->nentries);
    
    if(users_array) delete[] users_array;
    if( myu)  delete myu;

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
        movie_index = movie-1;
        temp = myu->recommend + movie_index;
        len = temp->length;
        
        for(Int_t j=0; j<len; j++){
            //Cory: What should I do here?
            dist2 = Distance2(myu->movies + movie_index, myu->users + temp->user_index[j]);//Cory: check 2nd
            Float_t d0 = Ideal(temp->rating[j]);
            if       (myu->step_mode == 0){ //Cory: Formula
                //prediction = 5.0*myu->MAXDIST/(4.0*dist+myu->MAXDIST); 
                U += -1.0*(temp->rating[j]-myu->RTHRESH)/sqrt(dist2);
            }else if (myu->step_mode == 1){
                //prediction = 4.0*(1-dist/myu->MAXDIST)+1.0;
                U +=  1.0*(temp->rating[j]-myu->RTHRESH)*sqrt(dist2);
            }else{
                //prediction = -4.0*pow(dist/myu->MAXDIST,2) + 5.0;
                //U +=  1.0*(temp->rating[j]-myu->RTHRESH)*dist2;
                U += pow( (sqrt(dist2) - d0),2);
                //u = 1/2 k (x-x0)2;
                //f = -dU = -k(x-x0)    
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
    Particle *temp_mov=NULL, *temp_usr=NULL;
    Int_t len, movie_index;
    Int_t mov_adr, usr_adr;
    Int_t nslots = myu->ndim * myu->nentries;
    
    
    Float_t* loop_sum = new Float_t[nslots];
    Float_t* move = new Float_t[myu->ndim];
        
    for(Int_t i=0; i<myu->ndim; i++){
        move[i]     = 0;    
    }
    for(Int_t i=0; i<nslots; i++){
        loop_sum[i] = 0;
    }
    
    for(Int_t movie=/*1*/special; movie<special+numspecial/*myu->nmovies+1*/; movie++){//Cory: Make upper limit special
        //printf("Iterating with movie:%i\n", movie);
        movie_index = movie - 1;
        temp = myu->recommend + movie_index;
        len = temp->length;
        mov_adr = myu->ndim * movie_index;
        
        for(Int_t i=0; i<len; i++){
            //for 2 processors
            //make increment +=2
            /*
            arg1.part1 = myu->movies + movie_index;
            arg1.part2 = myu->users + temp->user_index[i];
            arg1.rating = temp->rating[i];
            arg1.move = move1;
            TThread *th1 = new TThread("step1", Step, (void*) &arg1);

            if(i+1<len){
                arg2.part1 = myu->movies + movie_index;
                arg2.part2 = myu->users + temp->user_index[i+1];
                arg2.rating = temp->rating[i+1];
                arg2.move = move2;
                TThread *th2 = new TThread("step2", Step, (void*) &arg2);
            }

            th1->Run();
            if(i+1<len){
                th2->Run();
            }
            
            usr_adr = myu->ndim*(myu->nmovies + temp->user_index[i]);
            for(Int_t core=0; core<2; core++){
                if(core == 0) move = move1;
                else          move = move2;
                for(Int_t dim=0; dim<myu->ndim; dim++){
                    loop_sum[ mov_adr + dim ] += 0.5*move[dim];
                    loop_sum[ usr_adr + dim ] -= 0.5*move[dim];
                }
            }
            STEP_COUNT++;

            delete th1;
            delete th2;
            */
            
            Step(myu->movies + movie_index, myu->users + temp->user_index[i], temp->rating[i], move);

            usr_adr = myu->ndim*(myu->nmovies + temp->user_index[i]);
            for(Int_t dim=0; dim<myu->ndim; dim++){
                loop_sum[ mov_adr + dim ] += 0.5*move[dim];
                loop_sum[ usr_adr + dim ] -= 0.5*move[dim];
            } 
            STEP_COUNT++;
        }
    }
    for(Int_t movie=/*1*/special; movie<special+numspecial/*myu->nmovies+1*/; movie++){
        movie_index = movie - 1;
        temp_mov = myu->movies + movie_index;
        mov_adr = myu->ndim * movie_index;
        for(Int_t dim=0; dim<myu->ndim; dim++){
            temp_mov->x[dim] += loop_sum[mov_adr+dim];
            temp_mov->x[dim]  = OutofBounds( temp_mov->x[dim] );
            if(movie == 1) printf("dim: %02i loop_sum: %f \n", dim, loop_sum[mov_adr+dim]);
        }
    }
    for(Int_t user_index=0; user_index<myu->nusers; user_index++){
        temp_usr = myu->users + user_index;
        usr_adr = myu->ndim*( myu->nmovies + user_index);
        for(Int_t dim=0; dim<myu->ndim; dim++){
            temp_usr->x[dim] += loop_sum[usr_adr+dim];
            temp_usr->x[dim]  = OutofBounds( temp_usr->x[dim] );
            if(temp_usr->id == 1488844) printf("dim: %i loop_sum: %f \n", dim, loop_sum[usr_adr+dim]);
            //printf("loop_sum[i]: %f \n", loop_sum[i]);
        }
    }
    
    delete[] move;
    delete[] loop_sum;
}

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
   
    return dir;
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

Float_t
Ideal(Short_t rating){
    //return (myu->MAXDIST-2.0)*(5.0-rating)/4.0 + 1.0; //Formula
    //Cory:should this be max_dist??? or depend on it squared?
    return myu->MAXDIST*(5.0-rating)/4.0; //Formula
}

void
Step(Particle* part1, Particle* part2, Short_t rating, Float_t* move){
/*should two objects move closer or further based on rating
 * move = STEP_SIZE*(sign*mag)*(dist_x2/dist2)
 * should this be devided by 1/dist2 or 1/distx2
 * 
 * step_mode: 0=1/r
 *            1=random
 *            2=Hooke
 */
    Int_t sign;
    Float_t mag_move;
    Float_t* dist_x2 = new Float_t[myu->ndim];
    Float_t dist2 = Distance2(part1, part2, dist_x2);
    Float_t dist  = sqrt(dist2);
    Float_t scale=0;
    Float_t d0 = Ideal(rating);
    
    //Cory: don't know about this!
    if      (myu->step_mode == 0) scale = (myu->STEP_SIZE/dist2)/(sqrt(dist2)-d0);
    else if (myu->step_mode == 1) scale = (myu->STEP_SIZE/dist2); 
    else{
        scale = (myu->STEP_SIZE/dist2) * TMath::Abs(dist - d0) * (dist - d0);
    }
        
    if(part1->id == 1 && part2->id == 1488844){
        printf("Particle 1:");
        part1->Print(myu->ndim);
        printf("Particle 2:");
        part2->Print(myu->ndim);
        printf("d: %f\n", dist);
        printf("d0: %f\n", d0);
    }
    for(Int_t dim=0; dim<myu->ndim; dim++){
        sign    = Sign(part1, part2, rating, dim);
        mag_move = scale*dist_x2[dim];
        
        move[dim] = sign*mag_move;
        /*
        part1->x[dim] += 0.5*move[dim];  //Cory:!!!!!!!!!!!!!!!!!!!!!!!
        part2->x[dim] -= 0.5*move[dim];
        
        part1->x[dim] = OutofBounds( part1->x[dim] );
        part2->x[dim] = OutofBounds( part2->x[dim] );

        if(part1->id == 1 & part2->id == 1488844){
            printf("dim: %i Sign: %i mag_move: %f\n", dim, sign, mag_move);
        }
        */
        //if(STEP_COUNT %20000000 == 0) printf("d2=%.4lf dx2=%.4lf |move| is %.4lf\n", dist2, dist_x2[dim], TMath::Abs(move[dim])); 
        MOVE_TTL += TMath::Abs(move[dim]);
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
                /*
                if (mov == special){
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
                */
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
        printf("Num of predictions made %i\n", countrmse);
        printf("RMSE for me is %.4f\n", sqrt(sum2/countrmse));

        gStyle->SetStatY(gstyle_orig);
        hmovies->Draw();
        canvas1->Update();

        gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
        husers->Draw("sames");
        canvas1->Update();
        
        gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
        hspecial->Draw("sames");
            
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
            /*
            for(Int_t i=0; i<myu->nusers; i++){ //Cory: Make upper limit nusers + 1
                
                pos = myu->users[i].x;
              
                if(i<myu->nusers){
                    //husers->Fill(pos[0], pos[1]);
                }else{
                    //hmovies->Fill(pos[0], pos[1]);
                }
            }
            */
        }
        printf("RMSE for me is %.6f\n", sqrt(sum2/countrmse));

        gStyle->SetStatY(gstyle_orig);
        hmovies->Draw();
        //hmovies->Write("",TObject::kOverwrite);
        canvas1->Update();

        gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
        husers->Draw("sames");
        //husers->Write("",TObject::kOverwrite);
        canvas1->Update();
        
        gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
        hspecial->Draw("sames");
        //hspecial->Write("",TObject::kOverwrite);
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
Fill_Postree(TTree* ptree){
    /*
    Particle* ptemp;
    TBranch *busers  = ptree->Branch("particles" ,  ptemp, Form("index/I:id/I:x[%i]/F:is_user/O",myu->ndim) );
    //TBranch *bmovies = ptree->Branch("movies",  ptemp, Form("index/I:id/I:x[%i]/F:is_user/O",myu->ndim) );

    for(Int_t i=0; i<myu->nusers; i++){
        ptemp = myu->users + i;
        ptree->Fill();
    }
    printf("# entries is %li\n", (long Int_t)ptree->GetEntries() );
    
    for(Int_t i=0; i<myu->nmovies; i++){
        ptemp = myu->movies + i;
        ptree->Fill();
    }
    printf("# entries is %li\n", (long Int_t)ptree->GetEntries() );
    */
}

void
Fill_Rectree(TTree* rtree){
    /*
    Info* temp;
    Int_t length;
    Int_t movie;
    Int_t *user_id;
    Int_t *user_index;
    Short_t *rating;
    Int_t trial;
    Int_t test[5];
        
    rtree->Branch("movie", &movie, "movie/I");
    rtree->Branch("length", &length, "length/I");
    rtree->Branch("trial", &trial, "trial/I");

    rtree->Branch("test", test, "test[5]/I");
    rtree->Branch("user_id", user_id, "user_id[5]/I");
    rtree->Branch("user_index", user_index, "user_index[5]/I");
    rtree->Branch("rating", rating, "rating[5]/S");
    for(Int_t i=0; i<myu->nmovies; i++){
        movie = i+1;
        temp = myu->recommend + i;
        length = temp->length;

        //test = new Int_t[length];
        if(i < 100) for(Int_t j=0; j<5; j++) test[j] = temp->user_id[j];
        
        user_id = temp->user_id;
        user_index = temp->user_index;
        rating = temp->rating;
        if(i==50) trial = user_id[0];

        
        rtree->Fill();
        //delete test;
        if(i == 50) printf("%i %i %i %i %i\n", temp->user_id[0], temp->user_id[1], temp->user_id[2], temp->user_id[3], temp->user_id[4]);
        if(i == 50) printf("%i %i %i %i %i\n", user_id[0], user_id[1], user_id[2], user_id[3], user_id[4]);
                        
    }
    printf("Done Fill_Rectree\n");
    */
}

Int_t
Count_Movies(){
    Char_t* filename;
    Int_t count = 0;
    while(1){
        filename = Filename(count+1); //movies start at 1
        ifstream infile (filename);
        if(!infile) break;
        count++;
    }
    printf("%i unique movies found\n", count);
    return count;
}

void
Write_Stats(Int_t nmovies, Int_t nusers, Int_t nentries){
    FILE* outfile;
    outfile = fopen (myu->stats_filename, "w");
    fprintf(outfile, "nentries: %i\n", nentries);
    fprintf(outfile, "nmovies: %i\n", nmovies);
    fprintf(outfile, "nusers: %i\n", nusers);
    fclose (outfile);
}

void
Load_Stats(Int_t* nmovies, Int_t* nusers, Int_t* nentries){
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
        //printf("head:%x\n", head);
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
                if( !Add_Piece(user_id) ) users_entered++;
                //else printf("Couldn't add %i\n", id);
            }
            fprintf(outfile, "%i: %i\n", movie, count);
            if(i >= special && i < special+5) printf("movie:%i length:%i\n", i, count); 
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
Calc_Fit_Param(){

    Int_t movie_index, len;
    Float_t dist;

    Float_t beta_sum_num=0, beta_sum_den=0;
    //Float_t x_sum=0, y_sum=0;
    Float_t x_mean, y_mean;
    //Int_t count=0;
    
    Info* temp = NULL;
    TCanvas ccomp("ccomp", "Compare");
    TH2F hcomp("hcomp", "compare rating to dist", (Int_t)(2*ceil(myu->MAXDIST)), 0, (Int_t)ceil(myu->MAXDIST), 5, 1, 6);
    for(Int_t movie=1; movie<myu->nmovies+1; movie++){
        movie_index = movie-1;
        temp = myu->recommend + movie_index;
        len = temp->length;
        for(Int_t j=0; j<len; j++){
            if(temp->in_probe[j]){
                //count++;
                dist = Distance(myu->movies + movie_index, myu->users + temp->user_index[j], NULL);
                //x_sum += dist;
                //y_sum += temp->rating[j];
                
                hcomp.Fill(dist, temp->rating[j]);
            }
        }
    }
    //printf("count is %i\n", count);
    //printf("x_sum: %.4f y_sum: %.4f\n", x_sum, y_sum);
    //printf("my x_mean = %.4f, root's x_mean = %.4f\n", x_sum/count, hcomp.GetMean(1));
    //printf("my y_mean = %.4f, root's y_mean = %.4f\n", y_sum/count, hcomp.GetMean(2));
         
    x_mean = hcomp.GetMean(1);//x_sum / count;
    y_mean = hcomp.GetMean(2);//y_sum / count;

    for(Int_t movie=1; movie<myu->nmovies+1; movie++){
        movie_index = movie-1;
        temp = myu->recommend + movie_index;
        len = temp->length;
        for(Int_t j=0; j<len; j++){
            if(temp->in_probe){
                dist = Distance(myu->movies + movie_index, myu->users + temp->user_index[j], NULL);
                beta_sum_num += ( (dist - x_mean) * (temp->rating[j] - y_mean) );  
                beta_sum_den += pow(dist - x_mean, 2);
            }
        }
    }
    myu->BETA = beta_sum_num / beta_sum_den;
    myu->ALPHA = y_mean - myu->BETA*x_mean;

    printf("prediction = %.4f + %.4f*dist\n", myu->ALPHA, myu->BETA);

    hcomp.Draw("COLZ");
    ccomp.SaveAs(myu->compare_filename);
}

void
Write_Predictions(Bool_t real, Bool_t mine, Int_t predict_mode){
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
                        guess = Make_Prediction(movie, user_index, predict_mode);
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
Make_Prediction(Int_t movie, Int_t user_index, Int_t predict_mode){
    static Int_t count=0;
    Int_t movie_index = movie-1;
    Float_t dist = Distance(myu->movies + movie_index, myu->users + user_index);
    Float_t prediction=-1;
    if(predict_mode == 0){
        if       (myu->step_mode == 0){
            prediction = 5.0*myu->MAXDIST/(4.0*dist+myu->MAXDIST); //Formula 
        }else if (myu->step_mode == 1){
            prediction = 4.0*(1.0-dist/myu->MAXDIST)+1.0;
        }else{//HERE
            prediction = 5.0 - 4.0*(dist-1.0)/(myu->MAXDIST-2.0);//-5.0*pow(dist/myu->MAXDIST,2) + 5.5;
        }
    }else if(predict_mode == 1){
        prediction = 5.0 - 4.0*(dist)/(myu->MAXDIST);//Formula
    }else{
        prediction = 5.0 - 4.0*(dist-1.0)/(myu->MAXDIST-1.0);//Formula
        //prediction = 5.0 - 4.0*dist/myu->MAXDIST;
    }

    if(1){
        prediction =  myu->ALPHA + myu->BETA*dist;
        if(count < 10) printf("dist: %.4f   pred: %.4f\n", dist, prediction);
        count++;
        return prediction;
    }
    //printf("ANS = %.2lf\n", prediction);
    //Float_t prediction = 5.0*myu->RSCALE;
    //if(myu->RSCALE == 1.0) return round(myu->RSCALE*prediction);
    Float_t m=1.0;
    return myu->RSCALE*(m*prediction + (1.0-m)*round(prediction));
}

Float_t
newCompare_Predictions(){
    Float_t prediction, dist, diff, rmse;
    Float_t sum2=0;
    Int_t count=0, len, movie_index;
    Info* temp=NULL;
    
    for(Int_t movie=1; movie<myu->nmovies+1; movie++){
        movie_index = movie-1;
        temp = myu->recommend + movie_index;
        len = temp->length;
        for(Int_t j=0; j<len; j++){
            if(temp->in_probe[j]){
                count++;
                dist = Distance(myu->movies + movie_index, myu->users + temp->user_index[j], NULL);
                prediction = myu->ALPHA + myu->BETA*dist;
                
                diff = prediction - temp->rating[j];
                
                if( TMath::Abs(diff) > 1.5 ) printf("movie:%i user_id:%i dist:%.4f pred:%.4f real:%i\n", movie, temp->user_id[j], dist, prediction, temp->rating[j]);

                sum2 += pow( diff, 2);
            }
        }
    }
    rmse = sqrt ( sum2 / count );
    return rmse;
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
            TH2F hreg("hreg", "Real vs Mine", 100, 0, 6, 5, 1, 6);
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

            TCanvas c_recommend("c_recommend", "Recommend Plot");
            hreal.SetLineColor(kGreen);
            hmine.SetLineColor(kRed);

            gStyle->SetStatY(gstyle_orig);
            hreal.Draw();
            //hreal.Write("",TObject::kOverwrite);
            c_recommend.Update();
            gStyle->SetStatY(gStyle->GetStatY() - gStyle->GetStatH());
            hmine.Draw("sames");
            //hmine.Write("",TObject::kOverwrite);
            c_recommend.SaveAs(myu->recdist_filename);

            gStyle->SetStatY(gstyle_orig);
            TCanvas c_reg("c_reg", "Regression Plot");
            hreg.Draw("COLZ");
            //hreg.Write("",TObject::kOverwrite);
            c_reg.SaveAs(myu->reg_filename);
            
            FILE* outfile;
            outfile = fopen (myu->rmse_filename, "w");
            fprintf(outfile, "RMSE: %.4lf\n", rmse);
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
    Bool_t in_probe;
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

    for(Int_t movie=special; movie<special+numspecial;movie++){ //Cory: Changed this too
        //open file 'movie'
        //movie_index = (movie-1);
        in_filename = Filename(movie);
        ifstream infile (in_filename);
        if (infile.is_open()){
            count = 0;
            if(movie%5000==0) printf("movie: %i\n", movie);
            if(movie >= special && movie < special+5) printf("movie: %i length is %i\n", movie, num_predict[movie-1]);
           
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
                //in_probe = Is_In_Probe(movie, user_id);
                myu->recommend[movie-1].user_id   [count] = user_id;
                myu->recommend[movie-1].rating    [count] = rating;
                myu->recommend[movie-1].user_index[count] = user_index;
                //myu->recommend[movie-1].in_probe  [count] = in_probe;
                //if((movie == 5 && in_probe) || user_id == 1488844) printf("movie:%i user_id:%i rating:%i in_probe:%i\n",
                //                              movie, user_id, rating, in_probe);
                count++;
                infile.ignore(256, '\n');
            }
                    
            infile.close();
        }else printf("Unable to open file %s\n", in_filename);
    }

    Int_t input;
    ifstream probefile (myu->probe_filename);
    Info* temp=NULL;
    if (probefile.is_open()){
        while(!probefile.eof()){
            while ( probefile>>input  && probefile.peek() != ':'){
                //printf("input:%i\n", input);
                Int_t len = temp->length;
                for(Int_t j=0; j<len; j++){
                    if( temp->user_id[j] == input){
                        temp->in_probe[j] = true;
                        break;
                    }
                }
            }
            if (probefile.peek() == ':'){
                //printf("movie:%i user_id: %i     input:%i\n", movie, user_id, input);
                temp = myu->recommend + input - 1;
                probefile.ignore(256, '\n');
            }
        }
    }else printf("Couldn't open %s file\n", myu->probe_filename);
      
    delete[] num_predict;
    printf("\tLeaving Recommend_Tables\n");
    return 0;
}

Bool_t
Is_In_Probe(Int_t movie, Int_t user_id){
    Bool_t right_movie = false;
    Int_t input;
    ifstream probefile (myu->probe_filename);
    if (probefile.is_open()){
        while(!probefile.eof()){
            while ( probefile>>input  && probefile.peek() != ':'){
                //printf("input:%i\n", input);
                if(right_movie && input == user_id) return true;
            }
            if (probefile.peek() == ':'){
                //printf("movie:%i user_id: %i     input:%i\n", movie, user_id, input);
                if(right_movie) return false;
                if( input == movie ) right_movie = true;
                probefile.ignore(256, '\n');
            }
        }
    }else printf("Couldn't open %s file\n", myu->probe_filename);
    return false;
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
    time_t start1, start2, end1, end2;
    Double_t dif;
    Char_t c;
    Int_t ivalue;
    Float_t fvalue;
    Char_t* pos_filename = ("./output/positions/positions_2_10.txt");
    Int_t obj_entered;
    Char_t line[25];

    time(&start2);
    obj_entered = 0;
    ifstream infile2 (pos_filename);
    if (infile2.is_open()){
        //Particle* temp=myu->part;
        while ( !infile2.eof()){
            if(obj_entered%100000 == 0) printf("obj_entered:%i\n", obj_entered);
            infile2.getline(line, 20, ','); //get index 
            ivalue = atoi(line);

            infile2.getline(line, 20, ','); //get ID
            ivalue = atoi(line);
            
            for(Int_t dim=0; dim<10; dim++){
                infile2.getline(line, 20, ','); //get x[dim]
                fvalue = atof(line); //needs to be a float
            }
            infile2.getline(line, 20, '\n'); //get ID
            ivalue = atoi(line);
            
            obj_entered++;
        
        }
        infile2.close();
    }
    else printf("Unable to open file %s\n", pos_filename);

    time(&end2);
    dif = difftime(end2, start2);
    printf("Program took %.2lf seconds or %.2lf hours to run.\n", dif, dif/3600.0);


    time(&start1);
    obj_entered = 0;
    ifstream infile (pos_filename);
    if (infile.is_open()){
        while ( infile>>ivalue>>c>>ivalue>>c){
            //infile.ignore(1);
            //infile>>ivalue;
            //infile.ignore(1);
            //index, id, x[0], ..., x[ndim-1], is_user 
            if(obj_entered%100000 == 0) printf("obj_entered:%i\n", obj_entered);
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
    else printf("Unable to open file %s\n", pos_filename);
    time(&end1);
    dif = difftime(end1, start1);
    printf("Program took %.2lf seconds or %.2lf hours to run.\n", dif, dif/3600.0);
    */
    /*
    Universe* temp_u = new Universe();
    printf("%.2f %s\n", temp_u->RTHRESH, temp_u->pos_filename);
    delete temp_u;
    */
    //printf("%f, %f %f %f\n", round(5.2), round(5.5), round(5.9), round(6.0));    
    /*
    Char_t x = 5;
    printf("as a char:%c as a int:%i with size:%i\n", x,x, sizeof(x) );
    Short_t y = 5;
    printf("as a char:%c as a int:%i with size:%i\n", y,y, sizeof(y) );
    Int_t z = 5;
    printf("as a char:%c as a int:%i with size:%i\n", z,z, sizeof(z) );

    Info a(400);
    //printf("info size:%i len:%i k:%i l:%i\n", sizeof(a), a.length, sizeof(a.k), sizeof(a.l));
    */
    /*
    Char_t* old_name = Form("%s%i", "000000", 1);
    Char_t* new_name = Form("%07i", 1);
    printf("%s\n", old_name);
    printf("%s\n", new_name);
    */
    /*
//    Fill_Rectree()
    TFile* recfile = new TFile("recommendations.root","update");
    myu = new Universe();
    Int_t step_mode = 2, ndim = 2;
    Int_t nmovies, nusers, nentries;
    Float_t step_size=0.1;
    TTree* pos_tree= new TTree("pos_tree", "");
    TTree* rec_tree= new TTree("rec_tree", "");
    
    Load_Stats(&nmovies, &nusers, &nentries);
    myu->Set_Prop(nmovies, nusers, step_mode, ndim);
    myu->Set_Stepsize(step_size);
    Load_Userlist(nusers);
    Load_Positions();
    Recommend_Tables();
    Initialize(); 
    Fill_Rectree(rec_tree);
    if (1){
        TH1F hpx("hpx","This is the px distribution",100,-50,50);
        //TH2F *hpxpy = new TH2F("hpxpy","py ps px",40,-4,4,40,-4,4);
        for(Int_t i=0; i<1000; i++) hpx.Fill( gRandom->Gaus(50,10) );
        hpx.Write("",TObject::kOverwrite);
        //delete hpx;
    }
    
    recfile->Write("",TObject::kOverwrite);
    //recfile->Write();

   

    recfile->Close();

    delete recfile;
    delete myu;
    delete rec_tree;
    delete pos_tree;
    */
    //for 2 processors
    //make increment +=2

    //check();
    return;
}
/*
Double_t gPi = 0.0;
Double_t gIntervals = 1.0;
void* process(void *arg);

void
check(){
    gSystem->Load("CalcPiThread.so");
    
    // create the TThread instances
    
    gIntervals = 10;
    printf("Intervals %f\n", gIntervals);
    
    TThread *th1 = new TThread("pi1",process,(void*)"0");
    TThread *th2 = new TThread("pi2",process,(void*)"1");
    
    th1->Run();
    th2->Run();
}

void *process(void *arg)
{
    Double_t width, localsum;
    Int_t i;
    Int_t iproc = (*((char *) arg) - '0');
    
    printf("Intervals %f\n", gIntervals);
    
    // Set Width
    width = 1.0/gIntervals;
    
    // Do the local computations 
    
    localsum = 0.;
    for(i = iproc; i < gIntervals; i += 2)
    {
        Double_t x = (i + 0.5) * width;
        localsum += 4.0 / (1.0 + x * x);
    }
    
    localsum *= width;
    
    TThread::Lock();
    gPi += localsum;
    printf("Current Estimation of pi is %20.18f\n", gPi);
    TThread::UnLock();
    
    return (0);
}
*/
