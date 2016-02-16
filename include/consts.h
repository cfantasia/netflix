Int_t NMOVIES=17770;
const Int_t NUSERS_GIVEN=480189;
Int_t NUSERS=0;
const Int_t NRATINGS=5;
Int_t nratings=NMOVIES+NUSERS;

Int_t NDIM = 11;
const Int_t NPROP = 3;
const Char_t* PROP_NAME[NPROP] = {"E","Njets", "H"};
Int_t NENTRIES = NMOVIES + NUSERS; //500;
const Float_t MINX = -10.0;
const Float_t MAXX = 10.0;

const Float_t STEP_SIZE = 0.01;
const Float_t RTHRESH=3.0;
Float_t MAXDIST = sqrt(NDIM)*(MAXX-MINX)/2;
Float_t RSCALE = 1.0;
Float_t MOVE_TTL = 0;
Int_t STEP_MODE = 2;

Int_t STEP_COUNT;

Char_t* posbase_filename = ("/home/fantasia/work/netflix/output/positions/positions");
Char_t* pos_filename;

const Char_t* user_filename =("/home/fantasia/work/netflix/output/user_list.txt");
const Char_t* predict_filename=("/home/fantasia/work/netflix/output/my_predictions.txt");
const Char_t* rating_filename =("/home/fantasia/work/netflix/output/real_ratings.txt");
const Char_t* numpredict_filename=("/home/fantasia/work/netflix/output/numpredict.txt");
const Char_t* rmse_filename=("/home/fantasia/work/netflix/output/rmse.txt");
const Char_t* pic_filename =("/home/fantasia/work/netflix/output/Universe_Image.pic");
const Char_t* recdist_filename =("/home/fantasia/work/netflix/output/Recommend_Distribution.pic");
const Char_t* reg_filename =("/home/fantasia/work/netflix/output/Regression.pic");
const Char_t* stats_filename =("/home/fantasia/work/netflix/output/stats.txt");

const Char_t* probe_filename=("/home/fantasia/work/netflix/data/probe.txt");
const Char_t* qual_filename =("/home/fantasia/work/netflix/data/qualifying.txt");

const Char_t* EMPTYLINE =("0000000000000000000000000");
