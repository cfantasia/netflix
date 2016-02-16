class Info{
public:            
    Int_t* user_id;
    Int_t* user_index;
    Short_t* rating;
    Bool_t* in_probe;
    Int_t length;
           
    Info();
    Info(Int_t);
    ~Info();

    void Prep_Info(Int_t);
};

Info::Info(){
    length = 0;
    user_id = NULL;
    user_index = NULL;
    rating = NULL;
    in_probe = NULL;
}

Info::Info(Int_t len){
    length = len;
    user_id = new Int_t[len];
    user_index = new Int_t[len];
    rating = new Short_t[len];
    in_probe = new Bool_t[len];
}

Info::~Info(){
    length = 0;
    delete user_id;
    delete user_index;
    delete rating;
    delete in_probe;
}

void
Info::Prep_Info(Int_t len){
    length = len;
    user_id = new Int_t[len];
    user_index = new Int_t[len];
    rating = new Short_t[len];
    in_probe = new Bool_t[len];
}
