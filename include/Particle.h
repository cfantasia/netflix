class Particle{
  public:
    Int_t index;
    Int_t id;
    Float_t* x;
    //Float_t prop[NPROP];
    //Float_t rating;
    //Float_t U;
    Bool_t is_user;

    Particle();
    Particle(Int_t);
    //Particle(Int_t, Int_t, Float_t, Float_t, Int_t, Float_t, Bool_t);
    ~Particle();

    void Set_Ndim(Int_t);
    void Print(Int_t);
};

Particle::Particle(){
    index = -1;
    id = -1;
    //x = new Float_t[NDIM];
    //for(Int_t dim=0; dim<NDIM; dim++) x[dim]=0;
    //for(Int_t i=0; i<NPROP; i++) prop[i] = 0;
    //rating = -1;
    //U = 0;
    is_user = false;
}

Particle::Particle(Int_t ndim){
    index = -1;
    id = -1;
    x = new Float_t[ndim];
    for(Int_t dim=0; dim<ndim; dim++) x[dim]=0;
    is_user = false;
}
/*
Particle::Particle(Int_t in_index, Int_t in_id, Float_t in_x, Float_t in_prop,
                   Int_t in_rating, Float_t in_U, Bool_t in_is_user){
    index = in_index;
    id = in_id;
    for(Int_t dim=0; dim<NDIM; dim++) x[dim]=in_x[dim];
    for(Int_t i=0; i<NPROP; i++) prop[i] = in_prop[i];
    rating = in_rating;
    U = in_U;
    is_user = in_is_user;
}
*/
void Particle::Set_Ndim(Int_t ndim){
    x = new Float_t[ndim];
}

Particle::~Particle(){
    if(x) delete[] x;
}

void Particle::Print(Int_t ndim){
    printf("------------\n");
    printf("index: %i\n", index);
    printf("id: %i\n", id);
    printf("x:");
    for(Int_t dim=0; dim<ndim; dim++) printf("  %.4f\n", x[dim]);
    printf("is_user: %i\n", is_user);
}

