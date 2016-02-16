class Piece{
public:
    Int_t id;
    Piece* less;
    Piece* more;

    Piece();
    Piece(Int_t);
    ~Piece();
};
Piece::Piece(){
    id = -1;
    less = NULL;
    more = NULL;
}
Piece::Piece(Int_t input){
    id = input;
    less = NULL;
    more = NULL;
}
Piece::~Piece(){
    id = -1;
    delete this->less;// = NULL;
    delete this->more;// = NULL;
}
