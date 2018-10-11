#ifndef DEF_SIMPCELL
#define DEF_SIMPCELL

struct simpCell{
    int k, tum, ves;
    double r, rr;
};

inline bool compK (simpCell a, simpCell b){
    return (a.k < b.k);
}


inline bool compR (simpCell a, simpCell b){
    return (a.r < b.r);
}


inline bool compRR (simpCell a, simpCell b){
    return (a.rr < b.rr);
}

#endif
