#include <iostream>
#include <string>
#include <fstream>
void text_combine()
{
    int mod, hand, ch, status, a0, b; 
    double Q1, a00, a1, a5;
    ifstream f1 ("gaindata_t1.txt");
    ifstream f2 ("gain_par.txt");
    ofstream ff ("gain_data.txt");
    while((f1 >> mod >> hand >> ch >> Q1 >> status )&&(f2 >> a0 >> a00 >> a1 >> a5 >> b))
    {
        ff << mod <<"  "<< hand <<"  "<<  ch <<"  "<<  Q1 <<"  "<<  status <<"  "<<  a0 <<"  "<<  a00 <<"  "<<  a1 <<"  "<<  a5 <<"  "<<  b <<endl;
    }
    f1.close();
    f2.close();
    ff.close();
}