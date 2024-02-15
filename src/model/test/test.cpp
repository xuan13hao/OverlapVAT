#include <iostream>
#include "../minimap2.h"
//g++ -g -O2 test.cpp ../minimap2/libminimap2.a -lz -lm -o minimap-lite
using namespace std;
using std::cout;
using std::endl;
int main()
{

    string ref = "../../../data/ref.fa";
    string qry = "../../../data/query.fa";
    CVector* cv = DNA_MM_module(ref.c_str(),qry.c_str());
    int size = get_size(cv);
    cout<<size<<endl;
    for (size_t i = 0; i < size; i++)
    {
        VATMappedRec mr_ =  get_element(cv, i);
        cout<<mr_.qry_name<<"\t"<<mr_.ref_name<<"\t"<<mr_.strand<<"\t"<<mr_.cigar_<<endl;
    }
    destroy_vector(cv);
    return 0;
}