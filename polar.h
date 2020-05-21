#ifndef POLAR_H
#define POLAR_H

#include <vector>
#include <stack>
#include <math.h>

using namespace std;


extern int* bitreversal;
extern int* information_location;
extern int* frozen_location;
void readtxt(int block_length, int info_length, int frozen_length);
int sgn(double x);
vector<int> bpsk(vector<int> code_bits);
vector<int> crcgen(vector<int> data);
extern vector<int> crc;


class PolarCode
{
public:

    PolarCode(int n, int block_length, int info_length, int frozen_length, int crc_length) :
        _n(n)
    {
        _info_length = info_length;
        _block_length = block_length;
        _frozen_length = frozen_length;
        _crc_length = crc_length;
    }


    vector<int> encode(vector<int> block_bits);
    vector<int> decode_scl_llr(std::vector<double> llr, int list_size);

    void frozen();
    vector<int> _frozen_bits;


private:



    int _n;
    int _info_length;
    int _block_length;
    int _frozen_length;
    int _crc_length;


    //    vector<int> _frozen_bits;
    vector<int> _channel_order_descending;

    vector<double> _pathMetric_LLR;
    vector<vector<double*> > _arrayPointer_LLR;

    int _list_size;

    stack<int> _inactivePathIndices;
    vector<int > _activePath;
    vector<vector<int*> > _arrayPointer_C;
    vector<int*> _arrayPointer_Info;
    vector<vector<int> > _pathIndexToArrayIndex;
    vector<stack<int> > _inactiveArrayIndices;
    vector<vector<int> > _arrayReferenceCount;

    void initializeDataStructures();
    int assignInitialPath();
    int clonePath(int l);
    void killPath(int l);

    double* getArrayPointer_LLR(int lambda, int  l);
    int* getArrayPointer_C(int lambda, int  l);

    void recursivelyCalcLLR(int lambda, int phi);
    void recursivelyUpdateC(int lambda, int phi);

    void continuePaths_FrozenBit(int phi);
    void continuePaths_UnfrozenBit(int phi);
    vector<int> channel_order();


    //
    int findMostProbablePath(bool crc_check);

    bool crccheck(int* crc_data);





};

#endif // POLAR_H
