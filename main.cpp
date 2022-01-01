//
// Created by ywm on 2021/6/27.
//
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <complex>
#include "music_loader.h"
#include "mfshs.h"
#include "_m_cqt.h"
#include <time.h>
using namespace std;

template <class In, class In2, class T>
T _m_inner_product(In in1, In in2, In2 in3, T init)
{
    for (; in1 != in2; in1++, in3++)
    {
        init = init + (*in1) * (*in3);
    }
    return init;
}

void testCQT()
{
    CQT::_m_cqt<double> mc(44100, 27.5, 36, 267, 512);
    int length = countlines("./output.txt");
    double *y;
    y = readlines("./output.txt", y);
    int times = 1000;
    for (int i = 0; i < times; i++)
    {
        time_t start = clock();
        mc.calc(y, length, 44100);
        time_t endd = clock();
        cout << "[" << i << "] 耗时：" << double(endd - start) / CLOCKS_PER_SEC << endl;
    }
    delete[] y;
}

void testMusicLoader(){
    int times = 1000;
    double sum = 0.0;
    char *filename = "../9323.mp3";
    char *outfile = "./outfile.txt";
    vector<double> amp_y;
    for (int i = 0; i < times; i++){
        time_t start = clock();
        
        mp3_loader(filename, amp_y);
        time_t endd = clock();
        sum += double(endd - start) / CLOCKS_PER_SEC;
        cout << "[" << i << "]" << double(endd - start) / CLOCKS_PER_SEC << endl;

    }
    // fstream ofs(outfile, ios::out);
    // for(auto x: amp_y){
    //     ofs << x << endl;
    // }
    // ofs.close();
    cout << "avg cost:" << sum / times << endl;
}

int main()
{
    testMusicLoader();
    // testCQT();
}