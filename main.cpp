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

void music_decoder(){
    time_t start, end;
    char* filename = "../9323.mp3";
    char* outfilename = "../9323.txt";
    time(&start);
    vector<float> amp_y;
    mp3_loader(filename, amp_y);
//    mp3_loader(filename, outfilename);
//    audio_decode_example("audio_decode.txt", filename);
//    fstream ifs("amp_y.txt", ios::out);
//    pcm_loader(ifs, amp_y);
    time(&end);
    cout << "总耗时" <<setprecision(20) <<  difftime(end, start) << endl;
    cout << "successfully!" << endl;
}

//void mfshs_use(){
//    mfshs mf(outfilename);
//    mf.run(0, 100);
//    cout << ywm2 << endl;
//    cout << ywm_define << endl;
//}
template<class In, class In2, class T>
T _m_inner_product(In in1, In in2, In2 in3, T init){
    for(; in1 != in2; in1++, in3 ++){
        init = init + (*in1) * (*in3);
    }
    return init;
}
//// 测试STFT
//int main(){
//    const double fs=44100;
//    const double f=1000;
//    const double duration=30;
//    const int length = fs*duration;
//    double * signal = (double*)malloc(sizeof(double)*(length));
//    for(int t=0;t<length;++t){
////        signal[t]=sin(2*PI*f/fs*t);
//        signal[t]=t%1024;
//    }
//    ofstream ofs;
//    time_t start = clock();
//    CQT::_m_stft(signal, length, 2048, 512);
//    time_t endd = clock();
//    cout << double(endd - start) / CLOCKS_PER_SEC << endl;
//    delete[] signal;
//}


int main(){
//    vector<complex<double>> vec(1000, {1,2});
//    vector<complex<double>> vec2(1000, {2,3});
////    complex<double> res = inner_product<complex<double>>(vec.begin(), vec.end(), vec2.begin(), complex<double>(0,0));
//    time_t start = clock();
//
//    auto res = _m_inner_product(vec.begin(), vec.end(), vec2.begin(), complex<double>(0,0));
//    time_t end = clock();
//    cout << fixed << setprecision(10) << double(end - start) / CLOCKS_PER_SEC << endl;
//    cout << res << endl;
    CQT::_m_cqt<double> mc(44100, 27.5, 36, 267, 512);
    int length = countlines("./output.txt");
    double *y;
    y = readlines("./output.txt", y);
    mc.calc(y, length, 44100);
    auto res_frames = mc.res_frames;
    int n_frames = mc.get_n_frames();
    fstream res_frames_imag("cqt_out_imag.txt", ios::out);
    fstream res_frames_real("cqt_out_real.txt", ios::out);
    for (auto x: res_frames){
        for (int i = 0; i < n_frames; i++, x++){
            if (i != n_frames - 1){
                res_frames_imag << x->imag() << ",";
                res_frames_real << x->real() << ",";
            }else{
                res_frames_imag << x->imag() << endl;
                res_frames_real << x->real() << endl;
            }
        }
    }
    res_frames_imag.close();
    res_frames_real.close();
    //    CQT::_m_stft(y, length, 512, 512);
    //    double *res = new double[length / 2];
    //
    //    mc.resample_f(y, length, res);
    //    auto sums = accumulate(res, res + length / 2, 0.0);
    //    cout << sums << endl;
    //    float sr=22050, fmin=2349.32;
    //    time_t start = clock();
    //    mc.constant_q(sr, fmin, 36);
    //
    //    mc.cqt_filter_fft(sr, fmin);
    //    time_t _end = clock();
    //    cout << "constant_q 耗时：" << double(_end - start) / CLOCKS_PER_SEC << endl;
    //
    //    start = clock();
    //    int new_length = length / 2;
    //    int hop_length = 256;
    //    mc.cqt_response(res, new_length, hop_length);
    //    _end = clock();
    //    cout << "cqt_response 耗时：" << double(_end - start) / CLOCKS_PER_SEC << endl;
    delete[] y;
}