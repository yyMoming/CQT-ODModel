//
// Created by yyMoming on 2021/5/10.
//

#ifndef SIGHTSING_MFSHS_H
#define SIGHTSING_MFSHS_H
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string.h>

#include "utils.h"
using namespace std;

class mfshs{
public:

    int sampleRate = 44100;
    float *pitch = nullptr;
    float *energe = nullptr;
    float *data = nullptr;
    int fftLength = 8192;
    int windowLength = 2048;
    int hopSize = 512;
    int frameSize = 2048;
    float H = 5;
    float h = 0.8;
    int n_frames;
private:
    float *hammingwindow = nullptr;
    float *xFrame = nullptr;
    void calcHammingwindow(){
        hammingwindow = new float[windowLength];
        for(int i = 0; i < windowLength; i++){
            hammingwindow[i] = float(0.54 - 0.46 * cos(2 * M_PI * i / (windowLength - 1)));
        }
    }
public:
    /*
    * 构造函数，
    * param：wav_file  读取的wav文件名
   * */
    mfshs(const string);
    ~mfshs();
    /*
     * 多线程中的一个线程函数
     * */
    void run(const int, const int);

    /*
     * 获取一帧的note
     * */
    float getNode(float *);
    /*
     * 计算每一帧的MFSHS音高
     * param: frameData: 一帧长度的数据样本点
     * */
    float calculatePitch(const float* );

    /*
     * calculateMFSHSPitch
     * 计算MFSHS音高
     * param: float *fftResult, fftResult计算结果
     * param: int maxResultIndex, fftResult中的最大值的索引
     * */
    float calculateMFSHSPitch(float *, int );

};

#endif //SIGHTSING_MFSHS_H
