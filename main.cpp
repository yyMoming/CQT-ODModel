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
// #include "opencv2/opencv.hpp"
#include "Model.h"
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

void testMusicLoader()
{
    int times = 1000;
    double sum = 0.0;
    char *filename = "../9323.mp3";
    char *outfile = "./outfile.txt";
    vector<double> amp_y;

    for (int i = 0; i < times; i++)
    {
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

// void imshowSpectrogram(){
//     vector<double> amp_y;
//     char *filename = "../9323.mp3";
//     mp3_loader(filename, amp_y);
//     int length = amp_y.size();
//     cv::Mat canvas(800, 100, CV_32FC3, cv::Scalar(255, 255, 255));
//     if(length < 1){
//         cout << "图像数据不足" << endl;
//         return;
//     }

//     for (int i = 1; i < length; i++){
//         cv::Point2i p1(i - 1, amp_y[i - 1]);
//         cv::Point2i p2(i, amp_y[i]);
//         cv::line(canvas, p1, p2, cv::Scalar(0, 0, 255), 1, cv::LINE_8, 0);
//     }
//     cv::namedWindow("img", cv::WINDOW_FREERATIO);
//     cv::imshow("img", canvas);
//     cv::waitKey(0);

// }
void testModel()
{
    CQT::_m_cqt<double> mc(44100, 27.5, 36, 267, 512);
    int length = countlines("./output.txt");
    double *y;
    y = readlines("./output.txt", y);
    mc.calc(y, length, 44100);
    auto res_frames = mc.cqt_res();
    CQT::_m_OnsetDetector mOD;
    mOD.load_model("/home/ywm/MUSIC/Codes/onset_train/model.pt", "cpu");
    vector<float> res;
    mOD.predict(res_frames, mc.get_n_frames(), res);
    // CQT::_m_input(res_frames, mc.get_n_frames(), res);
    fstream ofs("predict.txt", ios::out);
    for (auto x : res)
    {
        ofs << x << endl;
    }
    ofs.close();
}

int main()
{
    testModel();
    // testMusicLoader();
    // testCQT();
}