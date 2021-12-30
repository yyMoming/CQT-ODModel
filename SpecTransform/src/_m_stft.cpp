//
// Created by yyMoming on 2021/7/8.
//

#include <complex>
#include <vector>
#include <fstream>
// #include "torch/torch.h"
#include "_m_cqt.h"
using namespace std;
const int countlines(const string &filename) {
    fstream fs(filename);
    if (!fs.is_open()) {
        cout << "此文件不存在！" << endl;
        pthread_exit(NULL);
    }
    string num;
    int lines = 0;
    while (getline(fs, num, '\n')) {
        lines++;
    }
    fs.close();
    return lines;
}

//
//int main(){
//    CQT::_m_cqt<double> mc(44100, 27.5, 36, 267);
//    int length = countlines("../../build/output.txt");
//    double *y = new double[length];
//
//    fstream ifs("../../build/output.txt", ios::in);
//    int i = 0;
//    while(ifs.good()){
//        ifs >> y[i];
//    }
//    ifs.close();
//    cout << length << endl;
//    double *res = new double[length / 2 + 1];
//    mc.resample_f(y, length, res);
//    vector<double> vec(res, res + length / 2 + 1);
//    cout << accumulate(vec.begin(), vec.end(), 0.0) <<endl;
//    cout << length << endl;
////    int a[9] = {1, 2, 3,4,5,6,7,8,9};
////    int *b = new int[20];
////    for(int i = 0; i < 20; i++){
////        b[i] = i;
////    }
////    vector<int> vec = {1,2,3,4,5,6};
////
////    torch::Tensor tmp1 = torch::tensor(vec);
////    cout << tmp1 << endl;
////    auto tmp2 = torch::from_blob(b, {2,10}, torch::kInt);
////    cout << tmp2 << endl;
//}

