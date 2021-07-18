//
// Created by yyMoming on 2021/5/19.
//

#ifndef SIGHTSING_UTILS_H
#define SIGHTSING_UTILS_H
#include <fstream>
#include <unistd.h>
#include <iostream>
using namespace std;
const int countlines(const string &filename);

template <class C>
C* readlines(const string &filename, C*){
    fstream fs(filename, ios::in);
    if(!fs.is_open()){
        cout << "文件【" << filename << "】 "<< "打开失败" << endl;
        pthread_exit(NULL);
    }
    int lines = countlines(filename);
    C* res = new C[lines];
//    C temp;
//    int i = 0;
//    while(fs.good()){
//        fs >> temp;
//        res[i] = temp;
//        i++;
//    }
    for(int i = 0; i < lines; i++){
        fs >> res[i];
    }
    fs.close();
    return res;
}


#endif //SIGHTSING_UTILS_H
