//
// Created by yyMoming on 2021/5/20.
//
#include "utils.h"

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