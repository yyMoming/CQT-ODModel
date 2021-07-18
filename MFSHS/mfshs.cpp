//
// Created by yyMoming on 2021/5/9.
//
#include "mfshs.h"
#ifdef __cplusplus
extern "C"
{
#endif

#include "fft.h"

#ifdef __cplusplus
};
#endif
int N;
mfshs::mfshs(const string wav_file){
//        fstream fs(wav_file, ios::in|ios::binary);
//        pcm_loader(fs);
    string data_file = wav_file;
    fstream df(data_file, ios::in | ios::binary);
//    df.seekg(0, df.end);
    int length = countlines(data_file);
//    int length = 2000;
    cout << df.tellg() << endl;
    int padSize = this->frameSize / 2;
    df.seekg(0, df.beg);
    cout << df.tellg() << endl;
    this->data = new float[length + 2 * padSize];
    n_frames = floor(length / this->hopSize) + 1;
    this->pitch = new float[n_frames];
    this->energe = new float[n_frames];
    xFrame = new float[this->n_frames * this->frameSize];
    
    float sample = 0.0;
    int i = padSize;
    while(df.good()){
        df >> sample;
        data[i] = sample;
        i++;
    }
    for(i = 0; i < this->n_frames; i++){
        vector<float> tmp(this->data + i * this->hopSize, this->data + i * this->hopSize + this->frameSize);
//            cout << "data (" << i << ")" << accumulate(tmp.begin(), tmp.end(), 0.0) << endl;
        memcpy(this->xFrame + i * frameSize, this->data + i * this->hopSize, this->frameSize * sizeof(float));
    }
    cout << (float*)(this->xFrame + (i - 1) * frameSize) << endl;
    N = this->fftLength;
    calcHammingwindow();
    delete[] data;

}
mfshs::~mfshs(){
    if(pitch)
        delete[] pitch;
    if(energe)
        delete[] energe;
    if(hammingwindow)
        delete[] hammingwindow;
    if(xFrame)
        delete[] xFrame;
//        delete pitch;
//        delete energe;
//        delete hammingwindow;
//        delete xFrame;
}
/*
 * 多线程中的一个线程函数
 * */
void mfshs::run(const int start, const int end){
    ofstream ofs("note.txt", ios::out);
    for(int i = start; i < end; i++){
        vector<float> oneFrame(this->xFrame + i * this->frameSize, this->xFrame + (i + 1) * this->frameSize);

        float sampleSum = 0.0;
        for(auto x: oneFrame){
            sampleSum += abs(x);
        }
        float meanAmp = sampleSum / this->frameSize;
        float note = 0.0;
        if(meanAmp > 0.0001)
            note = this->getNode(this->xFrame + i * this->frameSize);
        this->pitch[i] = note;
        this->energe[i] = meanAmp;
//        cout <<"[" << i << "]"<< note << endl;
        ofs << note << endl;

    }
    ofs.close();
}
/*
 * 获取一帧的note
 * */
float mfshs::getNode(float *data){
    float fpitchResult = this->calculatePitch(data);
    fpitchResult = fpitchResult >= 50 ? fpitchResult : 0;
    float fNote = fpitchResult != 0 ? (69 + 12 * log(fpitchResult / 440) / log(2)) - 20 : 0;
    if(fpitchResult <= 0) fNote = 0;
    return fNote;
}
/*
 * 计算每一帧的MFSHS音高
 * param: frameData: 一帧长度的数据样本点
 * */
float mfshs::calculatePitch(const float* frameData){
    float *fftResult = new float[fftLength / 2];
    memset(fftResult, 0, fftLength / 2);
    float *allFFTResultReal = new float[fftLength], *allFFTResultImg = new float[fftLength];
    for(int i = 0; i < this->windowLength; i++){
        allFFTResultReal[i] = frameData[i] * hammingwindow[i];

    }
    complexWithPhase *allFFTResultAbs = fftMain(this->fftLength, allFFTResultReal);
    int maxResultIndex = 0;
    for(int i = 0; i < fftLength * 1 / 8; i++){
        maxResultIndex = allFFTResultAbs[i].absolute > allFFTResultAbs[maxResultIndex].absolute ? i : maxResultIndex;
        fftResult[i] = allFFTResultAbs[i].absolute;
    }
    delete[] allFFTResultReal;
    delete[] allFFTResultImg;
    delete[] allFFTResultAbs;
    return calculateMFSHSPitch(fftResult, maxResultIndex);
}

/*
 * calculateMFSHSPitch
 * 计算MFSHS音高
 * param: float *fftResult, fftResult计算结果
 * param: int maxResultIndex, fftResult中的最大值的索引
 * */
float mfshs::calculateMFSHSPitch(float *fftResult, int maxResultIndex){

    vector<float> p(this->H);
    for(int i = 0; i < this->H; i++){
        p[i] = 0;
        int L = min(10, int(floor((i + 1) * this->fftLength / 8 / (maxResultIndex + 1))));
        for(int j = 0; j < L - 1; j++){
            int index = int(round((j + 1) * (maxResultIndex + 1) / (i + 1.0)));
            if (index != 0){
                float power_val = pow(this->h, j);
                p[i] += fftResult[index - 1] * power_val;
            }
        }
    }
    int maxPIndex = int(max_element(p.begin(), p.end()) - p.begin());
    float f0 = round((maxResultIndex + 1) / (maxPIndex + 1)) / this->fftLength;
    if (f0 > (1100.0 / this->sampleRate)){
        f0 = 0;
    }
    delete[] fftResult;
    return f0 * sampleRate;
}

int main(){
//    char* filename = "./9323.mp3";
//    char* outfilename = "./9323.txt";
    mfshs mf("./output.txt");
    mf.run(0, mf.n_frames - 1);
}
