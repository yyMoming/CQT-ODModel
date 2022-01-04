//
// Created by yyMoming on 2021/7/8.
//

#ifndef SINGTALENT__M_CQT_H
#define SINGTALENT__M_CQT_H
#include "unistd.h"
#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <memory>
#include <fftw3.h>
#include <numeric>
#include <iomanip>
#include "_m_stft.h"
using namespace std;
namespace CQT{
    template<class dtype>
    class _m_cqt{
    private:
        float sr; // 采样率
        float fmin;   //最小频率
        int n_filters;  //滤波器个数
        int bins_per_octave;    //每倍程频点数
        int n_bins; //总共的频点数
        double Q = 51.43862596945148;    //
        int n_fft;  //每次生成滤波器组要进行fft变换的长度
        vector<std::shared_ptr<std::complex<dtype>>> fft_basis;
        vector<std::shared_ptr<std::complex<dtype>>> filters;
        vector<float> inter_win;    // 采样点采样滤波器
        vector<float> inter_delta;
        int inter_win_len = 8193;
        int precision =512;
        float sample_ratio = 0.5;   //每次都进行下采样
        vector<float> lengths; // float *lengths;
        int n_frames;
        int hop_length;

        ofstream real;
        ofstream imag;
    public:
        vector<std::shared_ptr<std::complex<dtype>>> res_frames;
    public:
        _m_cqt(const float &sample_rate, const float &fmin, const int &bins_per_octave, const int &bins, const int& hop_length);
        ~_m_cqt();
        void calc(dtype* y, const int & y_len, const float& sr=44100);
        vector<std::shared_ptr<std::complex<dtype> >> cqt_res() const;
        int get_n_frames() const;
    private:
        void resample_f(const dtype *signal, const int &signal_len, dtype *resample_res) const;
        void constant_q(const float &my_sr, const float &fmin_t, const int& n_filters);
        void cqt_filter_fft(float &my_sr, float &fmin_t);
        void cqt_response(const dtype* y, const int& y_len, const int & hop_length);
    };
    template<class dtype>
    _m_cqt<dtype>::~_m_cqt() {
        // for(auto x: filters){
        //     if(x)
        //         delete[] x;

        // }
        // for(auto x: fft_basis){
        //     if(x)
        //         delete[] x;

        // }
        // for(auto x: res_frames){
        //     if(x)
        //         delete[] x;
        // }
        real.close();
        imag.close();
        cout << "_m_cqt has been destroyed!" << endl;
    }
    /*
     * cqt 模板类构造函数
     * sample_rate: 输入的音频信号采样率
     * fmin: CQT 变换中最小频点的频率（Hz）
     * bins_per_octave: 每个倍程的频点数（resolution）
     * bins: 变换结果的总频点数
     * */
    template<class dtype>
    _m_cqt<dtype>::_m_cqt(const float &sample_rate, const float &fmin, const int &bins_per_octave,
                          const int &bins, const int &hop_length):sr(sample_rate), fmin(fmin), bins_per_octave(bins_per_octave),
                          n_bins(bins), hop_length(hop_length){
        fstream f("../files/interp_win_param.txt");
        Q = double(1.0) / (pow(2.0, 1.0 / bins_per_octave) - 1);
        inter_win.resize(inter_win_len);
        inter_delta.resize(inter_win_len - 1);
        if(f.is_open()){
            for (int i = 0; i < inter_win_len; i++) {
                f >> inter_win[i];
                if (i > 0) {
                    inter_delta[i - 1] = inter_win[i] - inter_win[i - 1];
                }
            }
        }
        else{
            cout << "文件打开失败" << endl;
            f.close();
            pthread_exit(NULL);
        }

        f.close();
        real.open("cqt_real.txt", ios::out);
        imag.open("cqt_imag.txt", ios::out);
        cout << "template class _m_cqt has initialized" << endl;
    }
    template<class dtype>
    inline vector<std::shared_ptr<std::complex<dtype> >> _m_cqt<dtype>::cqt_res() const{
        return res_frames;
    }

    template<class dtype>
    inline int _m_cqt<dtype>::get_n_frames() const {
        return n_frames;
    }
    /*
     * y: 输入的采样信号
     * y_len: 输入信号的长度
     * sr: cqt 使用的采样率
     * */
    template<class dtype>
    void _m_cqt<dtype>::calc(dtype* y, const int & y_len, const float& sr){
        /*
            删除fft_basis, filters 和 res_frames
        */
       {
           vector<std::shared_ptr<std::complex<dtype>>> v3, v1, v2;
           fft_basis.swap(v1);
           filters.swap(v2);
           res_frames.swap(v3);
       }
        // 总帧数
        this->n_frames = floor(double(y_len) / hop_length);
        vector<dtype> freqs(n_bins);
        for(int i = 0; i < n_bins; i++){
            freqs[i] = fmin * pow(2.0, i / dtype(bins_per_octave));
        }
        // 滤波器个数
        this->n_filters = n_bins > bins_per_octave ? bins_per_octave : n_bins;
        float fmin_t = freqs[n_bins - bins_per_octave];
        float fmax_t = freqs[n_bins - 1];
        int n_octave = floor(log2(n_bins));
        int n_samples = ceil(double(y_len) * sample_ratio);

        unique_ptr<dtype> resample(new dtype[n_samples]);
        dtype *resample_res = resample.get();

        int my_len = n_samples;
        int my_hop_length = hop_length * sample_ratio;
        float my_sr = sr * sample_ratio;
        n_fft = pow(2.0, (double)ceil(log(hop_length)/log(2.0))) > n_fft ? static_cast<int>(pow(2.0, (double)ceil(log(hop_length)/log(2.0)))) :n_fft;
        lengths.resize(n_filters);
        constant_q(my_sr, fmin_t, n_filters);
//
//        //进行 行稀疏 处理之后的结果fft_basis
        cqt_filter_fft(my_sr, fmin_t);
        
        resample_f(y, y_len, resample_res);

        cqt_response(resample_res, my_len, my_hop_length);

        unique_ptr<dtype> copy_samples(new dtype[my_len]);
        dtype *copy_samples_ptr = copy_samples.get();
        std::complex<dtype> mul_num(sqrt(2.0), 0);
        dtype *t;
        for(int i = 0; i < n_octave - 1; i ++){
            for(auto x: fft_basis){
                transform(x.get(), x.get() + n_fft / 2 + 1, x.get(), std::bind(multiplies<complex<dtype>>(), std::placeholders::_1, mul_num));
            }

            if(i % 2 == 0){
                memset(copy_samples_ptr, 0, my_len * sizeof(dtype));
                resample_f(resample_res, my_len, copy_samples_ptr);
                t = copy_samples_ptr;
            }
            else{
                memset(resample_res, 0, my_len * sizeof(dtype));
                resample_f(copy_samples_ptr, my_len, resample_res);
                t = resample_res;
            }

            my_len = ceil(float(my_len) * sample_ratio);
            my_sr /= 2.0;
            my_hop_length = int(my_hop_length / 2);
//            auto tsum = accumulate(t, t + my_len, dtype(0));
            cqt_response(t, my_len, my_hop_length);
        }

        /*
         * 调整顺序
         * */
        for(int i = 0; i < n_octave / 2; i++){
            int high_freq = bins_per_octave * i;
            int low_freq = bins_per_octave * (n_octave - 1 - i);
            for(int j = 0; j < bins_per_octave; j++){
                auto t = res_frames[high_freq + j];
                res_frames[high_freq + j] = res_frames[low_freq +j];
                res_frames[low_freq +j] = t;
            }
        }
        res_frames.assign(res_frames.begin() + bins_per_octave * n_octave - n_bins, res_frames.end());
    }

    /*
     * cqt 采样函数
     * signal: 待采样的信号
     * signal_len: 带采样的信号长度
     * resample_res: 采样结果保存
     * */
    template<class dtype>
    void _m_cqt<dtype>::resample_f(const dtype *signal, const int &signal_len, dtype *resample_res) const{
        float scale = sample_ratio > 1 ? 1 : sample_ratio;
        float time_increment = 1.0f / sample_ratio;
        int index_step = static_cast<int>(scale * precision);
        float time_register = 0.0;

        int n = 0;
        float frac = 0.0;
        float index_frac = 0.0;
        int offset = 0;
        float eta = 0.0;
        float weight = 0.0;
//    vector<float> sig_in(signal, signal+signal_len);
        int n_out = int(sample_ratio * signal_len);
        int n_samples = ceil(sample_ratio * signal_len);
//        cout << "n_out:" << n_out << endl;
//        cout << "n_samples:" << n_samples << endl;
//    vector<float> res_sam(resample_res, resample_res+n_out);
        for(int t = 0; t < n_out; t++){
            // Grab the top bits as an index to the input buffer
            n = static_cast<int>(time_register);
            // Grab the fractional component of the time index
            frac = scale * (time_register - n);
            //Offset into the filter
            index_frac = frac * precision;
            offset = floor(index_frac);
            //Interpolation factor
            eta = index_frac - offset;
            // Compute the left wing of the filter response
            int imax = (n+1) > floor((inter_win_len - offset)/index_step) ? floor((inter_win_len - offset)/index_step): (n+1);
            for (int i = 0; i < imax; i++) {
                weight = (inter_win[offset + i * index_step] + eta * inter_delta[offset + i * index_step]);
                resample_res[t] += weight * signal[n-i];
            }
            //Invert P
            frac = scale - frac;
            //Offset into the filter
            index_frac = frac * precision;
            offset = floor(index_frac);
            //Interpolation factor
            eta = index_frac - offset;
            int k_max = (signal_len - n - 1) > floor((inter_win_len-offset)/index_step) ? floor((inter_win_len-offset)/index_step) : (signal_len - n - 1);
            for (int k = 0; k < k_max; k++) {
                weight = (inter_win[offset + k * index_step] + eta * inter_delta[offset + k * index_step]);
                resample_res[t] += weight * signal[n+k+1];
            }
            time_register += time_increment;
        }
        if (n_samples > n_out) {
            resample_res[n_samples - 1]=0.0;
        }
        float div = 2.0 * sqrt(sample_ratio); //特殊2
        for(int i = 0; i < n_samples; i++){
            resample_res[i] = resample_res[i] / div;
        }

    }

    /*
     * cqt 滤波器
     * my_sr: 当前sample rate（可能经过了采样）
     * fmin_t: 当前倍程的最低频率
     * n_filters: 每倍程滤波器个数
     * */
    template<class dtype>
    void _m_cqt<dtype>::constant_q(const float &my_sr, const float &fmin_t, const int& n_filters){
        /*
         * 临时 lengths
         * */
            lengths.resize(n_filters);
        for(int i=0; i < n_filters; i++){
            auto temp_f = fmin_t * pow(2.0, i / (float)bins_per_octave);
            this->lengths[i] = this->Q * my_sr / temp_f;
        }
        float max_in_length = *max_element(lengths.begin(), lengths.end());
        int max_len = static_cast<int>( pow(2,ceil(log(max_in_length) / log(2))));
        this->n_fft = max_len;    //每次生成滤波器组要进行fft变换的长度
        unique_ptr<std::complex<dtype>> exps(new complex<dtype>[max_len]);
        unique_ptr<std::complex<dtype>> window(new complex<dtype>[max_len]);


        this->filters.resize(n_filters);
        auto exps_ptr = exps.get(), window_ptr = window.get();
        for(int i = 0; i < n_filters; i++){

            int win_size = -1 * floor(-1*lengths[i]/2) + int(lengths[i]/2);
            std::complex<dtype> *complexPointer = new std::complex<dtype>[max_len];
            hanningWindow(window_ptr, win_size);
            for (int j = floor(-1 * lengths[i]/2); j < lengths[i]/2; j++) {
                double inner = 2 * M_PI * Q * (double)j / lengths[i];
                exps_ptr[j + win_size / 2 + 1] = std::complex<dtype>(cos(inner), sin(inner) );
////                exps_ptr[j + win_size / 2 + 1].imag = sin(inner);
//                exps_ptr[j + win_size / 2 + 1] = std::exp(inner);

            }
//            auto tmp = complexPointer;
            transform(exps_ptr, exps_ptr + win_size, window_ptr, complexPointer, multiplies<complex<dtype>>());
            vector<dtype> mags(win_size);
            for(int j = 0; j < win_size; j++){
                mags[j] = abs(complexPointer[j]);
            }
            dtype mags_sum = accumulate(mags.begin(), mags.end(), 0.0);
            for(int j = 0; j < win_size; j++){
                complexPointer[j] /= mags_sum;
            }
            uint pad_left = floor((max_len - win_size) / 2);
            memmove(complexPointer + pad_left, complexPointer, win_size * sizeof(std::complex<dtype>));
//            copy_backward(complexPointer, complexPointer + win_size, complexPointer + pad_left);
            fill_n(complexPointer, pad_left, std::complex<dtype>());
            const float nFloat = lengths[i];
//            const float n_max_len = nFloat / float(max_len);
            const complex<dtype> n_max_len(nFloat / float(max_len), 0);
            for(int j = 0; j < max_len; j++){
                complexPointer[j] = complexPointer[j] * n_max_len;
            }
            this->filters[i] = std::shared_ptr<std::complex<dtype>>(complexPointer);
            // this->filters[i] = std::shared_ptr<std::complex<dtype>>(complexPointer, [](std::complex<dtype> *ptr)
            //                                                         { delete[] ptr; });
        }

    }

    /*
     * cqt 滤波器变稀疏矩阵
     * my_sr: 当前sample rate（可能经过了采样）
     * fmin_t: 当前倍程的最低频率
     * */
    template<class dtype>
    void _m_cqt<dtype>::cqt_filter_fft(float &my_sr, float &fmin_t){
        unique_ptr<std::complex<dtype>> fft_res(new std::complex<dtype>[n_fft / 2 + 1]);
        vector<dtype> cumulative_mag(n_fft / 2 + 1);
        vector<dtype> mags(n_fft / 2 + 1), mags_sort(n_fft / 2 + 1);
        fftw_complex *in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft), *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft);
        fftw_plan p;
        auto fft_res_ptr = fft_res.get();
        int count = 0;
        for(auto x: filters){

            fftw_complex *tmp = reinterpret_cast<fftw_complex *>(x.get());
            memcpy(in, tmp, n_fft * sizeof(fftw_complex));
            p = fftw_plan_dft_1d(n_fft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p);
            fft_res_ptr = reinterpret_cast<std::complex<dtype> *>(out);

            std::complex<dtype> *sparse_vec = new std::complex<dtype>[n_fft / 2 + 1];
            // auto f = std::bind(std::abs<dtype>, std::complex<dtype>(1.0, 3.0));
            transform(fft_res_ptr, fft_res_ptr + n_fft / 2 + 1, mags.begin(), std::bind(std::abs<dtype>, std::placeholders::_1));
            dtype norm = accumulate(mags.begin(), mags.end(), 0.0);

            copy(mags.begin(), mags.end(), mags_sort.begin());
            sort(mags_sort.begin(), mags_sort.end());
            cumulative_mag[0] = mags_sort[0] / norm;
            for(int j = 1; j < n_fft / 2 + 1; j++){
                cumulative_mag[j] = mags_sort[j] / norm + cumulative_mag[j - 1];
            }
            int threshold_idx;
            for (int k = 0; k < n_fft / 2 + 1; k++){
                if(cumulative_mag[k] >= 0.01){
                    threshold_idx = k;
//                    k_vec.push_back(k);
                    break;
                }
            }
            for (int k = 0; k < n_fft / 2 + 1; k++){
                if(mags[k] >= mags_sort[threshold_idx]){
                    sparse_vec[k] = fft_res_ptr[k];
                }
            }
            fft_basis.push_back(std::shared_ptr<std::complex<dtype>>(sparse_vec));
            // fft_basis.push_back(std::shared_ptr<std::complex<dtype>>(sparse_vec, [](std::complex<dtype> *ptr)
            //                                                          { delete[] ptr; }));
            fftw_destroy_plan(p);
        }
        fftw_free(in);
        fftw_free(out);
    }

    /*
     * y: 当前采样信号
     * y_len: y的长度
     * hop_length: stft变换的hop size
     * */
    template<class dtype>
    void _m_cqt<dtype>::cqt_response(const dtype* y, const int &y_len, const int &hop_length){
//        auto sum = accumulate(y, y + y_len, dtype(0));
        // vector<std::complex<dtype> *> stft_ret;
        vector<std::shared_ptr<std::complex<dtype>>> stft_ret;
        _m_stft(y, y_len, n_fft, hop_length, stft_ret);
        /*
        ofstream _imag("cqt_stft_imag.csv", ios::out);
        ofstream _real("cqt_stft_real.csv", ios::out);
        for(auto x: stft_ret){
            for(int i = 0; i < n_fft / 2 + 1; i++){
                auto t = *(x + i);
                _real << t.real() << ",";
                _imag << t.imag() << ",";
            }
            _real << endl;
            _imag << endl;
        }
        exit(1);*/
        int n_frames = floor(double(y_len) / hop_length);
        cout <<"[" << y_len << ":" << hop_length << "]" << n_frames << endl;
        std::complex<dtype> *c = new std::complex<dtype>[n_fft / 2 + 1];
        for(auto x: fft_basis){
            std::complex<dtype> *mul_res = new complex<dtype>[stft_ret.size()];
            auto t = mul_res;
            for(auto frame: stft_ret){

                transform(x.get(), x.get() + n_fft / 2 + 1, frame.get(), c, multiplies<std::complex<dtype>>());
                auto sum_c = accumulate(c, c + n_fft / 2 + 1, std::complex<dtype>(0, 0));
                *t = sum_c;
//                real << sum_c.real() << ",";
//                imag << sum_c.imag() << ",";
                t++;
            }
//            real << endl;
//            imag << endl;
            auto mul_res_sum = accumulate(mul_res, mul_res + n_frames, std::complex<dtype>(0, 0));
            // res_frames.push_back(std::shared_ptr<std::complex<dtype>>(mul_res, [](std::complex<dtype> *ptr)
            //                                                           { delete[] ptr; }));
            res_frames.push_back(std::shared_ptr<std::complex<dtype>>(mul_res));
        }
        // 释放堆内存
        delete[] c;
    }
};
#endif //SINGTALENT__M_CQT_H
