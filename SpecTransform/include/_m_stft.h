//
// Created by yyMoming on 2021/7/7.
//

#ifndef SINGTALENT__M_STFT_H
#define SINGTALENT__M_STFT_H
#include <fftw3.h>
#include <complex>
#include <vector>
#include <fstream>
#include <functional>
#include <algorithm>
#include <memory>
namespace CQT{

    template<class dtype>
    inline void hanningWindow(std::vector<std::complex<dtype>> &win, int win_size){

        for(int i = 0; i < win_size; i++){
            win[i] = 0.5 - 0.5 * cos(2 * M_PI * i / (win_size - 1));

        }
    }

    template<class dtype>
    inline void hanningWindow(dtype *win, int win_size){

        for(int i = 0; i < win_size; i++){
            win[i] = 0.5 - 0.5 * cos(2 * M_PI * i / (win_size - 1));
        }
    }

    template<class dtype>
    inline void hanningWindow(std::complex<dtype> *win, int win_size){

        for(int i = 0; i < win_size; i++){
            dtype tmp = 0.5 - 0.5 * cos(2 * M_PI * i / (win_size - 1));
            win[i] = tmp;
//            cout << M_PI;
//            win[i].imag = 0;
        }
    }

    template<class dtype>
    void _m_stft(const dtype* y, const int& y_len, const int & n_fft, const int & hop_length, std::vector<std::shared_ptr<std::complex<dtype>>> &stft_ret){
        dtype *y_padded = new dtype[y_len + n_fft];
        dtype *pad = new dtype[n_fft / 2];(y + 1, y + 1 + n_fft / 2);
        copy(y + 1, y + 1 + n_fft / 2, pad);
        reverse(pad, pad + n_fft / 2);
        copy(pad, pad + n_fft / 2, y_padded);

        copy(y + y_len - n_fft / 2, y + y_len, pad);
        reverse(pad, pad + n_fft / 2);
        copy(pad, pad + n_fft / 2, y_padded + y_len + n_fft / 2);
        
        copy(y, y + y_len, y_padded + n_fft / 2);
        int n_frames = floor(double(y_len) / hop_length);

        fftw_complex *in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft), *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft);
        fftw_plan p;

        std::complex<dtype> *ptr = nullptr;
        for(int i = 0; i < n_frames; i++){
            std::complex<dtype> *res = new std::complex<dtype>[n_fft / 2 + 1];
            memset(res, 0, sizeof(std::complex<dtype>) * (n_fft / 2 + 1));
            for(int k = 0; k < n_fft; k++){
                in[k][0] = y_padded[i * hop_length + k];
                in[k][1] = 0;
            }
            ptr = reinterpret_cast<std::complex<dtype> *>(in);
            p = fftw_plan_dft_1d(n_fft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p);
            memcpy(res, out, (n_fft / 2 + 1) * sizeof(std::complex<dtype>));
            stft_ret.push_back(std::shared_ptr<std::complex<dtype>>(res));
            fftw_destroy_plan(p);
        }
        
        fftw_free(in);
        fftw_free(out);
        delete[] y_padded;
        delete[] pad;
    }

    /*
     * 计算并输出STFT结果, 加窗 hanning
     * */
    template<class dtype>
    void _m_stft(const dtype* y, const int& y_len, const int & n_fft, const int & hop_length){
        std::fstream fimag("_stft_imag.txt", std::ios::out);
        std::fstream freal("_stft_real.txt", std::ios::out);
        dtype *y_padded = new dtype[y_len + n_fft];
        copy(y, y + y_len, y_padded + n_fft / 2);
        std::complex<dtype> *window = new std::complex<dtype>[n_fft];
        hanningWindow(window, n_fft);
        int n_frames = floor(double(y_len) / hop_length);

        //        vector<std::complex<dtype> *> stft_res;
        fftw_complex *in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft), *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft);
        fftw_plan p;
        std::complex<dtype> *ptr = nullptr;
        for(int i = 0; i < n_frames; i++){
            for(int k = 0; k < n_fft; k++){
                in[k][0] = y_padded[i * hop_length + k];
                in[k][1] = 0;
            }
            ptr = reinterpret_cast<std::complex<dtype> *>(in);
            std::transform(ptr, ptr + n_fft, window, ptr, std::multiplies<std::complex<dtype>>());
            p = fftw_plan_dft_1d(n_fft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p);
            for(int k = 0; k < n_fft; k++){
                auto t = *(out + k);
                if(k == n_fft - 1){
                    fimag << t[1] << std::endl;
                    freal << t[0] << std::endl;
                    break;
                }else{
                    fimag << t[1] << ",";
                    freal << t[0] << ",";
                }


            }
//            fimag << std::endl;
//            freal << std::endl;
            fftw_destroy_plan(p);
        }
        fimag.close();
        freal.close();
        
        fftw_free(in);
        fftw_free(out);
        delete[] y_padded;
        delete[] window;
    }
}
#endif //SINGTALENT__M_STFT_H
