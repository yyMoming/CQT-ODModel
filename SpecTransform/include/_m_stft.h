//
// Created by yyMoming on 2021/7/7.
//

#ifndef SINGTALENT__M_STFT_H
#define SINGTALENT__M_STFT_H
#include <fftw3.h>
#include <complex>
#include <vector>
#include <fstream>
namespace CQT{
    template<class dtype>
    void _m_stft(const dtype* y, const int& y_len, const int & n_fft, const int & hop_length, std::vector<std::complex<dtype> *> &stft_ret){
        dtype *y_padded = new dtype[y_len + n_fft];
        copy(y, y + y_len, y_padded + n_fft / 2);
        int n_frames = floor(double(y_len + n_fft) / hop_length);
//        vector<std::complex<dtype> *> stft_res;
        fftw_complex *in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft), *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft);
        fftw_plan p;

        for(int i = 0; i < n_frames; i++){
            std::complex<dtype> *res = new std::complex<dtype>[n_fft / 2 + 1];
            for(int k = 0; k < n_fft; k++){
                in[k][0] = y_padded[i * hop_length + k];
                in[k][1] = 0;
            }
            p = fftw_plan_dft_1d(n_fft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p);
//            res = reinterpret_cast<std::complex<double> *>(out);
            memcpy(res, out, (n_fft / 2 + 1) * sizeof(std::complex<dtype>));
            stft_ret.push_back(res);
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        delete[] y_padded;
    }
    /*
     * 计算并输出STFT结果
     * */
    template<class dtype>
    void _m_stft(const dtype* y, const int& y_len, const int & n_fft, const int & hop_length){
        std::fstream fimag("_stft_imag.txt", std::ios::out);
        std::fstream freal("_stft_real.txt", std::ios::out);
        dtype *y_padded = new dtype[y_len + n_fft];
        copy(y, y + y_len, y_padded + n_fft / 2);
        int n_frames = floor(double(y_len + n_fft) / hop_length);
        //        vector<std::complex<dtype> *> stft_res;
        fftw_complex *in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft), *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft);
        fftw_plan p;
        for(int i = 0; i < n_frames; i++){
        //        copy(y_padded + i * n_fft, y_padded + (i + 1) * n_fft, in);
        //            memcpy(in, y_padded + i * hop_length, n_fft * sizeof(fftw_complex));
            for(int k = 0; k < n_fft; k++){
                in[k][0] = y_padded[i * hop_length + k];
                in[k][1] = 0;
            }
            p = fftw_plan_dft_1d(n_fft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p);
            for(int k = 0; k < n_fft; k++){
                auto t = *(out + k);
                if(k == n_fft - 1){
                    fimag << t[1] << std::endl;
                    freal << t[0] << std::endl;
                    break;
                }
                fimag << t[1] << ",";
                freal << t[0] << ",";

            }
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        delete[] y_padded;
    }


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
}
#endif //SINGTALENT__M_STFT_H
