//
// Created by ywm on 2021/7/26.
//

#ifndef SINGTALENT_MODEL_H
#define SINGTALENT_MODEL_H
#include <complex>
#include <vector>
#include <torch/torch.h>
#include "torch/script.h"
#include <string>
namespace CQT{
    void _m_input(const std::vector<std::shared_ptr<std::complex<double>>> &res_frames, const int &n_frames, std::vector<float> &res);

    class _m_OnsetDetector{
    public:
        _m_OnsetDetector(const int &pad_length=4):pad_length(pad_length){}
        /*  模型加载  */
        void load_model(const std::string &modelpath, const std::string &device_type);
        void load_model(const char *modelpath, const char *device_type);
        /*  模型预测  */
        void predict(const std::vector<std::shared_ptr<std::complex<double>>> &data, const int &n_frames, std::vector<float> &res);
    private:
        /* model of onset detector */
        torch::jit::script::Module module;
        /*  预测时组帧大小 */
        int pad_length;
    };
}

#endif //SINGTALENT_MODEL_H
