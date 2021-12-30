//
// Created by ywm on 2021/7/26.
//

#ifndef SINGTALENT_MODEL_H
#define SINGTALENT_MODEL_H
#include <complex>
#include <vector>
#include <torch/torch.h>
#include "torch/script.h"
namespace CQT{
    void _m_input(const std::vector<std::complex<double> *> &res_frames, const int &n_frames, std::vector<float> &res);

    class _m_OnsetDetector{
    private:
        /* model of onset detector */
        torch::jit::script::Module module;

    };
}

#endif //SINGTALENT_MODEL_H
