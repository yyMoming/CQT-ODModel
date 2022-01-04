//
// Created by ywm on 2021/7/20.
//
#include <torch/torch.h>
#include "torch/script.h"
#include <iostream>
#include <vector>
#include <time.h>
#include "Model.h"
#include <functional>
#include <algorithm>
#include <sys/stat.h>
namespace CQT
{

    void _m_OnsetDetector::load_model(const std::string &modelpath, const std::string &device_type)
    {
        time_t start = clock();
        try
        {
            /* code */
            torch::DeviceType device;
            if (device_type == "cpu")
            {
                device = at::kCPU;
            }
            else
            {
                device = at::kCUDA;
            }
            this->module = torch::jit::load("/home/ywm/MUSIC/Codes/onset_train/model.pt", device);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            exit(1);
        }
        time_t endd = clock();

        std::cout << "Load model successful! cost: " << double(endd - start) / CLOCKS_PER_SEC << std::endl;
    }
    void _m_OnsetDetector::load_model(const char *modelpath, const char *device_type)
    {
        time_t start = clock();
        try
        {
            /* code */
            torch::DeviceType device;
            if (strcmp(device_type, "cpu") == 0)
            {
                device = at::kCPU;
            }
            else
            {
                device = at::kCUDA;
            }
            this->module = torch::jit::load("/home/ywm/MUSIC/Codes/onset_train/model.pt", device);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            exit(1);
        }
        time_t endd = clock();
        std::cout << "Load model successful! cost: " << double(endd - start) / CLOCKS_PER_SEC << std::endl;
    }

    void _m_OnsetDetector::predict(const std::vector<std::shared_ptr<std::complex<double>>> &data, const int &n_frames, std::vector<float> &res){
        this->module.eval();
        int n_bins = data.size();
        std::cout << "=================预测数据开始加载===================" << std::endl;
        time_t start = clock();
        std::cout << "start:" <<  start << std::endl;
        torch::Tensor T = torch::ones({1, n_frames});
        for (auto x_ptr : data)
        {
            auto x = x_ptr.get();
            std::vector<double> dvec(n_frames);
            for (int i = 0; i < n_frames; i++, x++)
            {
                dvec[i] = abs(*(x));
            }
            torch::Tensor tensor = torch::tensor(dvec);
            tensor = tensor.unsqueeze(0);
            T = torch::cat({T, tensor}, 0);
        }
        T = T.slice(0, 1, n_bins + 1);
        std::cout << "T.size():" << T.sizes() << std::endl;
        torch::Tensor padding_tensor = torch::zeros({n_bins, this->pad_length});
        T = torch::cat({padding_tensor, T, padding_tensor}, 1);
        
        std::cout << "=================预测数据加载完毕===================" << std::endl;

        std::cout << "=================模型预测开始===================" << std::endl;
       
        for (int i = 0; i < n_frames; i++)
        {
            torch::Tensor tensor = T.slice(1, i, i + 9);
            auto min_val = torch::min(tensor);
            auto max_val = torch::max(tensor);
            tensor = (tensor - min_val) / (max_val - min_val + 1e-9);
            tensor = tensor.unsqueeze(0).unsqueeze(0);
            at::Tensor output = module.forward({tensor}).toTensor();
            // frames.push_back(output);
            res.push_back(output[0][0].item().toFloat());
        }
        time_t endd = clock();
        std::cout << "endd:" <<  endd << std::endl;
        std::cout << "=================模型预测完毕===================" << std::endl;
        std::cout << "================= 耗时： " << double(endd - start) / CLOCKS_PER_SEC << "===================" << std::endl;

}

    int padding_length = 4;
    void TorchTest(torch::DeviceType device_type)
    {
        torch::jit::script::Module module = torch::jit::load("/home/ywm/MUSIC/Codes/onset_train/model.pt", device_type);
        //    assert(module != nullptr);
        std::cout << "Load model successful!" << std::endl;
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(torch::ones({1, 1, 267, 9}).cuda());
        at::Tensor output = module.forward(inputs).toTensor();
        //    auto max_result = output.max(1, true);
        //    auto max_index = std::get<1>(max_result).item<float>();
        std::cout << output << std::endl;
    }

    void _m_input(const std::vector<std::shared_ptr<std::complex<double>>> &res_frames, const int &n_frames, std::vector<float> &res)
    {
        //        auto f = std::bind(std::abs<double>, std::placeholders::_1);
        //        auto f = std::bind(std::multiplies<int>(), std::placeholders::_1, std::placeholders::_2);
        //        std::cout << f(3, 4) << std::endl;
        time_t start = clock();
        torch::jit::script::Module module = torch::jit::load("/home/ywm/MUSIC/Codes/onset_train/model.pt", at::kCPU);
        time_t endd = clock();
        std::cout << "读取模型耗时：" << double(endd - start) / CLOCKS_PER_SEC << std::endl;
        //    assert(module != nullptr);
        module.eval();
        std::cout << "Load model successful!" << std::endl;

        torch::Tensor T = torch::ones({1, n_frames});
        int n_bins = 0;
        std::cout << "=================预测数据开始加载===================" << std::endl;
        start = clock();
        for (auto x_ptr : res_frames)
        {
            auto x = x_ptr.get();
            n_bins++;
            std::vector<double> dvec(n_frames);
            for (int i = 0; i < n_frames; i++, x++)
            {
                dvec[i] = abs(*(x));
            }
            torch::Tensor tensor = torch::tensor(dvec);
            tensor = tensor.unsqueeze(0);
            T = torch::cat({T, tensor}, 0);
        }
        T = T.slice(0, 1, n_bins + 1);
        std::cout << "T.size():" << T.sizes() << std::endl;
        torch::Tensor padding_tensor = torch::zeros({n_bins, padding_length});
        T = torch::cat({padding_tensor, T, padding_tensor}, 1);
        std::cout << "padded T.size()" << T.sizes() << std::endl;
        //        endd = clock();
        //        std::cout << "================= 预测数据加载完成, 耗时:" << double(endd - start) / CLOCKS_PER_SEC << "===================" << std::endl;
        std::vector<torch::Tensor> frames;
        std::vector<torch::jit::IValue> inputs;

        //        std::cout << "=================模型开始预测===================" << std::endl;
        //        start = clock();
        
        for (int i = 0; i < n_frames; i++)
        {
            torch::Tensor tensor = T.slice(1, i, i + 9);
            auto min_val = torch::min(tensor);
            auto max_val = torch::max(tensor);
            tensor = (tensor - min_val) / (max_val - min_val + 1e-9);
            tensor = tensor.unsqueeze(0).unsqueeze(0);
            at::Tensor output = module.forward({tensor}).toTensor();
            frames.push_back(output);
            res.push_back(output[0][0].item().toFloat());
        }
        endd = clock();
        std::cout << "================= 模型预测完成, 耗时:" << double(endd - start) / CLOCKS_PER_SEC << "===================" << std::endl;
        std::cout << frames.size() << std::endl;
    }
    /*
    int main() {
    //    torch::DeviceType device_type = at::kCPU;
    //    if (torch::cuda::is_available()) {
    //        std::cout << "cuda!" << std::endl;
    //        std::cout << "Device count : " << torch::cuda::device_count() << std::endl;
    //        torch::DeviceType device_type = at::kCUDA;
    //        TorchTest(device_type);
    //    }


        std::vector<double> dvec(1000);
        for(int i = 0; i < 10; i++){
            dvec[i] = i;
        }
    //    for(auto x: dvec){
    //        std::cout << x << ",";
    //    }
        std::cout << std::endl;
    //    torch::Tensor tensor = torch::from_blob(dvec, {10});
        std::time_t start = clock();
        torch::Tensor tensor = torch::tensor(dvec);
        std::time_t endd = clock();
        std::cout << "耗时：" << double(endd - start) / CLOCKS_PER_SEC << std::endl;
    //    std::cout << tensor << std::endl;
    }
    */

}
