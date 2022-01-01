//
// Created by ywm on 2021/4/26.
//

#ifndef SIGHTSING_MUSIC_LOADER_H
#define SIGHTSING_MUSIC_LOADER_H
#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __cplusplus
extern "C"
{
#endif
#include <libavutil/frame.h>
#include <libavutil/mem.h>
#include <libavcodec/avcodec.h>
#ifdef __cplusplus
};
#endif


#define INBUF_SIZE 4096

#define AUDIO_INBUF_SIZE 20480
#define AUDIO_REFILL_THRESH 4096
using namespace std;
// decode for mp3 file
static void decode(AVCodecContext *dec_ctx, AVPacket *pkt, AVFrame *frame,
                   FILE *outfile);
void audio_decode_example(const char *outfilename, const char *filename);

void mp3_loader(char *filename, char *outfilename);

void mp3_loader(char *filename, vector<float> &pcm_ptr);

void mp3_loader(char *filename, vector<double> &pcm_ptr);

void pcm_loader(std::fstream &fs, std::vector<float> &res);

#endif //SIGHTSING_MUSIC_LOADER_H
