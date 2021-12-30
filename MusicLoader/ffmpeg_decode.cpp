//
// Created by yyMoming on 2021/6/28.
//

//***************************************************************
// @file:    test.c
// @author:  dingfang
// @date    2019-07-24 18:55:16
//***************************************************************

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#ifdef __cplusplus
};
#endif

int openCodecContext(const AVFormatContext *pFormatCtx, int *pStreamIndex, enum AVMediaType type, AVCodecContext **ppCodecCtx)
{
    int streamIdx = -1;
    // 获取流下标
    for (int i = 0; i < pFormatCtx->nb_streams; i++)
    {
        if (pFormatCtx->streams[i]->codec->codec_type == type)
        {
            streamIdx = i;
            break;
        }
    }
    if (streamIdx == -1)
    {
        printf("find video stream failed!\n");
        exit(-2);
    }
    // 寻找解码器
    AVCodecContext  *pCodecCtx = pFormatCtx->streams[streamIdx]->codec;
    AVCodec *pCodec = avcodec_find_decoder(pCodecCtx->codec_id);
    if (NULL == pCodec)
    {
        printf("avcode find decoder failed!\n");
        exit(-2);
    }

    //打开解码器
    if (avcodec_open2(pCodecCtx, pCodec, NULL) < 0)
    {
        printf("avcode open failed!\n");
        exit(-2);
    }
    *ppCodecCtx        = pCodecCtx;
    *pStreamIndex    = streamIdx;

    return 0;
}

int main(void)
{
    AVFormatContext    *pInFormatCtx    = NULL;
    AVCodecContext  *pVideoCodecCtx    = NULL;
    AVCodecContext  *pAudioCodecCtx    = NULL;
    AVPacket *pPacket    = NULL;
    AVFrame *pFrame        = NULL;
    int ret;
    /* 支持本地文件和网络url */
    const char streamUrl[] = "/home/ywm/cpp/SingTalent/9323.mp3";

    /* 1. 注册 */
    av_register_all();

    pInFormatCtx = avformat_alloc_context();

    /* 2. 打开流 */
    if(avformat_open_input(&pInFormatCtx, streamUrl, NULL, NULL) != 0)
    {
        printf("Couldn't open input stream.\n");
        return -1;
    }

    /* 3. 获取流的信息 */
    if(avformat_find_stream_info(pInFormatCtx, NULL) < 0)
    {
        printf("Couldn't find stream information.\n");
        return -1;
    }

    int videoStreamIdx = -1;
    int audioStreamIdx = -1;
    /* 4. 寻找并打开解码器 */
    openCodecContext(pInFormatCtx, &videoStreamIdx, AVMEDIA_TYPE_VIDEO, &pVideoCodecCtx);
    openCodecContext(pInFormatCtx, &audioStreamIdx, AVMEDIA_TYPE_AUDIO, &pAudioCodecCtx);

    pPacket    = av_packet_alloc();
    pFrame    = av_frame_alloc();

    int cnt = 30;
    while (cnt--)
    {
        /* 5. 读流数据, 未解码的数据存放于pPacket */
        ret = av_read_frame(pInFormatCtx, pPacket);
        if (ret < 0)
        {
            printf("av_read_frame error\n");
            break;
        }

        /* 6. 解码, 解码后的数据存放于pFrame */
        /* 视频解码 */
        if (pPacket->stream_index == videoStreamIdx)
        {
            avcodec_decode_video2(pVideoCodecCtx, pFrame, &ret, pPacket);
            if (ret == 0)
            {
                printf("video decodec error!\n");
                continue;
            }
            printf("* * * * * * video * * * * * * * * *\n");
            printf("___height: [%d]\n", pFrame->height);
            printf("____width: [%d]\n", pFrame->width);
            printf("pict_type: [%d]\n", pFrame->pict_type);
            printf("___format: [%d]\n", pFrame->format);
            printf("* * * * * * * * * * * * * * * * * * *\n\n");
        }

        /* 音频解码 */
        if (pPacket->stream_index == audioStreamIdx)
        {
            avcodec_decode_audio4(pAudioCodecCtx, pFrame, &ret, pPacket);
            if (ret < 0)
            {
                printf("audio decodec error!\n");
                continue;
            }
            printf("* * * * * * audio * * * * * * * * * *\n");
            printf("____nb_samples: [%d]\n", pFrame->nb_samples);
            printf("__samples_rate: [%d]\n", pFrame->sample_rate);
            printf("channel_layout: [%lu]\n", pFrame->channel_layout);
            printf("________format: [%d]\n", pFrame->format);
            printf("* * * * * * * * * * * * * * * * * * *\n\n");
        }
        av_packet_unref(pPacket);
    }

    av_frame_free(&pFrame);
    av_packet_free(&pPacket);
    avcodec_close(pVideoCodecCtx);
    avcodec_close(pAudioCodecCtx);
    avformat_close_input(&pInFormatCtx);

    return 0;
}