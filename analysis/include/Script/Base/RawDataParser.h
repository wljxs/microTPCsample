/**
 * @file RawDataParser.h
 * @brief 原始数据解析器 - 完全封装ROOT文件和TTree操作
 * @author Huang Qixuan
 */

#pragma once

#include "DataModel.h"
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <tuple>
#include <vector>

class RawDataParser {
   public:
    RawDataParser(const std::string& rawFile);
    ~RawDataParser();

    bool Initialize();

    std::unordered_map<int, std::vector<RawData>> LoadEvent(int eventID);

    Long64_t GetTotalEvents() const { return m_numOfEvents; };

    std::vector<float> GetPedstd(int detectorID,int type) const;

    bool LoadPedestals();

   private:
    std::string m_rawFile;

    Long64_t m_numOfEvents{0};

    TFile* m_file{nullptr};
    TTree* m_tree{nullptr};
    TTree* m_pedTree{nullptr};// 用于pedestal数据

    std::vector<unsigned int>* m_apv_id{nullptr};
    std::vector<unsigned int>* m_apv_ch{nullptr};
    std::vector<unsigned int>* m_mm_strip{nullptr};
    std::vector<std::vector<short>>* m_apv_q{nullptr};
    unsigned int m_apv_evt{0};
    std::map<int, std::map<int, std::vector<float>>> m_pedestalMap; // [detID][type] -> stddev列表

    // 硬件通道映射
    std::tuple<int, int, int> MapBoardChannel(
        unsigned int boardID,
        unsigned int channelID,
        unsigned int mm_strip) const;
};
