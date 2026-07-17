#include "Script/BadStripScript.h"
#include "Algorithm/AnalysisUtils.h"
#include "Detector/DetectorFactory.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"
#include "TH1.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iostream>

void BadStripScript::LoadConfig(const json &config) {
  m_maxEvents = config.value("maxEvents", 0); // 0 = 全部事例
  m_progressInterval = config.value("progressInterval", 1000);
  m_detnum = config.value("detnumber", 4); // 默认分析1700
  m_stripnum = config.value("stripnum", 256);
  m_typenum = config.value("typenum", 2);
  m_badStripThreshold =
      config.value("badStripThreshold", 0.1); // 坏道定义：hit占比超过10%
  // accept both 'threshold' and legacy 'Threshold'
  m_Threshold =  config.value("threshold", config.value("Threshold", 300.0)); // ADC
}

void BadStripScript::Print() const {
  std::cout << "  BadStripScript Configuration:" << std::endl;
  std::cout << "  Max Events: "
            << (m_maxEvents == 0 ? "All" : std::to_string(m_maxEvents))
            << std::endl;
  std::cout << "  Progress Interval: " << m_progressInterval << std::endl;
  std::cout << "  Threshold: " << m_Threshold << " ADC" << std::endl;
  std::cout << "  Detector Number: " << m_detnum << std::endl;
  std::cout << "  Strip Number: " << m_stripnum << std::endl;
  std::cout << "  Type Number: " << m_typenum << std::endl;
  std::cout << "  Bad Strip Threshold: " << m_badStripThreshold << std::endl;
}

bool BadStripScript::Execute() {
  std::cout << "\n========================================" << std::endl;
  std::cout << " Channel Noise Analysis Script" << std::endl;
  std::cout << "========================================" << std::endl;

  auto parser = GetParser();
  if (!parser) {
    std::cerr << "Error: Parser not set!" << std::endl;
    return false;
  }

  // 分析通道噪声
  AnalyzeChannelNoise();

  // 保存结果
  if (!SaveResults()) {
      std::cerr << "Error: Failed to save results!" << std::endl;
      return false;
  }

  std::cout << "\n========================================" << std::endl;
  std::cout << " Channel Noise Analysis Completed" << std::endl;
  std::cout << "========================================" << std::endl;

  return true;
}

void BadStripScript::AnalyzeChannelNoise() {
  auto parser = GetParser();
  Long64_t totalEvents = parser->GetTotalEvents();
  Long64_t eventsToProcess = (m_maxEvents > 0 && m_maxEvents < totalEvents)
                                 ? m_maxEvents
                                 : totalEvents;

  std::cout << "\nAnalyzing " << eventsToProcess << " events..." << std::endl;

  // 临时存储：[detID][type][stripID] -> baseline_list
  std::map<int, std::map<int, std::map<int, std::vector<double>>>> tempData;
  std::map<int, std::map<int, std::map<int, int>>>
      hitOverThresholdCountMap; // [detID][type][stripID] -> hit count
  for (int detID = 1; detID < m_detnum + 1; ++detID) {
    for (int type = 0; type < m_typenum; ++type) {
      for (int stripID = 0; stripID < m_stripnum; ++stripID) {
        hitOverThresholdCountMap[detID][type][stripID] = 0;
      }
    }
  }
  // 遍历所有事例
  for (Long64_t iEvt = 0; iEvt < eventsToProcess; ++iEvt) {
    if ((iEvt % m_progressInterval) == 0 || iEvt == eventsToProcess - 1) {
      std::cout << "\r[Channel Noise] Processed " << (iEvt + 1) << "/"
                << eventsToProcess << std::flush;
    }

    auto rawHitsMap = parser->LoadEvent(iEvt);
    if (rawHitsMap.empty())
      continue;

    // 遍历所有探测器
    for (const auto &[detID, rawHits] : rawHitsMap) {
      // 遍历该探测器的所有hit
      for (const auto &raw : rawHits) {
        int stripID = raw.stripID;
        int type = raw.type;
        const auto &waveform = raw.adc;
        int peak = *std::max_element(waveform.begin(), waveform.end());
        double baseline = static_cast<double>(waveform[0]);
        if ((peak - baseline) > m_Threshold) {
          hitOverThresholdCountMap[detID][type][stripID]++;
        }
        if (waveform.empty())
          continue;
        // baseline只取第一个采样点
        tempData[detID][type][stripID].push_back(baseline);
      }
    }
  }
  TFile *file =
      new TFile((GetOutputDir() + "BadStripAnalysis.root").c_str(), "RECREATE");
  json badStripsJson;
  badStripsJson["version"] = 1;
  badStripsJson["badStrips"] = json::object();
  for (int detID = 1; detID < 5; ++detID) { // 只是分析1700的
    std::string detName = Form("Det%d", detID);
    file->mkdir(detName.c_str());
    file->cd(detName.c_str());
    badStripsJson["badStrips"][std::to_string(detID)] = json::object();
    for (int type = 0; type < 2; ++type) {
      std::string typeName = (type == 0) ? "X" : "Y";
      file->mkdir((detName + "/" + typeName).c_str());
      file->cd((detName + "/" + typeName).c_str());
      float meanHitOverThresholdCount = 0;
      for (const auto &[stripID, hitOverThresholdCount] :
           hitOverThresholdCountMap[detID][type]) {
        meanHitOverThresholdCount += hitOverThresholdCount;
      }
      meanHitOverThresholdCount /= m_stripnum; // strip数量

      TH1D *hitmap =
          new TH1D(Form("hitmap_det%d_type%s", detID, typeName.c_str()),
                   Form("Hit Over Threshold Count for Det %d Type %s", detID,
                        typeName.c_str()),
                   m_stripnum, 0.5, m_stripnum + 0.5);
      TH1D *badStripDist =
          new TH1D(Form("badStripDist_det%d_type%s", detID, typeName.c_str()),
                   Form("Bad Strip Distribution for Det %d Type %s", detID,
                        typeName.c_str()),
                   m_stripnum, 0.5, m_stripnum + 0.5);
      std::vector<int> badList;
      for (int stripID = 1; stripID < m_stripnum + 1; ++stripID) {
        int hitCount = hitOverThresholdCountMap[detID][type][stripID];
        hitmap->SetBinContent(stripID - 1, hitCount);
        if (hitCount < meanHitOverThresholdCount * m_badStripThreshold) {
          badStripDist->SetBinContent(stripID - 1, 1);
          badList.push_back(stripID);
        }
      }
      badStripsJson["badStrips"][std::to_string(detID)][std::to_string(type)] =
          badList;
      hitmap->Write();
      badStripDist->Write();
      delete hitmap;
      delete badStripDist;
      file->cd("..");
    }
    file->cd("..");
  }
  file->Close();

  // Also export machine-readable JSON for later reconstruction stages
  const std::string jsonPath = GetOutputDir() + "bad_strips.json";
  try {
    std::ofstream out(jsonPath);
    out << badStripsJson.dump(2);
    std::cout << "\n[BadStripScript] Exported bad strips JSON: " << jsonPath
              << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "\n[BadStripScript] Failed to write bad strips JSON: "
              << e.what() << std::endl;
  }

  // // 计算每个通道的统计量
  // m_channelNoiseMap.clear();
  // int totalChannels = 0;
  // for(int detID = 1; detID < m_detnum + 1; ++detID) {// 只是分析1700的
  //     for(int type = 0; type < m_typenum; ++type) {
  //         int meanHItOverThresholdCount = 0;
  //         for(const auto& [stripID, hitOverThresholdCount] :
  //         hitOverThresholdCountMap[detID][type]) {
  //             meanHItOverThresholdCount += hitOverThresholdCount;
  //         }
  //         meanHItOverThresholdCount /= m_stripnum;//strip数量
  //         for(int stripID = 0; stripID < m_stripnum; ++stripID) {
  //             auto it_strip = tempData[detID][type].find(stripID);
  //             if (it_strip == tempData[detID][type].end())
  //                 continue;
  //             const auto& baselineList = it_strip->second;
  //             if (baselineList.empty())
  //                 continue;
  //             ChannelNoise cn;
  //             cn.detectorID = detID;
  //             cn.stripID = stripID;
  //             cn.type = type;
  //             cn.nSamples = baselineList.size();
  //             cn.baselineList = baselineList;  // 保存原始数据用于绘制直方图

  //             // 使用 AnalysisUtils 计算平均基线和RMS
  //             cn.meanBaseline = AnalysisUtils::CalculateMean(baselineList);
  //             cn.rmsBaseline = AnalysisUtils::CalculateRMS(baselineList);

  //             // meanNoise 和 maxNoise 在这个实现中设为基线的RMS
  //             cn.meanNoise = cn.rmsBaseline;
  //             cn.maxNoise = cn.rmsBaseline;

  //             // 判断通道状态
  //             if (hitOverThresholdCountMap[detID][type][stripID] <
  //             meanHItOverThresholdCount * 0.2) {
  //                 cn.isBadChannel = true;  // 异常通道
  //             } else {
  //                 cn.isBadChannel = false; // 正常通道
  //             }

  //             m_channelNoiseMap[detID][type][stripID] = cn;

  //             totalChannels++;
  //         }
  //     }
  // }

  // for (const auto& [detID, typeMap] : tempData) {
  //     for (const auto& [type, stripMap] : typeMap) {
  //         int meanHItOverThresholdCount = 0;
  //         for(const auto& [stripID, hitOverThresholdCount] :
  //         hitOverThresholdCountMap[detID][type]) {
  //             meanHItOverThresholdCount += hitOverThresholdCount;
  //         }
  //         meanHItOverThresholdCount /= stripMap.size();
  //         for (const auto& [stripID, baselineList] : stripMap) {
  //             if (baselineList.empty())
  //                 continue;

  //             ChannelNoise cn;
  //             cn.detectorID = detID;
  //             cn.stripID = stripID;
  //             cn.type = type;
  //             cn.nSamples = baselineList.size();
  //             cn.baselineList = baselineList;  // 保存原始数据用于绘制直方图

  //             // 使用 AnalysisUtils 计算平均基线和RMS
  //             cn.meanBaseline = AnalysisUtils::CalculateMean(baselineList);
  //             cn.rmsBaseline = AnalysisUtils::CalculateRMS(baselineList);

  //             // meanNoise 和 maxNoise 在这个实现中设为基线的RMS
  //             cn.meanNoise = cn.rmsBaseline;
  //             cn.maxNoise = cn.rmsBaseline;

  //             // 判断通道状态
  //             if (hitOverThresholdCountMap[detID][type][stripID] <
  //             meanHItOverThresholdCount * 0.2) {
  //                 cn.isBadChannel = true;  // 异常通道
  //             } else {
  //                 cn.isBadChannel = false; // 正常通道
  //             }

  //             m_channelNoiseMap[detID][type][stripID] = cn;

  //             totalChannels++;
  //         }
  //     }
  // }
}

bool BadStripScript::SaveResults() {
  std::string outputFile = GetOutputDir() + "BadStripAnalysis.root";
  std::cout << "\nSaving results to: " << outputFile << std::endl;

  TFile *file = new TFile(outputFile.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot create output file: " << outputFile
              << std::endl;
    return false;
  }

  // 创建TTree保存通道噪声统计
  TTree *tree = new TTree("ChannelNoise", "Channel Noise Statistics");

  Int_t detectorID, stripID, type, nSamples;
  Double_t meanBaseline, rmsBaseline, meanNoise, maxNoise;
  Bool_t isBadChannel;

  tree->Branch("detectorID", &detectorID);
  tree->Branch("stripID", &stripID);
  tree->Branch("type", &type);
  tree->Branch("meanBaseline", &meanBaseline);
  tree->Branch("rmsBaseline", &rmsBaseline);
  tree->Branch("meanNoise", &meanNoise);
  tree->Branch("maxNoise", &maxNoise);
  tree->Branch("nSamples", &nSamples);
  tree->Branch("isBadChannel", &isBadChannel);

  // 填充TTree
  for (const auto &[detID, typeMap] : m_channelNoiseMap) {
    for (const auto &[tp, stripMap] : typeMap) {
      for (const auto &[strip, cn] : stripMap) {
        detectorID = cn.detectorID;
        stripID = cn.stripID;
        type = cn.type;
        meanBaseline = cn.meanBaseline;
        rmsBaseline = cn.rmsBaseline;
        meanNoise = cn.meanNoise;
        maxNoise = cn.maxNoise;
        nSamples = cn.nSamples;
        isBadChannel = cn.isBadChannel;
        tree->Fill();
      }
    }
  }

  tree->Write();

  // 创建每个探测器的通道噪声分布直方图
  auto &factory = DetectorFactory::GetInstance();
  const auto &detectors = factory.GetAllDetectors();

  for (const auto &[detID, det] : detectors) {
    if (m_channelNoiseMap.find(detID) == m_channelNoiseMap.end())
      continue;

    file->mkdir(Form("Detector_%d", detID));
    file->cd(Form("Detector_%d", detID));

    // 为X和Y方向分别创建直方图
    for (int type = 0; type < 2; ++type) {
      if (m_channelNoiseMap[detID].find(type) == m_channelNoiseMap[detID].end())
        continue;

      const auto &stripMap = m_channelNoiseMap[detID][type];

      // 找出strip范围
      int minStrip = 9999, maxStrip = -1;
      for (const auto &[strip, cn] : stripMap) {
        if (strip < minStrip)
          minStrip = strip;
        if (strip > maxStrip)
          maxStrip = strip;
      }

      int nBins = maxStrip - minStrip + 1;
      std::string typeName = (type == 0) ? "X" : "Y";

      // 平均噪声 vs 通道号
      TH1D *hNoise = new TH1D(
          Form("hNoise_%s", typeName.c_str()),
          Form("Detector %d - %s Channels Mean Noise;Strip ID;Mean Noise (ADC)",
               detID, typeName.c_str()),
          nBins, minStrip - 0.5, maxStrip + 0.5);

      // 平均基线 vs 通道号
      TH1D *hBaseline = new TH1D(Form("hBaseline_%s", typeName.c_str()),
                                 Form("Detector %d - %s Channels Mean "
                                      "Baseline;Strip ID;Mean Baseline (ADC)",
                                      detID, typeName.c_str()),
                                 nBins, minStrip - 0.5, maxStrip + 0.5);

      // 基线RMS vs 通道号
      TH1D *hBaselineRMS = new TH1D(Form("hBaselineRMS_%s", typeName.c_str()),
                                    Form("Detector %d - %s Channels Baseline "
                                         "RMS;Strip ID;Baseline RMS (ADC)",
                                         detID, typeName.c_str()),
                                    nBins, minStrip - 0.5, maxStrip + 0.5);
      TH1D *hBadChannel = new TH1D(Form("hBadChannel_%s", typeName.c_str()),
                                   Form("Detector %d - %s Channels Bad Channel "
                                        "Flag;Strip ID;Is Bad Channel",
                                        detID, typeName.c_str()),
                                   nBins, minStrip - 0.5, maxStrip + 0.5);

      // 填充直方图
      for (const auto &[strip, cn] : stripMap) {
        hNoise->SetBinContent(strip - minStrip + 1, cn.meanNoise);
        hBaseline->SetBinContent(strip - minStrip + 1, cn.meanBaseline);
        hBaselineRMS->SetBinContent(strip - minStrip + 1, cn.rmsBaseline);
        hBadChannel->SetBinContent(strip - minStrip + 1,
                                   cn.isBadChannel ? 1 : 0);
      }

      hNoise->Write();
      hBaseline->Write();
      hBaselineRMS->Write();
      hBadChannel->Write();

      delete hNoise;
      delete hBaseline;
      delete hBaselineRMS;
      delete hBadChannel;

      // 为异常通道创建详细的基线分布直方图
      file->mkdir(Form("Detector_%d/BadChannels_%s", detID, typeName.c_str()));
      file->cd(Form("Detector_%d/BadChannels_%s", detID, typeName.c_str()));

      for (const auto &[strip, cn] : stripMap) {
        if (cn.isBadChannel) {
          // 为每个异常通道创建基线分布直方图
          TH1D *hChannelBaseline =
              new TH1D(Form("hBaseline_Det%d_%s_Strip%d", detID,
                            typeName.c_str(), strip),
                       Form("Detector %d - %s Strip %d Baseline "
                            "Distribution;Baseline (ADC);Entries",
                            detID, typeName.c_str(), strip),
                       100, -100, 100);
          hChannelBaseline->StatOverflows(kTRUE);

          // 填充原始基线数据
          for (double baseline : cn.baselineList) {
            hChannelBaseline->Fill(baseline);
          }
          std::cout << AnalysisUtils::CalculateMean(cn.baselineList) << " "
                    << hChannelBaseline->GetMean() << std::endl;
          std::cout << cn.baselineList.size() << " "
                    << AnalysisUtils::CalculateRMS(cn.baselineList) << " "
                    << hChannelBaseline->GetStdDev() << std::endl;

          std::cout << "\n  Created histogram for bad channel: Det" << detID
                    << " " << typeName << " Strip " << strip
                    << " (RMS=" << cn.rmsBaseline << ")" << std::flush;

          hChannelBaseline->Write();
          delete hChannelBaseline;
        }
      }
    }

    file->cd();
  }

  // 创建全局噪声分布直方图
  file->cd();
  TH1D *hNoiseDistribution =
      new TH1D("hNoiseDistribution",
               "Global Channel Noise Distribution;Mean Noise (ADC);Channels",
               100, 0, 20);

  TH1D *hBaselineDistribution = new TH1D(
      "hBaselineDistribution",
      "Global Channel Baseline Distribution;Mean Baseline (ADC);Channels", 100,
      0, 500);

  for (const auto &[detID, typeMap] : m_channelNoiseMap) {
    for (const auto &[tp, stripMap] : typeMap) {
      for (const auto &[strip, cn] : stripMap) {
        hNoiseDistribution->Fill(cn.meanNoise);
        hBaselineDistribution->Fill(cn.meanBaseline);
      }
    }
  }

  hNoiseDistribution->Write();
  hBaselineDistribution->Write();

  delete hNoiseDistribution;
  delete hBaselineDistribution;

  file->Close();
  delete file;

  std::cout << "Results saved successfully!" << std::endl;
  return true;
}

// 注册脚本
REGISTER_SCRIPT("BadStripScript", BadStripScript);
