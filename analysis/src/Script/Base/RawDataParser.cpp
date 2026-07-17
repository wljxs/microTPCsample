
#include "Script/Base/RawDataParser.h"
#include "Detector/DetectorFactory.h"
#include <array>
#include <iostream>

RawDataParser::RawDataParser(const std::string &rawFile) : m_rawFile(rawFile) {}

RawDataParser::~RawDataParser() {
  if (m_file) {
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }
}

// 该代码是为了读取和解析存储在ROOT文件中的原始数据，然后变成为直接用的东西。

bool RawDataParser::Initialize() {
  // 打开ROOT文件
  m_file = TFile::Open(m_rawFile.c_str(), "READ");
  if (!m_file || m_file->IsZombie()) {
    std::cerr << "[RawDataParser] Failed to open file: " << m_rawFile
              << std::endl;
    return false;
  }

  // 获取TTree
  m_tree = (TTree *)m_file->Get("raw");
  if (!m_tree) {
    std::cerr << "[RawDataParser] raw TTree not found in " << m_rawFile
              << std::endl;
    return false;
  }
  m_pedTree = (TTree *)m_file->Get("pedestals");
  if (!m_pedTree) {
    std::cerr << "[RawDataParser] pedestal TTree not found in " << m_rawFile
              << std::endl;
    return false;
  }

  // 设置分支地址
  if (m_tree) {
    m_tree->SetBranchAddress("apv_id", &m_apv_id);
    m_tree->SetBranchAddress("apv_ch", &m_apv_ch);
    m_tree->SetBranchAddress("mm_strip", &m_mm_strip);
    m_tree->SetBranchAddress("apv_q", &m_apv_q);
    m_tree->SetBranchAddress("apv_evt", &m_apv_evt);
    m_numOfEvents = m_tree->GetEntries();
  }

  if(!LoadPedestals()){
      std::cerr << "[RawDataParser] Failed to load pedestals from " << m_rawFile
              << std::endl;
      return false;
  }

  std::cout << "[RawDataParser] Initialized with " << m_file->GetName()
            << " containing " << m_numOfEvents << " events" << std::endl;
  return true;
}

bool RawDataParser::LoadPedestals() {
    if (!m_pedTree) {
        std::cerr << "[RawDataParser] Pedestal TTree not available" << std::endl;
        return false;
    }

  // pedestals 树里的分支通常是 std::vector<...>，ROOT 需要传入“指针的地址”
  // 也就是：std::vector<T>* ptr = nullptr; SetBranchAddress(name, &ptr);
  std::vector<unsigned int>* apv_id_temp = nullptr;
  std::vector<unsigned int>* apv_ch_temp = nullptr;
  std::vector<unsigned int>* mm_strip_temp = nullptr;
  std::vector<float>* stddev = nullptr;

  if (!m_pedTree->GetBranch("apv_id") || !m_pedTree->GetBranch("apv_ch") ||
    !m_pedTree->GetBranch("mm_strip") || !m_pedTree->GetBranch("apv_pedstd")) {
    std::cerr << "[RawDataParser] Pedestal branches missing (need apv_id, apv_ch, mm_strip, apv_pedstd)" << std::endl;
    return false;
  }

  m_pedTree->SetBranchAddress("apv_id", &apv_id_temp);
  m_pedTree->SetBranchAddress("apv_ch", &apv_ch_temp);
  m_pedTree->SetBranchAddress("mm_strip", &mm_strip_temp);
  m_pedTree->SetBranchAddress("apv_pedstd", &stddev);

  const Long64_t nEntries = m_pedTree->GetEntries();
  if (nEntries <= 0) {
    std::cerr << "[RawDataParser] pedestals tree is empty" << std::endl;
    return false;
  }

  m_pedestalMap.clear();
  long long filled = 0;

  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    m_pedTree->GetEntry(entry);

    if (!apv_id_temp || !apv_ch_temp || !mm_strip_temp || !stddev) {
      std::cerr << "[RawDataParser] Pedestal branch pointers are null at entry " << entry << std::endl;
      continue;
    }

    const size_t n = apv_id_temp->size();
    if (apv_ch_temp->size() != n || mm_strip_temp->size() != n || stddev->size() != n) {
      std::cerr << "[RawDataParser] Pedestal vector sizes mismatch at entry " << entry
            << " (apv_id=" << apv_id_temp->size()
            << ", apv_ch=" << apv_ch_temp->size()
            << ", mm_strip=" << mm_strip_temp->size()
            << ", apv_pedstd=" << stddev->size() << ")" << std::endl;
      continue;
    }

    for (size_t i = 0; i < n; ++i) {
      auto [detID, stripID, type] =
        MapBoardChannel((*apv_id_temp)[i], (*apv_ch_temp)[i], (*mm_strip_temp)[i]);

      if (stripID <= 0) {
        continue;
      }

      auto &vec = m_pedestalMap[detID][type];
      if (vec.size() < static_cast<size_t>(stripID)) {
        vec.resize(static_cast<size_t>(stripID), 0.0f);
      }
      vec[static_cast<size_t>(stripID - 1)] = (*stddev)[i];
      filled++;
    }
  }

  std::cout << "[RawDataParser] Loaded pedestal stddev for " << filled << " strips" << std::endl;
    return true;
}

std::vector<float> RawDataParser::GetPedstd(int detectorID,int type) const {
    if (m_pedestalMap.find(detectorID) != m_pedestalMap.end() &&
        m_pedestalMap.at(detectorID).find(type) != m_pedestalMap.at(detectorID).end()) {
        return m_pedestalMap.at(detectorID).at(type);
    } else {
        std::cerr << "[RawDataParser] Pedestal stddev not found for detID=" << detectorID
                  << ", type=" << type << std::endl;
        return {};
    }
}

std::unordered_map<int, std::vector<RawData>>
RawDataParser::LoadEvent(int eventID) {
  std::unordered_map<int, std::vector<RawData>> result;

  if (!m_tree) {
    std::cerr << "[RawDataParser] ERROR: TTree not initialized\n";
    return result;
  }

  if (eventID < 0 || eventID >= GetTotalEvents()) {
    std::cerr << "[RawDataParser] ERROR: Invalid event index " << eventID
              << "\n";
    return result;
  }

  // ---- Load TTree Entry ----
  m_tree->GetEntry(eventID);

  const size_t nHits = m_apv_id->size();
  result.reserve(16);

  auto &factory = DetectorFactory::GetInstance();
  const auto &detectors = factory.GetAllDetectors();

  // ---- Loop all APV hits ----
  for (size_t j = 0; j < nHits; ++j) {

    auto [detID, stripID, type] =
        MapBoardChannel((*m_apv_id)[j], (*m_apv_ch)[j], (*m_mm_strip)[j]);

    if (detectors.find(detID) == detectors.end())
      continue;

    RawData raw{stripID, type, (*m_apv_q)[j]};

    result[detID].push_back(std::move(raw));
  }

  return result;
}

// 硬件映射常量
constexpr std::array<int, 16> kBoardToRawIndex = {0, 0, 1, 1, 2, 2, 3, 3,
                                                  4, 4, 5, 5, 6, 6, 7, 7};

std::tuple<int, int, int>
RawDataParser::MapBoardChannel(unsigned int boardID, unsigned int channelID,
                               unsigned int mm_strip) const {//这里面如果有探测器还需要设置一下

  int rawDataIndex = (boardID < kBoardToRawIndex.size())
                         ? kBoardToRawIndex[boardID]
                         : static_cast<int>(boardID) / 2;

  int type = (rawDataIndex % 2 == 0) ? 0 : 1;
  int detID = (rawDataIndex / 2) + 1;
  int stripID = static_cast<int>(mm_strip);

//   // 1884
//   if (boardID == 12)
//     stripID = 256 - channelID;
//   else if (boardID == 13)
//     stripID = 128 - channelID;
//   else if (boardID == 14)
//     stripID = 256 - channelID;
//   else if (boardID == 15)
//     stripID = 128 - channelID;
 //1700                               
  if (boardID == 12)
    stripID = channelID + 1;
  else if (boardID == 13)
    stripID = 129 + channelID;
  else if (boardID == 14)
    stripID = 256 - channelID;
  else if (boardID == 15)
    stripID = 128 - channelID;
  // 1978
  // if (boardID == 14) {
  //     if (channelID % 2 == 0)
  //         stripID = channelID + 2;
  //     else
  //         stripID = channelID;
  // } else if (boardID == 15) {
  //     if (channelID % 2 == 0)
  //         stripID = 130 + channelID;
  //     else
  //         stripID = 128 + channelID;
  // }

  return {detID, stripID, type};
}
