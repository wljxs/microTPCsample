#pragma once
#include "Algorithm/algorithms/ClusterBuilder.h"
#include "Algorithm/algorithms/ClusterReconstructor.h"
#include "Algorithm/algorithms/WaveformProcessor.h"
#include "DataModel.h"
#include "Detector.h"
#include <map>
#include <memory>
#include <stdexcept>

class DetectorFrame {
   public:
    explicit DetectorFrame(const Detector& det)
        : m_det(det) {}

    // Data Interfaces
    void SetRawData(const std::vector<RawData>& raw) { m_raw = raw; }
    void SetStripPedestalStd(std::shared_ptr<const std::map<int, std::vector<float>>> pedestalStd) {
        m_stripPedestalStd = std::move(pedestalStd);
    }
    void AddRawData(const RawData& raw) { m_raw.push_back(raw); }

    double GetT0() const { return t0; }
    const std::vector<RawData>& Raw() const { return m_raw; }
    const std::vector<StripHit>& StripHits() const { return m_stripHits; }
    const std::vector<Cluster>& Clusters() const { return m_clusters; }
    const std::vector<Cluster> Clusters(int type) const {
        std::vector<Cluster> result;
        result.reserve(m_clusters.size());
        for (const auto& clu : m_clusters)
            if (clu.type == type)
                result.push_back(clu);
        return result;
    };
    const std::vector<LocalHit>& LocalHits() const { return m_localHits; }
    const std::vector<GlobalHit>& GlobalHits() const { return m_globalHits; }
    const std::map<int, std::vector<float>>* StripPedestalStd() const { return m_stripPedestalStd.get(); }


    const std::vector<StripHitDeconv>& StripHitsDeconv() const { return m_stripHitsDeconv; }
    const std::vector<ClusterDeconv>& ClustersDeconv() const { return m_clustersDeconv; }
    const std::vector<ClusterDeconv>& ClustersDeconvBad() const { return m_clustersDeconv_bad; }
    // 可修改数据访问接口（供算法使用）
    std::vector<StripHit>& GetMutableStripHits() { return m_stripHits; }
    std::vector<Cluster>& GetMutableClusters() { return m_clusters; }
    std::vector<StripHitDeconv>& GetMutableStripHitsDeconv() { return m_stripHitsDeconv; }
    std::vector<ClusterDeconv>& GetMutableClustersDeconv() { return m_clustersDeconv; }
    std::vector<ClusterDeconv>& GetMutableClustersDeconvBad() { return m_clustersDeconv_bad; }

    // 根据索引获取单个StripHit（带边界检查）
    const StripHit& GetStripHit(int index) const {
        if (index < 0 || index >= static_cast<int>(m_stripHits.size())) {
            throw std::out_of_range("StripHit index out of range");
        }
        return m_stripHits[index];
    }

    std::vector<const StripHit*> GetStripHitsFromCluster(const Cluster& cluster) const {
        std::vector<const StripHit*> result;
        for (int idx : cluster.stripHitIndices) {
            if (idx >= 0 && idx < static_cast<int>(m_stripHits.size()))
                result.push_back(&m_stripHits[idx]);
        }
        return result;
    }

    const RawData* GetRawFromStrip(const StripHit& sh) const {
        if (sh.rawIndices < 0 || sh.rawIndices >= static_cast<int>(m_raw.size())) return nullptr;
        return &m_raw[sh.rawIndices];
    }

    bool Process(double t0);
    bool Process();

    // Extract StripHit Info from raw Data
    bool AnalyzeRaw();
    // Clustering StripHit
    bool Clustering();
    // Reconstruction: 将 Cluster 匹配并转换为 LocalHit
    bool Reconstruct();

    bool TransformToGlobal();

    //反卷积算法
    bool AnalyzeRaw(double t0);

    void clear() {
        m_raw.clear();
        m_stripHits.clear();
        m_clusters.clear();
        m_localHits.clear();
        m_globalHits.clear();
        m_stripPedestalStd.reset();
    }

    // Acess
    const Detector& det() const { return m_det; }

   private:
    double t0;
    //   Raw Level (from DAQ)
    std::vector<RawData> m_raw;

    // pedestals.apv_pedstd: (type -> [strip1..N])
    // 使用 shared_ptr 避免每个 event 拷贝整张表
    std::shared_ptr<const std::map<int, std::vector<float>>> m_stripPedestalStd;

    //   Step 1. Strip Level
    std::vector<StripHit> m_stripHits;

    //   Step 2. Cluster Level
    std::vector<Cluster> m_clusters;

    //   Step 3. Local Reconstruction
    std::vector<LocalHit> m_localHits;

    //   Step 4. Global Coordinates
    std::vector<GlobalHit> m_globalHits;

    // Deconvolution结果，前面几个是用电荷中心法得到的，这个是单独的
    std::vector<StripHitDeconv> m_stripHitsDeconv;
    std::vector<ClusterDeconv> m_clustersDeconv;
    std::vector<ClusterDeconv> m_clustersDeconv_bad; // 只保存cluster里有坏条的cluster，便于后续分析坏条对microTPC重建的影响

   private:
    const Detector& m_det;
};
