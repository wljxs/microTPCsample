#include "DetectorFrame.h"

bool DetectorFrame::AnalyzeRaw() {
    if (m_raw.empty()) return false;

    auto waveformProcessor = m_det.GetAlgorithm<WaveformProcessor>("WaveformProcessor");
    // std::cout << "\nAnalyzing raw data for Detector " << m_det.GetID() << ", " << m_raw.size() << " raws" << std::endl;
    
    return waveformProcessor->Process(*this);
}

// Clustering StripHit
bool DetectorFrame::Clustering() {
    if (m_stripHits.empty()) return false;

    auto clusterBuilder = m_det.GetAlgorithm<ClusterBuilder>("ClusterBuilder");
    // std::cout << "\nClustering strip hits for Detector " << m_det.GetID() << ", " <<m_stripHits.size() << " strip hits" << std::endl;
    return clusterBuilder->Process(*this);
}

// Reconstruction: 重建 Cluster 位置并生成 LocalHit
bool DetectorFrame::Reconstruct() {
    if (m_clusters.empty()) return false;
    // Step 1: 调用 ClusterReconstructor 计算每个 Cluster 的 pos
    auto clusterRecon = m_det.GetAlgorithm<ClusterReconstructor>("ClusterReconstructor");
    clusterRecon->Process(*this);

    // Step 2: 调用Detector几何逻辑生成LocalHits
    m_localHits = m_det.CalcLocalHitsFromClusters(m_clusters);
    // std::cout << "\nReconstructed " << m_clusters.size() << " clusters hits for Detector " << m_det.GetID() << std::endl;
    // std::cout << "Generated " << m_localHits.size() << " local hits for Detector " << m_det.GetID() << std::endl;
    return !m_localHits.empty();
}

bool DetectorFrame::TransformToGlobal() {
    if (m_localHits.empty()) {return false;}
    // 预分配内存优化
    m_globalHits.reserve(m_localHits.size());

    // LocalHits -> GlobalHits
    for (const auto& lh : m_localHits) {
        m_globalHits.push_back(m_det.LocalToGlobal(lh.localPos));
    }
    // std::cout << 1 << std::endl;
    return true;
}

bool DetectorFrame::Process(double t0) {
    this->t0 = t0;
    return AnalyzeRaw() && Clustering() && Reconstruct() && TransformToGlobal();
}

bool DetectorFrame::Process() {
    return AnalyzeRaw() && Clustering() && Reconstruct() && TransformToGlobal();
}

//如果是反卷积算法的话，直接AnalyzeRaw就行了，里面有相关的分析逻辑。
bool DetectorFrame::AnalyzeRaw(double t0) {
    this->t0 = t0;
    if (m_raw.empty()) return false;

    auto waveformProcessor = m_det.GetAlgorithm<WaveformProcessor>("WaveformProcessor");
    // std::cout << "\nAnalyzing raw data for Detector " << m_det.GetID() << ", " << m_raw.size() << " raws" << std::endl;
    
    return waveformProcessor->Process(*this);
}
