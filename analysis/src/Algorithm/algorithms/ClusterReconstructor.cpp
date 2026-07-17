#include "algorithms/ClusterReconstructor.h"
#include "AlgorithmFactory.h"
#include "DataModel.h"
#include "DetectorFrame.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <cmath>

REGISTER_ALGORITHM("ClusterReconstructor", ClusterReconstructor)

bool ClusterReconstructor::Process(DetectorFrame& frame) {
    auto& clusters = frame.GetMutableClusters();
    if (clusters.empty()) return false;

    const auto& stripHits = frame.StripHits();

    // 遍历每个Cluster，调用内部重建逻辑更新pos字段
    for (auto& cluster : clusters) {
        reconstructChargeWeighted(cluster, stripHits);
    }
    return true;
}

void ClusterReconstructor::reconstructChargeWeighted(Cluster& cluster, const std::vector<StripHit>& stripHits) {
    if (cluster.stripHitIndices.empty()) {
        cluster.pos = 0.0;
        return;
    }

    // 电荷加权位置重建（以stripID为单位）
    double weightedSum = 0.0;
    double totalCharge = 0.0;

    for (int idx : cluster.stripHitIndices) {
        const auto& strip = stripHits[idx];
        if (strip.isBadChannel) continue;
        double weight = pow(strip.amp, 1);
        weightedSum += strip.ID * weight;
        totalCharge += weight;
    }

    if (totalCharge > 0) {
        cluster.pos = weightedSum / totalCharge;
    } else {
        int centerIdx = cluster.size / 2;
        cluster.pos = stripHits[cluster.stripHitIndices[centerIdx]].ID;
    }
}