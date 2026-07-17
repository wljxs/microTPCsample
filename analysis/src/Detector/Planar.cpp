#include "Detector/Planar.h"
#include "DataModel.h"
#include "TVector3.h"
#include "iostream"
#include <cmath>
#include <map>
#include <vector>

Planar::Planar(int id, const std::string& name, const json& config) : Detector(id, name, config) {

    if (config.contains("readoutPlanes")) {
        m_config.readoutPlaneAngle.clear();
        m_config.readoutPlanePitch.clear();
        m_config.readoutPlaneStripNumber.clear();
        for (auto& [planeIDStr, planeData] : config["readoutPlanes"].items()) {
            try {
                int type = std::stoi(planeIDStr);
                m_config.readoutPlaneAngle.emplace(type, planeData["angle"]);
                m_config.readoutPlanePitch.emplace(type, planeData["pitch"]);
                m_config.readoutPlaneStripNumber.emplace(type, int(planeData["num"]));
            } catch (std::exception& e) {
                std::cout << "Error: " << e.what() << " in readout plane configuration" << std::endl;
            }
        }
    }
}

GlobalHit Planar::CalcHitFromTrack(const Track& track) const {
    GlobalHit globalHit;

    TVector3 rot = GetRot();
    double rx = rot.X(), ry = rot.Y(), rz = rot.Z();

    double cx = cos(rx), sx = sin(rx);
    double cy = cos(ry), sy = sin(ry);
    double cz = cos(rz), sz = sin(rz);

    // 正确的法向量计算（旋转矩阵的第三列），因为初始的时候认为是z轴方向
    TVector3 normal(cz * sy * cx + sz * sx,
                    sz * sy * cx - cz * sx,
                    cy * cx);
    normal = normal.Unit();

    // 更合理的轨迹参数化
    TVector3 trackOrigin(track.bx, track.by, 0);

    // 假设track包含完整的方向向量
    TVector3 trackDirection(track.kx, track.ky, 1);
    trackDirection = trackDirection.Unit();

    TVector3 planePoint = GetPos();
    TVector3 diff = planePoint - trackOrigin;
    double denominator = trackDirection.Dot(normal);

    if (fabs(denominator) < 1e-12) {//通过说明track与平面平行
        return globalHit;  // 平行，无交点
    }

    double t = diff.Dot(normal) / denominator;
    if (t < 0) {
        return globalHit;  // 反向，无交点
    }

    globalHit = trackOrigin + t * trackDirection;
    return globalHit;
}


std::pair<std::vector<GlobalHit>,std::vector<TVector3>> Planar::CalcHitsFromTrack(const Track& track,int ntier,double thickness){
    std::vector<GlobalHit> hits;
    std::vector<TVector3> localHits;
    Track trackAtZ = track;    
    TVector3 rot = GetRot();
    double rx = rot.X(), ry = rot.Y(), rz = rot.Z();
    double cx = cos(rx), sx = sin(rx);
    double cy = cos(ry), sy = sin(ry);
    double cz = cos(rz), sz = sin(rz);

    TVector3 normal(cz * sy * cx + sz * sx,
                    sz * sy * cx - cz * sx,
                    cy * cx);
    normal = normal.Unit();
    TVector3 trackOrigin(track.bx, track.by, 0);
    TVector3 trackDirection(track.kx, track.ky, 1);
    trackDirection = trackDirection.Unit();

    TVector3 planePoint = GetPos();
    double denominator = trackDirection.Dot(normal);
    if (fabs(denominator) < 1e-12) {
        return std::pair<std::vector<GlobalHit>,std::vector<TVector3>>(std::vector<GlobalHit>(), std::vector<TVector3>());  // 平行，无交点
    }

    // 逆旋转矩阵 R^T，用于全局坐标 -> 局域坐标
    // R = Rz(rz) * Ry(ry) * Rx(rx), R^T = 转置
    double r00 = cz * cy;
    double r01 = sz * cy;
    double r02 = -sy;
    double r10 = cz * sy * sx - sz * cx;
    double r11 = sz * sy * sx + cz * cx;
    double r12 = cy * sx;
    // r20, r21, r22 即 normal 的三个分量

    for(int i=0;i<ntier;i++){
        double z_ref_position = thickness / ntier * (i + 0.5) ;
        TVector3 refPoint = planePoint + normal * z_ref_position;
        TVector3 diff = refPoint - trackOrigin;
        double t = diff.Dot(normal) / denominator;
        if (t < 0) {
            continue;  // 反向，无交点
        }
        GlobalHit hit = trackOrigin + t * trackDirection;
        hits.push_back(hit);

        // 全局坐标 -> 局域坐标: local = R^T * (hit - pos)
        TVector3 p = hit - refPoint;
        TVector3 local(
            r00 * p.X() + r01 * p.Y() + r02 * p.Z(),
            r10 * p.X() + r11 * p.Y() + r12 * p.Z(),
            normal.X() * p.X() + normal.Y() * p.Y() + normal.Z() * p.Z()
        );
        localHits.push_back(local);
    }
    return std::pair<std::vector<GlobalHit>,std::vector<TVector3>>(hits, localHits);
}

std::vector<LocalHit> Planar::CalcLocalHitsFromClusters(const std::vector<Cluster>& clusters) const {
    std::vector<LocalHit> localHits;
    if (clusters.empty()) {
        std::cout << "No clusters to process." << std::endl; 
        return localHits;}
    
    // 按type分组clusters
    std::map<int, std::vector<size_t>> clustersByType;  // type -> cluster indices
    //用cluster的maxAmp先拍个序
    std::vector<Cluster> sortedClusters = clusters;
    std::sort(sortedClusters.begin(), sortedClusters.end(), [](const Cluster& a, const Cluster& b) {
        return a.maxAmp > b.maxAmp;  // 按maxAmp降序排序
    });


    for (size_t i = 0; i < sortedClusters.size(); ++i) {
        clustersByType[sortedClusters[i].type].push_back(i);
    }

    // 获取所有type
    std::vector<int> types = m_config.readoutPlaneType;

    // 根据type数量决定匹配策略
    if (types.size() == 1) {
        // 单type：每个cluster直接生成LocalHit
        int type = types[0];
        for (size_t idx : clustersByType[type]) {
            LocalHit lh;
            if (type == 0) {
                if(clusters[idx].isBad) lh.isValid = false; // 如果cluster包含坏道，则认为这个hit无效
                lh.localPos.SetXYZ(sortedClusters[idx].pos * m_config.readoutPlanePitch.at(type), 0, 0);
            } else if (type == 1) {
                if(clusters[idx].isBad) lh.isValid = false; // 如果cluster包含坏道，则认为这个hit无效
                lh.localPos.SetXYZ(0, sortedClusters[idx].pos * m_config.readoutPlanePitch.at(type), 0);
            }
            lh.clusterIndices.push_back(idx);
            localHits.push_back(lh);
        }
    } else if (types.size() == 2) {
        // 双type：X-Y匹配
        int type0 = types[0];
        int type1 = types[1];

        const auto& indices0 = clustersByType[type0];
        const auto& indices1 = clustersByType[type1];

        // 简单匹配：所有组合（可根据时间、电荷等优化）
        for (size_t idx0 : indices0) {
            for (size_t idx1 : indices1) {
                LocalHit lh;
                double x = 0, y = 0;

                if (type0 == 0) {
                    if(clusters[idx0].isBad || clusters[idx1].isBad) lh.isValid = false; // 如果任一cluster包含坏道，则认为这个hit无效
                    x = sortedClusters[idx0].pos * m_config.readoutPlanePitch.at(type0);
                    y = sortedClusters[idx1].pos * m_config.readoutPlanePitch.at(type1);
                } else {
                    if(clusters[idx0].isBad || clusters[idx1].isBad) lh.isValid = false; // 如果任一cluster包含坏道，则认为这个hit无效
                    x = sortedClusters[idx1].pos * m_config.readoutPlanePitch.at(type1);
                    y = sortedClusters[idx0].pos * m_config.readoutPlanePitch.at(type0);
                }

                lh.localPos.SetXYZ(x, y, 0);
                lh.clusterIndices.push_back(idx0);
                lh.clusterIndices.push_back(idx1);
                localHits.push_back(lh);
            }
        }
    } else {
        // 多type：复杂匹配逻辑（暂不实现）
        std::cerr << "[Planar] Warning: More than 2 types not supported yet" << std::endl;
    }
    return localHits;
}
