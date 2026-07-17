#pragma once

#include "TVector3.h"
#include <map>
#include <memory>

class DetectorFrame;

struct RawData {
    int stripID;             // 读出板编号
    int type;                // 通道编号，就是x/y之类的
    std::vector<short> adc;  // 波形采样点 (ADC counts)
};

struct StripHit {
    int ID;    // 条号
    int type;  // 读出条类型 (X, Y, U/V…)
    int rawIndices;

    // ---- 提取的关键量 ----
    double amp;       // 峰值 - baseline
    double charge;    // 积分电荷
    int peakTime;     // 峰值时间 (采样点索引)
    double time;      // 信号时间 (前沿拟合+CFD)
    double riseTime;  // 上升时间 (10%-90%)
    double width;     // 信号宽度
    int baseline;    // 基线水平 (ADC counts)

    // ---- 误差信息 ----
    double timeError;  // 时间误差

    // ---- 标志位 ----
    bool isSaturated;  // 是否饱和
    bool isValid;      // 是否有效

    // ---- Bad strip flag ----
    bool isBadChannel = false;  // 是否坏道（用于屏蔽）
};

struct Cluster {
    int type;                          // 属于X/Y/U/V 哪个方向
    std::vector<int> stripHitIndices;  // StripHit在DetectorFrame::m_stripHits中的全局索引

    // ---- 聚类整体量 ----
    int size;       // 聚类条数
    int range;      // 聚类范围 (最大ID - 最小ID)
    double charge;  // 聚类总电荷
    double maxAmp;
    double time;      // 聚类时间 (最早)
    double centroid;  // 聚类重心
    double pos;       // cluster重建位置
    bool isBad;       // 是否包含坏道（用于屏蔽）
};

typedef TVector3 GlobalHit;

struct LocalHit {
    TVector3 localPos;
    std::vector<int> clusterIndices;
    bool isValid = true;
};

// 径迹数据结构
struct Track {
    double kx;    // x方向斜率
    double ky;    // y方向斜率
    double bx;    // x截距
    double by;    // y截距
    double chi2;  // 拟合质量
};

struct Event {
    int eventID;
    std::map<int, std::shared_ptr<DetectorFrame>> detectorFramesMap;
    Track track;
    double t0;
};

struct StripHitDeconv{//以后要是有速度问题，可以想想怎么弄
    int ID;   //条号
    int type; //读出条类型
    int rawIndices;
    bool isBad = false; //是否坏道
    std::vector<double> chargetier; //不同层的电荷量
    double FwHM; //半高宽
    double riseTime; //上升时间
    int peakTime; //达峰时间
    int Integral; //积分电荷
    int amp; //振幅
    double overthresholdTime; //过阈时刻
    double rnorm; //残差范数
    double rnorm_no; //没有归一化的残差范数
    double chi2; //拟合的χ²
    double t0;
    double scale; //拟合的缩放因子
    double r2; //拟合的R²
};

struct ClusterDeconv{//一下的斜率，截距的y方向都是基于条的
    int type;
    int size;
    double energy;
    double maxAmp;
    double k;
    double b;
    double chi2;
    double chi2_re; //   chi2/ndf
    std::vector<int> stripHitIndices;
    std::vector<int> stripIDs;
    std::vector<double> chargetier;//统计每一层的电荷
    int tier_hascharge; //有电荷的层数
    std::vector<double> pos; //统计每一层的重心位置
    bool isValid = false;
    bool isBad = false;
    double chargeposition; //电荷重心位置
    double microTPCposition; //得到的打到微TPC的重心位置
};