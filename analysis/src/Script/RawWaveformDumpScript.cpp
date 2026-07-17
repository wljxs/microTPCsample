#include "Script/RawWaveformDumpScript.h"

#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <unordered_set>

REGISTER_SCRIPT("RawWaveformDump", RawWaveformDumpScript);

namespace {

bool containsInt(const std::vector<int>& v, int x) {
    return std::find(v.begin(), v.end(), x) != v.end();
}

struct WaveformMetrics {
    double baseline = 0.0;
    double peak = 0.0;          // absolute (ADC)
    int peakIdx = -1;

    double peakAboveBaseline = 0.0;

    // FWHM based on (peak-baseline)/2
    double halfHeight = 0.0;    // value above baseline
    int leftIdx = -1;
    int rightIdx = -1;
    int widthSamples = -1;

    int failMask = 0;

    static constexpr int kFailEmpty = 1 << 0;
    static constexpr int kFailBelowNoise = 1 << 1;
    static constexpr int kFailFwhm = 1 << 2;
};

WaveformMetrics analyzeWaveform(const std::vector<short>& adc, int baselineNSamples, double noiseTh) {
    WaveformMetrics m;

    if (adc.empty()) {
        m.failMask |= WaveformMetrics::kFailEmpty;
        return m;
    }

    const int n = static_cast<int>(adc.size());
    const int nb = std::max(1, std::min(baselineNSamples, n));

    double bsum = 0.0;
    for (int i = 0; i < nb; ++i) bsum += static_cast<double>(adc[i]);
    m.baseline = bsum / nb;

    m.peak = static_cast<double>(adc[0]);
    m.peakIdx = 0;
    for (int i = 1; i < n; ++i) {
        const double v = static_cast<double>(adc[i]);
        if (v > m.peak) {
            m.peak = v;
            m.peakIdx = i;
        }
    }

    m.peakAboveBaseline = m.peak - m.baseline;
    if (m.peakAboveBaseline < noiseTh) {
        m.failMask |= WaveformMetrics::kFailBelowNoise;
    }

    // FWHM: find first indices on both sides where (adc-baseline) drops below halfHeight
    m.halfHeight = 0.5 * m.peakAboveBaseline;
    if (m.halfHeight <= 0.0) {
        m.failMask |= WaveformMetrics::kFailFwhm;
        return m;
    }

    int left = -1;
    for (int i = m.peakIdx; i >= 0; --i) {
        if ((static_cast<double>(adc[i]) - m.baseline) < m.halfHeight) {
            left = i;
            break;
        }
    }

    int right = -1;
    for (int i = m.peakIdx; i < n; ++i) {
        if ((static_cast<double>(adc[i]) - m.baseline) < m.halfHeight) {
            right = i;
            break;
        }
    }

    m.leftIdx = left;
    m.rightIdx = right;

    if (left < 0 || right < 0 || left >= right) {
        m.failMask |= WaveformMetrics::kFailFwhm;
    } else {
        m.widthSamples = right - left;
    }

    return m;
}

} // namespace

void RawWaveformDumpScript::LoadConfig(const json& config) {
    m_startEvent = config.value("startEvent", 0);
    m_maxEvents = config.value("maxEvents", 1000);
    m_progressInterval = config.value("progressInterval", 1000);

    if (config.contains("detIDs") && config["detIDs"].is_array()) {
        m_detIDs = config["detIDs"].get<std::vector<int>>();
    }
    if (config.contains("stripTypes") && config["stripTypes"].is_array()) {
        m_stripTypes = config["stripTypes"].get<std::vector<int>>();
    }
    if (config.contains("stripIDs") && config["stripIDs"].is_array()) {
        m_stripIDs = config["stripIDs"].get<std::vector<int>>();
    }

    m_baselineNSamples = config.value("baselineNSamples", 1);
    m_noiseThreshold = config.value("noiseThreshold", 0.0);
    m_dumpOnlyFailed = config.value("dumpOnlyFailed", false);

    m_outputFile = config.value("outputFile", std::string("RawWaveforms.root"));
}

void RawWaveformDumpScript::Print() const {
    std::cout << "[" << GetName() << "]" << std::endl;
    std::cout << "  startEvent: " << m_startEvent << std::endl;
    std::cout << "  maxEvents: " << m_maxEvents << std::endl;
    std::cout << "  progressInterval: " << m_progressInterval << std::endl;
    std::cout << "  baselineNSamples: " << m_baselineNSamples << std::endl;
    std::cout << "  noiseThreshold: " << m_noiseThreshold << std::endl;
    std::cout << "  dumpOnlyFailed: " << (m_dumpOnlyFailed ? "true" : "false") << std::endl;
    std::cout << "  outputFile: " << (GetOutputDir() + m_outputFile) << std::endl;

    if (!m_detIDs.empty()) {
        std::cout << "  detIDs filter size: " << m_detIDs.size() << std::endl;
    }
    if (!m_stripTypes.empty()) {
        std::cout << "  stripTypes filter size: " << m_stripTypes.size() << std::endl;
    }
    if (!m_stripIDs.empty()) {
        std::cout << "  stripIDs filter size: " << m_stripIDs.size() << std::endl;
    }
}

bool RawWaveformDumpScript::passFilter(int detID, int stripType, int stripID) const {
    if (!m_detIDs.empty() && !containsInt(m_detIDs, detID)) return false;
    if (!m_stripTypes.empty() && !containsInt(m_stripTypes, stripType)) return false;
    if (!m_stripIDs.empty() && !containsInt(m_stripIDs, stripID)) return false;
    return true;
}

bool RawWaveformDumpScript::Execute() {
    auto parser = GetParser();
    if (!parser) {
        std::cerr << "[RawWaveformDump] parser is null" << std::endl;
        return false;
    }

    const Long64_t totalEvents = parser->GetTotalEvents();
    if (totalEvents <= 0) {
        std::cerr << "[RawWaveformDump] totalEvents <= 0" << std::endl;
        return false;
    }

    int start = std::max(0, m_startEvent);
    if (start >= totalEvents) {
        std::cerr << "[RawWaveformDump] startEvent out of range: " << start << std::endl;
        return false;
    }

    int endExclusive;
    if (m_maxEvents <= 0) {
        endExclusive = static_cast<int>(totalEvents);
    } else {
        endExclusive = std::min(static_cast<int>(totalEvents), start + m_maxEvents);
    }

    const std::string outPath = GetOutputDir() + m_outputFile;
    TFile* fout = TFile::Open(outPath.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[RawWaveformDump] cannot create output file: " << outPath << std::endl;
        return false;
    }

    // Tree branches
    int b_eventID = 0;
    int b_detID = 0;
    int b_stripType = 0;
    int b_stripID = 0;

    int b_nSamples = 0;
    std::vector<short> b_adc;

    double b_baseline = 0.0;
    double b_peak = 0.0;
    int b_peakIdx = -1;
    double b_peakAboveBaseline = 0.0;

    double b_halfHeight = 0.0;
    int b_leftIdx = -1;
    int b_rightIdx = -1;
    int b_widthSamples = -1;

    int b_failMask = 0;

    TTree* t = new TTree("RawWaveforms", "Raw waveforms (adc) with simple metrics");
    t->Branch("eventID", &b_eventID);
    t->Branch("detID", &b_detID);
    t->Branch("stripType", &b_stripType);
    t->Branch("stripID", &b_stripID);
    t->Branch("nSamples", &b_nSamples);
    t->Branch("adc", &b_adc);

    t->Branch("baseline", &b_baseline);
    t->Branch("peak", &b_peak);
    t->Branch("peakIdx", &b_peakIdx);
    t->Branch("peakAboveBaseline", &b_peakAboveBaseline);

    t->Branch("halfHeight", &b_halfHeight);
    t->Branch("leftIdx", &b_leftIdx);
    t->Branch("rightIdx", &b_rightIdx);
    t->Branch("widthSamples", &b_widthSamples);

    t->Branch("failMask", &b_failMask);

    long long filled = 0;

    for (int evt = start; evt < endExclusive; ++evt) {
        if (m_progressInterval > 0 && ((evt - start) % m_progressInterval) == 0) {
            std::cout << "\r[RawWaveformDump] event " << evt << "/" << totalEvents
                      << " filled=" << filled << std::flush;
        }

        auto rawMap = parser->LoadEvent(evt);
        if (rawMap.empty()) continue;

        for (const auto& [detID, raws] : rawMap) {
            for (const auto& raw : raws) {
                if (!passFilter(detID, raw.type, raw.stripID)) continue;

                const auto m = analyzeWaveform(raw.adc, m_baselineNSamples, m_noiseThreshold);
                const bool isFailed = (m.failMask != 0);
                if (m_dumpOnlyFailed && !isFailed) continue;

                b_eventID = evt;
                b_detID = detID;
                b_stripType = raw.type;
                b_stripID = raw.stripID;

                b_nSamples = static_cast<int>(raw.adc.size());
                b_adc = raw.adc;

                b_baseline = m.baseline;
                b_peak = m.peak;
                b_peakIdx = m.peakIdx;
                b_peakAboveBaseline = m.peakAboveBaseline;

                b_halfHeight = m.halfHeight;
                b_leftIdx = m.leftIdx;
                b_rightIdx = m.rightIdx;
                b_widthSamples = m.widthSamples;

                b_failMask = m.failMask;

                t->Fill();
                ++filled;
            }
        }
    }

    std::cout << "\n[RawWaveformDump] Done. filled=" << filled << " output=" << outPath << std::endl;

    fout->cd();
    t->Write();
    fout->Close();
    delete fout;

    return true;
}
