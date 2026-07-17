#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

/**
 * @brief Bad strip mask for a single detector.
 *
 * Loaded from a JSON file produced by BadStripScript.
 * Supports two schemas:
 * 1) { "badStrips": { "<detID>": { "<type>": [stripIDs...] } } }
 * 2) { "badStrips": [ {"detectorID":1,"type":0,"strips":[...]}, ... ] }
 */
class BadStripMask {
   public:
    BadStripMask() = default;

    bool IsBad(int type, int stripID) const {
        auto it = m_badByType.find(type);
        if (it == m_badByType.end()) return false;
        return it->second.find(stripID) != it->second.end();
    }

    size_t Count(int type) const {
        auto it = m_badByType.find(type);
        if (it == m_badByType.end()) return 0;
        return it->second.size();
    }

    static std::shared_ptr<BadStripMask> LoadForDetector(const std::string& filePath, int detectorID) {
        std::ifstream in(filePath);
        if (!in.is_open()) {
            std::cerr << "[BadStripMask] Cannot open bad strips file: " << filePath << std::endl;
            return nullptr;
        }

        json j;
        try {
            in >> j;
        } catch (const std::exception& e) {
            std::cerr << "[BadStripMask] Failed to parse JSON '" << filePath << "': " << e.what() << std::endl;
            return nullptr;
        }

        auto mask = std::make_shared<BadStripMask>();

        if (!j.contains("badStrips")) {
            std::cerr << "[BadStripMask] JSON has no 'badStrips' field: " << filePath << std::endl;
            return nullptr;
        }

        const auto& bs = j["badStrips"];
        // Schema 1: object keyed by detector ID
        if (bs.is_object()) {
            const std::string detKey = std::to_string(detectorID);
            if (!bs.contains(detKey)) {
                return mask;  // valid but empty for this detector
            }
            const auto& detNode = bs[detKey];
            if (!detNode.is_object()) return mask;
            for (auto it = detNode.begin(); it != detNode.end(); ++it) {
                int type = 0;
                try {
                    type = std::stoi(it.key());
                } catch (...) {
                    continue;
                }
                const auto& strips = it.value();
                if (!strips.is_array()) continue;
                for (const auto& sid : strips) {
                    if (!sid.is_number_integer()) continue;
                    mask->m_badByType[type].insert(sid.get<int>());
                }
            }
            return mask;
        }

        // Schema 2: array of entries
        if (bs.is_array()) {
            for (const auto& entry : bs) {
                if (!entry.is_object()) continue;
                if (!entry.contains("detectorID") || !entry.contains("type") || !entry.contains("strips")) continue;
                if (entry["detectorID"].get<int>() != detectorID) continue;
                const int type = entry["type"].get<int>();
                const auto& strips = entry["strips"];
                if (!strips.is_array()) continue;
                for (const auto& sid : strips) {
                    if (!sid.is_number_integer()) continue;
                    mask->m_badByType[type].insert(sid.get<int>());
                }
            }
            return mask;
        }

        return mask;
    }

   private:
    std::unordered_map<int, std::unordered_set<int>> m_badByType;
};
