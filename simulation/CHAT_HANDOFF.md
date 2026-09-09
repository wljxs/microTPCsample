# simulation 新聊天交接

新聊天开始后先读这份文件，再查看实际代码。这里只记录当前已经实现并验证的内容，
不把其他目录和未确认的设想混进来。

## 当前流程

```text
config/apv25_1.jsonc
        |
        v
microTPCsimulation.C -> result/apv25_1/events.root
        |
        | 读取 simulation_metadata.tierCenter
        v
microTPCsample.C -> result/apv25_1/tiers/tierN.root
        |
        v
analysis/mean.C -> result/apv25_1/mean/meanN.root
        |
        v
analysis/matix.C -> result/apv25_1/matrix.root
```

`run.sh` 是一键入口。配置文件在脚本顶部选择，不继续增加复杂的命令行参数。

## 当前代码改动

### `Config.h`

建立了共用配置结构：

```cpp
struct Simulation {
    int runid;
    unsigned int nEvents;
    std::string detectorDir;
    std::string gasFile;
    std::string fieldFile;
    double theta;
    double phi;
    unsigned int seed;
    std::string mode;
    bool deltaElectron;
    double sampleInterval;  // ns，模拟采样间隔
    double gap;             // cm，漂移区上边界
};

struct Electronics {
    unsigned int n;
    double tau;             // ns，成形时间常数
    double sampleInterval;  // ns，电子学采样间隔
};

struct Config {
    std::string outputDir;
    double totalTime;       // ns
    Simulation simulation;
    Electronics electronics;
};
```

- 配置统一使用 JSONC，允许注释。
- 文件不存在时使用结构体中的默认值。
- 模拟采样间隔与电子学采样间隔分开。
- 总时间 `totalTime` 共用。
- 成形参数直接使用 `n` 和 `tau`。

### `config/apv25_1.jsonc`

当前配置为：

- 输出目录：`../result/apv25_1`
- 事例数：1000
- 探测器目录：`../microTPC5mmone`
- 气体文件：`ar_45_cf4_40_co2_15100-100000.gas`
- 电场文件：`field.txt`
- `theta = 20 degree`，`phi = 90 degree`
- `seed = 42`
- `mode = "center"`
- `deltaElectron = false`
- 模拟采样间隔：5 ns
- 漂移区 `gap = 0.5 cm`
- APV25：`n = 1`，`tau = 50 ns`，电子学采样间隔 25 ns
- 总时间：700 ns

### `microTPCsimulation.C`

- 读取共享 JSONC，输出 `<outputDir>/events.root`。
- 气体、探测器和电场文件来自配置。
- ROOT 与 Garfield 随机引擎都设置 `config.simulation.seed`。
- `mode` 控制入射位置，`deltaElectron` 控制 delta 电子输运。
- 使用 Garfield `Shaper(n, tau, ...)` 成形。
- 在漂移区中心取得电场，通过 Garfield `ElectronVelocity()` 得到 z 方向漂移速度。
- 根据漂移速度和电子学采样间隔确定层厚及层中心。
- 把配置、`gap`、`tierCenter` 和运行统计写入 `simulation_metadata`。

层中心计算为：

```cpp
const double dz =
    std::abs(vz) * config.electronics.sampleInterval;

std::vector<double> tierCenter;
for (unsigned int tier = 0; tier * dz < z0; ++tier) {
    const double lower = std::max(zBottom, tier * dz);
    const double upper = std::min((tier + 1) * dz, z0);
    if (upper > lower) {
        tierCenter.push_back(0.5 * (lower + upper));
    }
}
```

边界从 `z = 0` 开始按 `dz` 排列，再由 `[zBottom, gap]` 截断，因此不完整层取
截断后自身的中点。层数直接等于 `tierCenter.size()`。

### `microTPCsample.C`

- 读取同一个 JSONC。
- 从 `events.root/simulation_metadata` 读取 `tierCenter`，不重复计算层数。
- 一次运行生成全部 `tiers/tierN.root`。
- 每个 tier 文件写入该层的 `templateZ`。

### `analysis/mean.C`

- 对每层中心通道的原始波形和成形波形求事例平均。
- 输出到 `<outputDir>/mean/meanN.root`。
- 把输入 tier 文件中的 `templateZ` 原样写入 mean 文件。

### `analysis/matix.C`

- 从 `mean0.root` 开始按连续编号读取，遇到第一个不存在的文件停止。
- 把所有层的平均成形波形写入 `matrix.root/matrix_tree`。
- 把各层 x/y 平均波形叠加图写入 `matrix.root/meanWaveforms`。
- 不输出单独的图片文件。

### `run.sh`

- 配置并编译 simulation。
- 运行 `microTPCsimulation` 和 `microTPCsample`。
- 依次处理连续编号的 tier 文件。
- 生成 mean 文件和最终 `matrix.root`。

## 最新运行结果

已使用当前代码和 `apv25_1.jsonc` 完整重跑成功：

```text
events.root 事件数 = 1000
abs(vz)             = 0.00415333 cm/ns
dz                  = 0.103833 cm
层数                = 5
tierCenter           = 0.0569167, 0.1557500, 0.2595833,
                       0.3634167, 0.4576667 cm

primaryElectrons     = 29925
driftSuccess         = 29856
readoutInRange       = 29856
libraryOutOfRange    = 0
emptyLibraryRegion   = 0
sampledElectrons     = 29856
```

已核对：

- 新 `events.root` 中有 1000 个事件和 1 条 metadata；
- metadata 中已有最新的 `outputDir`、`detectorDir`、`gasFile`、`fieldFile`、
  `gap` 和 5 个 `tierCenter`；
- `tier0.root` 至 `tier4.root` 均已重新生成；
- `mean0.root` 至 `mean4.root` 均已重新生成；
- `matrix_tree` 有 5 条记录；
- `meanWaveforms` 已保存在 `matrix.root` 中。

## 代码风格

继续修改时按下面的习惯写：

1. 代码要能直接看出物理量如何计算，避免没有必要的封装。
2. 推导可以展开写清楚，不为了少几行把公式藏起来。
3. 常见物理量使用常见短名称，例如 `n`、`tau`、`z0`、`dz`。
4. 其他名称使用清楚的 camelCase，例如 `sampleInterval`、`totalTime`、
   `tierCenter`、`outputDir`、`detectorDir`、`gasFile`、`fieldFile`。
5. 单位写在结构体和 JSONC 注释中，不机械地加进变量名。
6. 配置成员能直接使用时就写 `config.simulation.xxx` 或
   `config.electronics.xxx`，不另建同义变量。
7. 局部变量应表示真正的中间物理量，例如 `z0`、`dz`、`lower`、`upper`。
8. 模拟和电子学参数分开，但结构体名称不加多余的 `Parameters`。
9. 使用 JSONC 统一配置并写清注释，不把参数拆成复杂的命令行选项。
10. 新增 metadata 的命名尽量与配置和 C++ 成员一致。
11. 连续编号文件使用简单的 `while (!gSystem->AccessPathName(...))` 方式发现。
12. 注释主要说明物理含义、单位和为什么这样计算。

## 新聊天阅读顺序

1. `simulation/CHAT_HANDOFF.md`
2. `simulation/Config.h`
3. `simulation/config/apv25_1.jsonc`
4. `simulation/microTPCsimulation.C`
5. `simulation/microTPCsample.C`
6. `simulation/analysis/mean.C`
7. `simulation/analysis/matix.C`
8. `simulation/run.sh`
9. `simulation/SIMULATION_WORKFLOW.md`

新聊天可直接发送：

```text
请先读 simulation/CHAT_HANDOFF.md，再看相关实际代码。它记录了当前 simulation
改动、最新验证结果和我的代码风格。请只修改我接下来指定的内容。
```
