# Simulation 工作流说明

本文说明 `simulation/` 目录中每个主要文件的职责，以及当前从配置到模板矩阵的完整流程。

旧文档中的探测器结构、气体输运参数、性能测量和历史重建结果已移到
[`SIMULATION_HISTORY.md`](SIMULATION_HISTORY.md) 保留。

## 一键运行

在 `simulation/` 目录执行：

```bash
./run.sh
```

当前脚本使用 `config/apv25_1.jsonc`。切换电子学方案时，只需修改
`run.sh` 开头的 `CONFIG_FILE`。

```text
config/apv25_1.jsonc
        ↓
microTPCsimulation.C
        ↓
result/apv25_1/events.root
        ↓ 读取 simulation_metadata.tierCenter
microTPCsample.C
        ↓
result/apv25_1/tiers/tierN.root
        ↓
analysis/mean.C
        ↓
result/apv25_1/mean/meanN.root
        ↓
analysis/matix.C
        ↓
result/apv25_1/matrix.root
```

## 每个文件在做什么

### 配置与运行

#### `Config.h`

定义 `Simulation`、`Electronics` 和 `Config` 三个结构，并提供公共的
`ReadConfig()`。`microTPCsimulation.C` 和 `microTPCsample.C` 使用同一套结构和命名。
解析器支持带 `//` 注释的 JSONC；文件不存在时使用结构体默认值。

#### `config/default.jsonc`

默认配置和字段说明。程序没有收到配置路径时读取该文件。

#### `config/apv25_1.jsonc`

当前 APV25 第1组配置。主模拟和每一层模板各运行1000个事例，输出目录为
`result/apv25_1/`。

#### `SIMULATION_HISTORY.md`

保留旧文档中的详细物理说明和历史验证结果。其中的旧路径和旧运行方式仅供追溯，
不作为当前操作指南。

#### `run.sh`

一键入口，依次完成：

1. CMake配置和编译；
2. 运行 `microTPCsimulation`；
3. 运行 `microTPCsample`；
4. 对每个tier调用 `mean.C`；
5. 调用 `matrix()` 生成矩阵。

脚本只控制执行顺序，物理参数全部来自JSONC。

### 模拟程序

#### `microTPCsimulation.C`

完整muon径迹模拟主程序：

1. 读取JSONC和单电子抽样库；
2. 按配置加载探测器目录、气体文件和COMSOL电场；
3. 用 `TrackHeed` 产生径迹和初级电离电子；
4. 用 `AvalancheMC` 将电子漂移到读出面；
5. 按终点位置从抽样库选择单电子响应；
6. 按到达时间平移并累加各strip信号；
7. 用配置中的 `n` 和 `tau` 做成形卷积；
8. 输出 `events.root`。

程序还会在漂移区中心取得电场，用 Garfield 的 `ElectronVelocity()` 计算z方向
漂移速度，并由电子学采样周期计算各层中心位置。

`deltaElectron` 控制delta电子输运：

- `true`：`EnableDeltaElectronTransport()`；
- `false`：`DisableDeltaElectronTransport()`。

#### `microTPCsample.C`

生成各深度层的模板事例。它读取同一JSONC，再从
`events.root/simulation_metadata` 读取 `tierCenter`。中心数组有几个元素，就依次
生成几个 `tiers/tierN.root`，不再需要手工传层号和层数。

每个tier运行 `config.simulation.nEvents` 个事例，`templateZ` 保存该层中心位置。

#### `microTPC.C`

生成单电子响应抽样库的原始数据。它使用 `microTPC5mmsmall/` 的实际场和
weighting potentials，通过 `AvalancheMicroscopic` 模拟放大区雪崩，输出8路单电子
感应波形和初始位置。它属于抽样库制作支线，不在日常一键流程中运行。

#### `microTPCsimulationxray.C`

X-ray版本的完整事例模拟。该文件是独立旧入口，目前未接入 `Config.h` 和 `run.sh`。

#### `gas.C`

Magboltz气体表生成示例。当前生成Xe/Ne/iC4H10的 `trd_xe.gas`，并不是APV25配置
使用的Ar/CF4/CO2气体。

### 分析宏

#### `analysis/merge.C`

从 `avalanche0.root` 开始查找连续编号的单电子文件，将 `tree`、
`tree_after_conv`、`tree2` 和 `tree3` 合并为 `merged1.root`。

#### `analysis/buildSampleLibrary.C`

读取 `merged1.root`，按单电子初始局部位置 `(x0, y0)` 分成 `10 × 10` 个区域，
生成 `analysis/samplelibrarysmall.root`。其中：

- `meta` 保存位置范围、分箱数和时间轴；
- `region_ix_iy` 保存对应区域的8路单电子波形。

主模拟和分层模板都依赖这个抽样库。

#### `analysis/mean.C`

读取一个 `tierN.root`，分别计算x/y原始波形和卷积波形的事例平均，输出：

```text
<outputDir>/mean/meanN.root
```

它复制 `templateZ`，但不计算层数或层中心。

#### `analysis/matix.C`

文件名是 `matix.C`，宏函数名是 `matrix()`。它读取连续编号的 `meanN.root`，
将各层平均卷积波形写入 `matrix.root`。ROOT文件包含：

- `matrix_tree`：各层中心位置和x/y模板；
- `meanWaveforms`：各层x/y平均波形的叠加画布。

画布只保存在ROOT文件内，不另外生成PNG。

### 构建和ROOT字典

#### `CMakeLists.txt`

查找Garfield++、ROOT和Eigen，生成ROOT字典并编译 `microTPC`、
`microTPCsimulation`、`microTPCsample`、`microTPCsimulationxray` 和 `gas`。

#### `Electron.h`

定义可由ROOT序列化的 `fCluster` 和 `fElectron`。

#### `ElectronLinkDef.h`

声明上述类型及其容器的ROOT字典。

### 探测器目录

#### `microTPC5mmone/`

主模拟和分层模板使用的单周期COMSOL模型。主要文件包括：

- `mesh.mphtxt`：有限元网格；
- `mplist.txt`：材料/domain映射；
- `field.txt`、`field250.txt`：不同电场文件；
- `*.gas`：Garfield/Magboltz气体输运表。

具体目录、气体和场文件由JSONC中的 `detectorDir`、`gasFile`、`fieldFile` 指定。

#### `microTPC5mmsmall/`

制作单电子响应时使用的多strip模型，还包含各读出条的weighting-potential文件。

#### `microTPC5mmone/Readme.md` 和 `microTPC5mmsmall/Readme.md`

记录相应COMSOL模型的几何、电压和气隙信息。

## JSONC配置字段

```jsonc
{
    "outputDir": "../result/apv25_1",
    "totalTime": 700.0,             // ns

    "simulation": {
        "runid": 0,
        "nEvents": 1000,
        "detectorDir": "../microTPC5mmone",
        "gasFile": "ar_45_cf4_40_co2_15100-100000.gas",
        "fieldFile": "field.txt",
        "theta": 20.0,              // degree
        "phi": 90.0,                // degree
        "seed": 42,
        "mode": "center",           // "center" 或 "uniform"
        "deltaElectron": false,
        "sampleInterval": 5.0,      // ns
        "gap": 0.5                  // 漂移区上边界z0，cm
    },

    "electronics": {
        "n": 1,
        "tau": 50.0,                // ns
        "sampleInterval": 25.0      // ns
    }
}
```

两个 `sampleInterval` 的用途不同：

- `simulation.sampleInterval`：模拟波形时间间隔和卷积时间步长；
- `electronics.sampleInterval`：电子学采样周期，只用于计算深度层。

当前代码不会按电子学采样周期重新抽样输出波形。

## 分层方法

主模拟先计算：

```text
dz = abs(vz) × electronics.sampleInterval
```

然后以 `z = 0` 为起点按 `dz` 分层，再用 `zBottom` 和 `gap` 截断：

```cpp
lower = max(zBottom, tier * dz);
upper = min((tier + 1) * dz, gap);
center = 0.5 * (lower + upper);
```

`zBottom = 0.01 cm`。最下层可能不完整，因此使用截断后该层自身的中点。
主模拟把最终的 `gap` 和 `tierCenter` 写入metadata，sample不重复计算。

## ROOT输出

### `events.root`

- `tree`：逐事例时间轴、strip编号、卷积前/后波形和径迹真值；
- `tree2`：漂移成功电子的终点、时间和strip编号；
- `tree3`：TrackHeed初级电子的位置、时间、能量和strip编号；
- `simulation_metadata`：实际配置、`gap`、`tierCenter` 和运行计数。

metadata主要分支：

```text
theta  phi  seed  nEvents  runid  mode  deltaElectron
detectorDir  gasFile  fieldFile
simulationSampleInterval  n  tau  electronicsSampleInterval
totalTime  outputDir  gap  tierCenter
nStrips  pitch  sampleLibrary
primaryElectrons  driftSuccess  readoutInRange
libraryOutOfRange  emptyLibraryRegion  sampledElectrons
```

### `tierN.root`

- `tree`：每个模板事例的10条x和10条y波形；
- `tree2`：电子漂移终点；
- `templateZ`：层中心，单位cm；
- `signalCoordinates`、`tierSchemaVersion`：格式说明。

### `meanN.root`

- `tree_mean`：一条平均波形；
- `mean_sigx`、`mean_sigy`：平均原始波形图；
- `mean_sigx_conv`、`mean_sigy_conv`：平均卷积波形图；
- `templateZ`：层中心。

### `matrix.root`

- `matrix_tree`：各层模板矩阵；
- `meanWaveforms`：所有层平均卷积波形的叠加画布。

## 输出目录

```text
result/apv25_1/
├── events.root
├── tiers/
│   ├── tier0.root
│   └── ...
├── mean/
│   ├── mean0.root
│   └── ...
└── matrix.root
```

## 注意事项

1. 程序通常从 `simulation/build/` 运行，因此JSONC相对路径也以该目录为基准。
2. `samplelibrarysmall.root` 的气体、场、几何和时间格点应与当前配置兼容。
3. `simulation.sampleInterval` 当前应与抽样库时间间隔一致。
4. `electronics.sampleInterval` 只控制深度分层，不改变波形时间采样。
5. `microTPCsample` 必须在主模拟之后运行，因为层中心来自 `events.root`。
6. tier和mean文件必须从0开始连续编号。
7. 更换气体、探测器或电场后，应重新生成events、tier、mean和matrix。
