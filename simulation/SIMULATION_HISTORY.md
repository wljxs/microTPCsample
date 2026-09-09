# Simulation 历史与验证记录

本文归档重写工作流文档前的详细内容，包括探测器结构、气体输运参数、抽样库约束、
性能测量以及3/4/5层重建结果。这里的路径、命令和“当前配置”可能已经过时；
实际运行方法以 [SIMULATION_WORKFLOW.md](SIMULATION_WORKFLOW.md) 为准。

---

# `microTPCsimulation.C` 模拟流程说明

本文说明当前 `microTPCsimulation.C` 的实际模拟链路。程序结合 Garfield++ 的轨迹和漂移模拟，以及预先建立的单电子波形抽样库，生成逐条波形 ROOT 文件。

## 总体流程

```text
读取单电子波形抽样库
        ↓
配置气体、COMSOL 电场和读出时间窗
        ↓
TrackHeed 生成 20° muon 轨迹与初级电离电子
        ↓
AvalancheMC 将每个初级电子漂移到读出面
        ↓
按漂移终点在抽样库中随机选择单电子响应
        ↓
将响应按到达时间平移并累加到 x/y strips
        ↓
用成形核卷积，写入 ROOT 的 tree/tree2/tree3
```

## 目录和文件职责

从用途上看，`simulation/` 中的内容可以分为四类：主模拟程序、单电子抽样库制作、
探测器输入数据，以及编译/模拟产生的结果。当前完整链路是：

```text
microTPC.C
  生成带位置标签的单电子雪崩及感应信号
        ↓
analysis/merge.C
  合并多批单电子 ROOT 文件
        ↓
analysis/buildSampleLibrary.C
  按局部 (x0, y0) 分区建立单电子抽样库
        ↓
analysis/samplelibrarysmall.root
        ↓
microTPCsimulation.C
  生成粒子径迹、漂移初级电子、抽取并叠加模板
        ↓
result/*.root
```

### 主模拟程序

- `microTPCsimulation.C`：当前完整事例模拟的主程序。它用 `TrackHeed` 生成
  150 GeV muon 的初级电离，用 `AvalancheMC` 将电子漂移到读出面，按终点的
  strip 和局部位置从 `samplelibrarysmall.root` 随机抽取单电子响应，最后进行时间
  平移、逐 strip 累加和 50 ns 成形卷积。它写出波形、初级电离真值、漂移终点真值
  和运行配置元数据。
- `microTPC.C`：单电子响应的原始生成器。它在读出面附近发射单电子，使用
  `AvalancheMicroscopic` 显式模拟高场放大区的雪崩，并利用 COMSOL weighting
  potential 计算中心及相邻 x/y strips 的感应信号。输出文件是建立抽样库的原始输入。
- `microTPCsample.C`：按漂移深度层生成直接模板的辅助程序。它在指定 tier 的中心
  高度放置一个电子，模拟漂移后抽取单电子响应，输出 10 条 x、10 条 y strip 波形及
  漂移终点信息。输出 strip 坐标以最终落点为中心，`sigx2/sigy2` 直接保存落点所在
  strip 的响应。`build/tiersimulation_*` 和后续 `mean*.root` 属于这条模板制作支线，
  而不是完整 muon 事例模拟。
- `gas.C`：独立的 Magboltz 气体表生成示例。当前代码生成的是
  Xe/Ne/iC4H10 并写出 `trd_xe.gas`，与主流程加载的 Ar/CF4/CO2 气体不同，不能把它
  当作当前主模拟所用 `.gas` 文件的生成配置。

### 抽样库和模板分析程序

- `analysis/merge.C`：查找 `build/run/result/` 下的 ROOT 文件，按文件名中的数字自然
  排序，合并其中的 `tree`、`tree_after_conv`、`tree2` 和 `tree3`，输出
  `merged1.root`。
- `analysis/buildSampleLibrary.C`：读取 `merged1.root` 中逐条对齐的 `tree` 和
  `tree3`，只接受 `z0 = 0.01 cm` 附近的单电子样本，将局部 `(x0, y0)` 划成
  `10 x 10` 个区域，输出 `samplelibrarysmall.root`。库中的 `meta` 保存分箱范围和
  时间轴，`region_ix_iy` 保存对应区域的 8 路单电子波形。
- `analysis/mean.C`：读取一个 `tierN.root`，对已经按最终落点居中的
  `sigx2/sigy2` 求平均，生成只有一个平均模板 entry 的 `meanN.root`。
- `analysis/matix.C`：读取连续编号的 `mean0.root`、`mean1.root` 等文件，将不同
  深度层的 x/y 平均波形组合成 `matrix.root`，供分析端的反卷积或 NNLS 重建使用。
  文件名拼作 `matix.C`，其中的 ROOT 宏函数名实际是 `matrix()`。

### 探测器和气体输入

- `microTPC5mmsmall/`：生成单电子响应时使用的多 strip COMSOL 模型。
  `mesh.mphtxt` 是有限元网格，`mplist.txt` 是 domain/材料映射，`field.txt` 是实际
  电势场；`xlabel*.txt` 和 `ylabel*.txt` 是各读出 strip 的 weighting potential。
- `microTPC5mmone/`：完整径迹漂移使用的单周期实际电场模型。此阶段不重新计算
  weighting signal，而是从抽样库取响应，因此该目录不需要 `xlabel/ylabel` 文件。
- `ar_45_cf4_40_co2_15100-100000.gas`：Ar 45%、CF4 40%、CO2 15% 的 Magboltz
  输运表。
- `field250.txt`：另一套电压配置的场文件；当前源码读取的是 `field.txt`，没有使用它。
- 两个模型目录的 `Readme.md`：简要记录 5 mm 气隙、DLC +500 V、drift -500 V
  和 U-groove 几何信息。

### 构建文件和数据产物

- `CMakeLists.txt`：查找 Garfield++、ROOT、Eigen，生成 ROOT 字典，并构建
  `microTPC`、`microTPCsimulation`、`microTPCsample` 和 `gas`。
- `Electron.h`：定义可由 ROOT 序列化的 `fCluster`、`fElectron` 数据类。
- `ElectronLinkDef.h`：声明这些类及其 `std::vector` 的 ROOT 字典。当前主要程序
  没有直接使用这两个类，属于可复用基础设施或历史设计遗留。
- `build/`：CMake 缓存、目标文件、ROOT 字典、可执行文件，以及部分 tier 模拟结果；
  其中程序二进制不一定与当前源码同步。
- `.cache/clangd/`：编辑器生成的代码索引，不参与模拟。
- `analysis/samplelibrarysmall.root`：按位置分区的单电子响应库，是主模拟的重要输入。
- `analysis/3direct*`、`analysis/4direct`、`analysis/5*`：不同层数或生成方式的平均
  模板与矩阵，主要供重建比较。
- `result/angle_scan/`：不同角度、随机起点模式下的完整模拟结果。
- `result/nodelta/`：禁用 delta-electron transport 的完整模拟结果。

### 当前源码需要注意的不一致

1. `analysis/mean.C` 默认读取 `../build/tiersimulation/tierN.root`，但仓库现有目录多为
   `tiersimulation_3direct`、`tiersimulation_4direct` 等；使用时需要匹配正确目录。
2. `CMakeLists.txt` 中的 `Eigen3_DIR` 是机器相关的绝对路径，换环境时应改成实际安装
   位置或让 CMake 从标准搜索路径寻找 Eigen。
3. `gas.C` 的混合气与当前主模拟不同，修改或重新生成气体表时应先确认目标气体配置。

## 坐标、单位和轨迹真值

Garfield 内部长度单位为 **cm**。`pitch = 0.04` cm，即 0.4 mm；气隙 z 范围约为 `0.01--0.51` cm，即 0.1--5.1 mm。分析端若使用 mm，应将坐标乘以 10。

每个事件从随机位置开始：

```cpp
x0 = 2 * pitch + Uniform(0, 2 * pitch);
y0 = 2 * pitch + Uniform(0, 2 * pitch);
z0 = 0.49;  // cm
t0 = 0.;
```

方向定义为：

```cpp
dx = 0.;
dy = sin(20 deg);
dz = -cos(20 deg);
```

因此轨迹在 **y-z 平面**倾斜，`kx = 0`、`ky = -tan(20 deg)`。真值采用 `x = kx * z + bx`、`y = ky * z + by`；斜率无量纲，截距的单位为 cm。程序禁用了 delta-electron transport。

## 模拟的探测器结构

这是一个带 **U-groove 读出结构**的微型气体 TPC / MicroMegas 类几何。`microTPC5mmone/Readme.md` 将场图描述为“气隙 5 mm、DLC 500、drift -500”的 U-groove 模型。COMSOL 工程名为 `microTPCunit3.mph`，网格和静电势均以 mm 为单位导出。

从 `field.txt` 的节点坐标与电势可恢复出下列主要结构；具体材料名称不能仅由导出文件完全确定，因此 DLC、mesh 和基板的命名以模型说明为准：

```text
z ≈ 5.06 mm : 漂移阴极，约 -500 V
      │
      │ 约 5 mm 漂移气隙，平均漂移场量级约 1 kV/cm
      │
z ≈ 0.05 mm : 约 0 V 的 mesh / 电极平面
      │
      │ 约 50 μm 高场间隙，0 V 到 +500 V，量级约 100 kV/cm
      │
z ≈ 0 mm    : 约 +500 V 的 DLC / 读出侧
      │
z < 0 mm    : U-groove 介质、金属与读出基板区域，场图延伸到约 -0.2 mm
```

COMSOL 网格包含 108146 个顶点，静电势文件包含 108147 个节点；电势范围为 -500 到 +500 V。`mplist.txt` 定义了四类材料：相对介电常数分别为 `1`、`1e10`、`4.5`、`3.5`。其中 `1e10` 在 Garfield 场图中表示导体，`epsilon = 1` 的 domain 是气体漂移介质；其余两类是不同介电常数的绝缘层。场图的 domain 编号到材料索引的映射也保存在 `mplist.txt`。

当前程序的 `Sensor` 区域从 `z = 0.01 cm = 0.1 mm` 开始。因此 TrackHeed 电子在约 5 mm 漂移区中被追踪至读出侧附近；程序并未在这份轨迹代码中调用 Garfield 的 avalanche 过程来显式模拟高场间隙内的雪崩倍增。读出侧的单电子响应由抽样库提供，而非由本 event 循环重新计算。

## 1. 单电子波形抽样库

程序打开 `../analysis/samplelibrarysmall.root`。该文件应包含：

- `meta`：`xmin/xmax/ymin/ymax/nbins/t`，描述局部位置分箱和时间轴；
- `region_ix_iy`：每个局部位置 bin 的单电子波形样本。

每个 region 包含 8 路信号：x 方向的 `sigxtestminus1`、`sigxtest1`、`sigxtestminus2`、`sigzero`、`sigxtest2`，以及 y 方向的 `sigytestminus1`、`sigy`、`sigytest1`。

每一个漂移成功的电子会在对应 `(ix, iy)` region 中随机选择一条样本。这样保留单电子响应随位置变化和样本涨落，同时不必对每个电子重新计算感应信号。

### 抽样库的来源与约束

抽样库的原始输入由 `simulation/microTPC.C` 产生，随后可用
`simulation/analysis/merge.C` 合并，再由
`simulation/analysis/buildSampleLibrary.C` 分区构建。`microTPC.C` 对每个
event 写入两棵逐条对应的树：

- `tree`：8 路单电子波形；
- `tree3`：该波形对应的初始条件 `x0, y0, z0, t0, e0`。

当前 `microTPC.C` 固定 `z0 = 0.01 cm`，并在
`x0 = 0.056--0.064 cm`、`y0 = 0.04--0.08 cm` 的范围内随机取点。建库程序
按 `x0 × y0` 划分为 `10 × 10` 个局部区域，并显式要求
`abs(z0 - 0.01 cm) <= 1e-6 cm`；同时检查 `tree` 与 `tree3` 的条目数一致。
这保证同一 region 中的模板来自相同初始高度，且波形与位置标签不会错位。

旧版建库代码的注释声称按 `z0` 筛选，但实际没有执行该检查。若输入确实完全
由当前 `microTPC.C` 产生，则所有条目的 `z0` 本来就是 `0.01 cm`，旧库内容
不会因此混入不同高度；但该逻辑不适用于未来混合高度的输入。建议用现版建库程序
重新生成抽样库，以保存 `z0` 分支并得到筛选统计。

当前库每个区域约有 400 条单电子模板。它足以描述常见响应的平均形状，但并不能
充分约束稀有的高增益或异常波形尾部：概率为 1% 的响应在 400 条中平均仅有 4 条，
概率为 0.1% 的响应通常不会被采到。若要报告增益尾部、阈值效率或高分位数，应增加
每区独立模板数，并用 100/200/400/800 条子库比较结果是否收敛。

## 2. 气体、电场和读出时间

气体由 `microTPC5mmone/ar_45_cf4_40_co2_15100-100000.gas` 加载，并启用 Ar Penning transfer。COMSOL 电场读取 `mesh.mphtxt`、`mplist.txt` 和 `field.txt`，x/y 方向设置为周期性边界。

气体文件是 Magboltz 生成的输运表，条件为：

- Ar 45%、CF4 40%、CO2 15%；
- 293.15 K、1 atm、无磁场；
- 30 个约化电场点，范围为 `0.1316--131.6 V/(cm Torr)`；在 1 atm 下约等价于 `100--100000 V/cm`，覆盖低场漂移和高场放大区；
- 程序额外启用 Ar 的 Penning transfer，概率 `rPenning = 0.20`。

因此，该混合气的漂移速度、扩散和电离/激发相关输运参数均随局部 COMSOL 电场插值取得；而不是在代码中使用一个固定漂移速度。

### 从当前气体表读出的代表性输运参数

下表使用 Garfield++ 直接加载当前 `.gas` 文件，并按程序相同设置启用 `rPenning = 0.20` 后查询。无磁场时，速度沿电场方向；`D_L` 和 `D_T` 的单位是 `μm / √cm`，Townsend 系数 `alpha` 与附着系数 `eta` 的单位是 `cm^-1`。

| 电场 (V/cm) | 漂移速度 (mm/ns) | `D_L` | `D_T` | `alpha` | `eta` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 0.00427 | 231.4 | 228.6 | 0 | 0 |
| 500 | 0.02061 | 124.9 | 126.0 | 0 | 0 |
| 1000 | 0.03922 | 101.8 | 100.3 | 0 | 0 |
| 5000 | 0.10108 | 79.4 | 164.6 | 0 | 0.364 |
| 10000 | 0.08934 | 97.7 | 208.6 | 0.003 | 16.0 |
| 50000 | 0.15870 | 121.7 | 155.8 | 513.1 | 31.4 |
| 100000 | 0.27507 | 102.8 | 119.1 | 1762.3 | 14.1 |

漂移区由约 -500 V 到 0 V 跨越约 5 mm，对应的平均场约为 `1 kV/cm`。所以本模型中最有代表性的低场参数是：

- 漂移速度约 `0.039 mm/ns`；电子穿越 5 mm 约需 `5 / 0.03922 ≈ 128 ns`；
- 对 5 mm（0.5 cm）漂移距离，扩散 RMS 的量级为 `D × √L`，即纵向约 `102 × √0.5 ≈ 72 μm`、横向约 `100 × √0.5 ≈ 71 μm`；
- 因此分析配置中的 `vdrift = 0.04 mm/ns` 与当前气体表和漂移场是一致的近似。

高场间隙的电场量级可达 `100 kV/cm`。表中可见此时 `alpha` 已很大，物理上应出现强雪崩倍增；而在漂移区 `alpha = 0`，初级电子只发生漂移和扩散。需要注意：当前 `microTPCsimulation.C` 调用的是 `AvalancheMC::DriftElectron`，没有在 event 循环中显式调用雪崩过程。因此高场增益不由该循环直接产生，若输出波形具有增益，其幅度和增益涨落应由单电子抽样库的生成过程承担。

CF4 使该混合气在 kV/cm 量级具有较快漂移速度；CO2 提供淬灭和输运调节；Ar 是主要电离组分。Penning transfer 设置会提高 Ar 激发态向可电离通道转化的概率，尤其会影响高场区的有效 Townsend 系数和预期增益。该物理效果应与生成抽样库时的设置保持一致。

读出时间窗为 140 个采样点：

```text
tmin  = 0 ns
tmax  = 700 ns
tstep = 5 ns
```

程序配置单极性一阶 shaper，时间常数为 50 ns。`kernel_vec` 直接由同一个
Garfield `Shaper::Shape(t)` 生成，用于后续的显式离散卷积；因此不会出现手写
传递函数与 Garfield shaper 参数不一致的问题。

## 3. TrackHeed 初级电离和电子漂移

每个 event 调用：

```cpp
track.NewTrack(x0, y0, z0, t0, dx, dy, dz);
```

随后用 `GetCluster` 遍历 Heed 的电离 cluster，用 `GetElectron` 取得每个初级电子的 `(xe0, ye0, ze0, te0, ee0)`。

每个电子调用：

```cpp
drift.DriftElectron(xe0, ye0, ze0, te0, ee0);
```

漂移成功后，终点 `(xe1, ye1, ze1, te1)` 写入 `tree2`；初始电离位置写入 `tree3`。因此 `tree3` 是电离真值，`tree2` 是漂移/扩散后的读出面终点；波形的 strip 分配使用终点，而不是初始位置。

## 4. 单电子响应的时间平移和 strip 累加

终点先映射到 x/y strip，再根据一个 pitch 内的局部位置计算抽样库 bin `(ix, iy)`。x 方向按位置选用二条或三条相邻 strip 的模板；y 方向累加三条相邻 strip 的模板。

到达时间 `te1` 决定库波形的起始位置：

```cpp
fracBegin = te1 / tstep;
beginindex = int(fracBegin);
alpha = fracBegin - beginindex;
```

每个库波形使用线性插值：

```text
previous = (j > 0) ? sample[j - 1] : 0
value = (1 - alpha) * sample[j] + alpha * previous
```

这里模板起点之前视为零，并使用前一个样本实现 `alpha` 个采样点的延迟。使用
`sample[j + 1]` 会反向把模板提前。随后将响应加入电子影响到的每条 x/y strip 的
未卷积波形。一个 event 的全部初级电子重复此过程，所以输出包含电离涨落、漂移
终点分布和抽样库响应涨落。

## 5. 成形卷积和输出

所有电子累加完成后，每条 strip 的未卷积响应与 `kernel_vec` 做离散卷积，得到 `vec_sigx_conv` / `vec_sigy_conv`。只有卷积后非零的 strip 写入输出。

当前输出文件为：

```text
../result/nodelta/20test.root
```

相对路径以程序运行目录（通常是 `simulation/build`）为准。

程序现在支持可复现的角度扫描命令行：

```text
./microTPCsimulation \
  --run-id N \
  --theta-deg V \
  --phi-deg V \
  --events N \
  --seed N \
  --position-mode center|uniform \
  --output PATH
```

`theta` 是相对向下z方向的极角，`phi` 从+x向+y计数；`phi=90 deg` 与原来的
Y-Z平面20度轨迹一致。旧的 `[runid] [output-root-file] [n-events]` 位置参数形式
仍兼容。随机种子必须非零，并同时设置本地抽样和Garfield随机引擎。

角度扫描不再把顶部起点固定在第2--4条附近，而是先指定轨迹在
`z=0.25 cm` 处穿过20条读出的中心，再反算顶部起点。`center` 固定在中心；
`uniform` 在中心附近一个pitch内均匀抽样，用于平均局部U-groove相位。输出新增
`simulation_metadata`，记录角度、方位角、seed、事件数、位置模式和电子接受统计；
`tree` 新增顶部起点及中点坐标。

运行结束会打印从初级电子、漂移成功电子、读出范围内电子到抽样库成功取样电子的计数；其中 `library (ix, iy) out of range` 是检查抽样区越界的直接指标。

以当前 20° 轨迹、当前随机起点范围、`microTPC5mmone` 场图和 `samplelibrarysmall.root` 做的 100-event 临时运行中：2884 个初级电子中 2881 个漂移成功且位于读出范围，`(ix, iy)` 越界数为 0，空 region 数也为 0。这说明在当前配置下，抽样库分箱覆盖了实际到达读出区域的电子；改变入射位置、角度、气体或场图后应重新检查该计数。

### 运行速度

在相同配置下的 100-event 性能测量中，事件循环约需 19.7 s（整段进程约 21.7 s），
约有 2800 个初级电子。分段计时显示 `AvalancheMC::DriftElectron` 约占 75%，
ROOT 模板随机读取与波形叠加约占 23%，140 点显式卷积不足 1%。因此当前速度主要
受有限元场中的逐电子数值漂移限制，而不是卷积；模板库很大时应使用有界的区域缓存
或批量读取，而不宜简单地把整个库常驻内存。

## ROOT 输出结构

### `tree`：供波形分析读取

| 分支 | 含义 |
| --- | --- |
| `event` | event 编号，包含 `runid * nEvents` 偏移 |
| `t` | 140 点时间轴，单位 ns |
| `types` | strip 类型：0 为 x，1 为 y |
| `stripIDs` | 对应方向内的 strip 编号，范围 0--19 |
| `waveforms` | 未卷积波形 |
| `waveforms_conv` | 成形卷积后的波形 |
| `kx`, `bx`, `ky`, `by` | 轨迹真值，坐标单位为 cm |

同一索引的 `types`、`stripIDs`、`waveforms` 和 `waveforms_conv` 对应同一条 strip。

### `tree2`：漂移终点真值

保存每个漂移成功、确实到达 `z = 0.01 cm` 读出面且位于读出范围内的电子终点位置、
到达时间以及 x/y strip 编号。Garfield 的 endpoint 接口不提供终点能量，因此不再写
语义不明确的 `e1` 分支。

`microTPCsample` 产生的 tier 文件带有 `signal_coordinates=endpoint_centered` 和
`tier_schema_version=2` 元数据。已有旧 tier 文件可用
`analysis/convertLegacyTier.C` 按 `tree2` 中的最终落点转换；转换会去掉旧
`template_sig*` 和 `e1` 分支。

### `tree3`：初始电离真值

保存每个 TrackHeed 初级电子的初始位置、产生时间、能量和初始位置对应的 x/y strip 编号。

## 与分析程序的衔接

`analysis/src/simulationanalysis.C` 读取 `tree.waveforms_conv`，每 5 个 5 ns 模拟点取一个点，形成 15 个、间隔 25 ns 的分析采样点，再进行解卷积和 microTPC 重建。分析端会检查模板/波形长度，并且仅使用通过 cluster quality cuts 的 `isValid` cluster 进入位置、角度和相关性统计。

模拟和分析目前都使用 `simulation/result/nodelta/20test.root`，可直接进行闭环分析。

### 原始 5 层、未做后期位置修正的结论

当前程序固定使用 5 个直接独立模板，读取
`simulation/analysis/5/matrix.root`。所有下列数值来自
`analysis/result/20test_nodelta_midestbig.root` 的 5000-event 原始重建：没有逐层位置
偏置扣除、没有角度零点扣除，也没有用结果样本反过来标定。

| 指标 | X | Y |
| --- | ---: | ---: |
| microTPC 位置残差均值 / RMS | -0.13 / 43.2 μm | +10.0 / 32.8 μm |
| 角残差均值 / RMS | +0.00003 / 0.0204 | +0.0139 / 0.0284 |

五层的逐层位置残差在 X 向均值接近零、RMS 为约 `73--140 μm`；Y 向中间三层均值约为
`+22、+2、-2 μm`，两端层约为 `-28、+59 μm`，RMS 为约 `78--101 μm`。因此原始
Y 角残差的正偏置主要来自端层相反的位置偏置在直线拟合中的传播；不能把经后期位置或
角度零点修正的数值当成原始算法分辨率。

下表列出所有层的**未修正**统计。`q16/q84` 是中心 68% 区间的边界；它比 RMS 更能
反映厚尾对“典型事件”分辨率的影响。单位均为 μm。

| 层中心 z (mm) | X: 均值 | X: RMS | X: q16 / q84 | Y: 均值 | Y: RMS | Y: q16 / q84 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 | +1.1 | 73.1 | -56.7 / +58.2 | -27.7 | 87.5 | -102.3 / +47.4 |
| 1.5 | -0.4 | 106.7 | -70.2 / +68.6 | +22.1 | 77.9 | -50.8 / +94.9 |
| 2.5 | +0.6 | 140.5 | -105.7 / +105.3 | +1.6 | 77.8 | -69.9 / +72.9 |
| 3.5 | +5.4 | 139.0 | -112.6 / +131.3 | -1.6 | 84.9 | -78.7 / +79.1 |
| 4.5 | -1.2 | 132.0 | -144.9 / +142.3 | +59.0 | 101.2 | -24.8 / +137.5 |

Y 的第 0 和第 4 层分别有负、正的长尾（偏度约 `-1.28`、`+2.00`），这也是原始 Y
角度正偏置的直接来源。X 的第 0--3 层主要表现为近零均值但 RMS 大于中心 68% 宽度，
说明低电荷事件形成了对称或近对称厚尾，而不是整体位置平移。

部分层的残差不是严格高斯。X 向低 tier-charge 样本有明显厚尾，来自高相关模板下 NNLS
在相邻层之间的电荷互相补偿；Y 向端层还叠加了模板/连续深度响应与离散 5 层基函数之间
的模型失配。五层模板矩阵的条件数约为 X/Y = `92/97`，说明这些小的波形或模板误差会被
放大。正式报告原始算法时，应同时给出各层 RMS、均值、中心 68% 区间和电荷切片，而不应
仅用单高斯拟合掩盖厚尾。

### 3 层与 4 层的原始对照结果（未修正）

下列 ROOT 文件是此前产生的历史对照，保留用于比较；当前 `SimulationAnalysis` 代码并
**不**运行这些配置，也没有对这些数值施加任何位置或角度后期校正。

| 层数与模板 | 原始结果文件 | X microTPC：均值 / RMS (μm) | Y microTPC：均值 / RMS (μm) |
| --- | --- | ---: | ---: |
| 3 层，5 层模板深度插值 | `20test_nodelta_midestbig_3tier.root` | -0.10 / 44.4 | +11.1 / 50.5 |
| 4 层，直接模板 | `20test_nodelta_midestbig_4tier_direct_diagnostics.root` | -0.06 / 44.2 | +7.3 / 37.7 |
| 5 层，当前原始配置 | `20test_nodelta_midestbig.root` | -0.13 / 43.2 | +10.0 / 32.8 |

3 层：

| 层中心 z (mm) | X：均值 / RMS (μm) | Y：均值 / RMS (μm) |
| ---: | ---: | ---: |
| 0.833 | -0.5 / 56.4 | -97.1 / 101.9 |
| 2.500 | +0.9 / 76.3 | +11.3 / 85.9 |
| 4.167 | -0.2 / 90.1 | +122.6 / 105.1 |

4 层：

| 层中心 z (mm) | X：均值 / RMS (μm) | Y：均值 / RMS (μm) |
| ---: | ---: | ---: |
| 0.625 | +0.4 / 60.5 | -58.8 / 92.6 |
| 1.875 | -0.6 / 86.2 | +2.8 / 79.6 |
| 3.125 | +3.0 / 105.1 | +5.1 / 81.6 |
| 4.375 | -0.3 / 111.6 | +84.8 / 99.2 |

这些逐层表保留不同层数的各自层中心，不能逐行当作相同深度的直接比较。它说明：
3 层的端层 Y 偏置最大；4 层减轻了该问题，但仍有端层偏置；5 层在 Y 的中间层最好，
代价是模板相关性更高、X 的分层 RMS 更宽。

### 原始分层残差协方差矩阵

以下矩阵来自输出 ROOT 文件的 `correlation/correlation_one_{x,y}/h_cov_type_{0,1}`；行和列
均按从低 z 到高 z 的 tier 编号排列，单位为 `mm²`。非零非对角元表示层间残差共同涨落或
电荷互相补偿，负值尤其表示一个层偏正时另一个层倾向偏负。

#### 3 层

```text
Cov_X = [[ 0.003181,  0.000909,  0.001816],
         [ 0.000909,  0.005822, -0.000685],
         [ 0.001816, -0.000685,  0.008120]]

Cov_Y = [[ 0.010373, -0.001531,  0.000723],
         [-0.001531,  0.007384, -0.001295],
         [ 0.000723, -0.001295,  0.011233]]
```

#### 4 层

```text
Cov_X = [[ 0.003659,  0.001272,  0.001221,  0.000925],
         [ 0.001272,  0.007433, -0.000258,  0.002817],
         [ 0.001221, -0.000258,  0.011038, -0.003297],
         [ 0.000925,  0.002817, -0.003297,  0.012445]]

Cov_Y = [[ 0.008975, -0.001140,  0.000319,  0.000285],
         [-0.001140,  0.006335, -0.001407,  0.000396],
         [ 0.000319, -0.001407,  0.006664, -0.000507],
         [ 0.000285,  0.000396, -0.000507,  0.011047]]
```

#### 5 层（当前原始配置）

```text
Cov_X = [[ 0.005339,  0.001334,  0.001321,  0.000883,  0.000833],
         [ 0.001334,  0.011393, -0.000369, -0.000042,  0.002555],
         [ 0.001321, -0.000369,  0.019738, -0.003118,  0.002944],
         [ 0.000883, -0.000042, -0.003118,  0.019332, -0.008246],
         [ 0.000833,  0.002555,  0.002944, -0.008246,  0.017436]]

Cov_Y = [[ 0.007944, -0.000590,  0.000449, -0.000201,  0.000198],
         [-0.000590,  0.006063, -0.000958,  0.000614,  0.000093],
         [ 0.000449, -0.000958,  0.006060, -0.000936,  0.000314],
         [-0.000201,  0.000614, -0.000936,  0.007211, -0.000194],
         [ 0.000198,  0.000093,  0.000314, -0.000194,  0.015887]]
```

5 层 X 的第 3/4 层（从 0 开始编号）协方差为 `-0.008246 mm²`，是最明显的层间反补偿
项；这与高条件数、该两层模板难以区分的现象一致。

## 运行前检查和注意事项

1. `samplelibrarysmall.root` 存在，且 `meta.nbins` 必须与固定数组大小 `kNbins = 10` 一致。
2. 抽样库的气体、电场、几何、时间窗和成形函数应与本程序一致。
3. `microTPC5mmone` 下的 COMSOL 与 gas 文件应齐全；输出目录也应已存在。
4. 相邻 strip 写入在边界 strip 附近需要数组范围保护。
5. `(ix, iy)` 映射越界的电子会被跳过，建议记录其比例。
6. 显式卷积使用包含 `j = k` 的标准因果求和；反卷积模板也应采用同一边界定义。
7. 建议后续在 ROOT metadata 中记录几何、gas、时间窗、shaper、随机种子和抽样库版本，保证结果可复现。
