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

程序配置单极性一阶 shaper，时间常数为 50 ns；同时按相同参数生成 `kernel_vec`，用于后续的显式离散卷积。

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
value = (1 - alpha) * sample[j] + alpha * sample[j + 1]
```

并加入电子影响到的每条 x/y strip 的未卷积波形。一个 event 的全部初级电子重复此过程，所以输出包含电离涨落、漂移终点分布和抽样库响应涨落。

## 5. 成形卷积和输出

所有电子累加完成后，每条 strip 的未卷积响应与 `kernel_vec` 做离散卷积，得到 `vec_sigx_conv` / `vec_sigy_conv`。只有卷积后非零的 strip 写入输出。

当前输出文件为：

```text
../result/nodelta/20test.root
```

相对路径以程序运行目录（通常是 `simulation/build`）为准。

程序支持可选命令行参数：

```text
./microTPCsimulation [runid] [output-root-file] [n-events]
```

省略后两个参数时，仍使用默认输出路径和默认 event 数。运行结束会打印从初级电子、漂移成功电子、读出范围内电子到抽样库成功取样电子的计数；其中 `library (ix, iy) out of range` 是检查抽样区越界的直接指标。

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

保存每个漂移成功且位于读出范围内的电子终点位置、到达时间、能量以及 x/y strip 编号。

### `tree3`：初始电离真值

保存每个 TrackHeed 初级电子的初始位置、产生时间、能量和初始位置对应的 x/y strip 编号。

## 与分析程序的衔接

`analysis/src/simulationanalysis.C` 读取 `tree.waveforms_conv`，每 5 个 5 ns 模拟点取一个点，形成 15 个、间隔 25 ns 的分析采样点，再进行解卷积和 microTPC 重建。分析端会检查模板/波形长度，并且仅使用通过 cluster quality cuts 的 `isValid` cluster 进入位置、角度和相关性统计。

模拟和分析目前都使用 `simulation/result/nodelta/20test.root`，可直接进行闭环分析。

## 运行前检查和注意事项

1. `samplelibrarysmall.root` 存在，且 `meta.nbins` 必须与固定数组大小 `kNbins = 10` 一致。
2. 抽样库的气体、电场、几何、时间窗和成形函数应与本程序一致。
3. `microTPC5mmone` 下的 COMSOL 与 gas 文件应齐全；输出目录也应已存在。
4. 相邻 strip 写入在边界 strip 附近需要数组范围保护。
5. `(ix, iy)` 映射越界的电子会被跳过，建议记录其比例。
6. 显式卷积使用包含 `j = k` 的标准因果求和；反卷积模板也应采用同一边界定义。
7. 建议后续在 ROOT metadata 中记录几何、gas、时间窗、shaper、随机种子和抽样库版本，保证结果可复现。
