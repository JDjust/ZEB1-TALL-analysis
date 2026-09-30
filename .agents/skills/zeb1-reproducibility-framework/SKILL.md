---
name: zeb1-reproducibility-framework
description: Optimize this ZEB1 T-ALL project through traceable data/QC, reproducible figure exports and a single current delivery. Use for changes to this project's datasets, scripts, figures, manuscripts or GitHub release; ground decisions in the attached inspected sources.
---

# ZEB1 可复现分析与科研图形框架

由女娲的主题工作流提炼，用中性专业表达，不模拟任何专家本人。本框架基于 8 个实际获取的一手来源及本项目代码/日志的定向审查。它不是完整人物蒸馏，也不是全项目统计效度认证。

## 执行顺序

1. 读取根 README、相关 panel README 和 data/QC_LIMITATIONS.md；确认是 Full_Master、Blood Advances 还是 Scientific Reports 的稿图编号。
2. 从问题反向找到输入表、处理代码、实际 QC 日志和图注。按具体阶段记录四状态：代码规定（无执行证据）；已有执行证据（链接原日志、输入、参数和当时环境，并说明适用性）；本轮执行（记录命令、输入 SHA256、参数、环境和验收）；未知（列出缺失记录）。绘图、缓存重算、原始数据处理分别标记，不能扩大通过范围。数据再分发权限另外标为已核对、受限或未知。
3. 选择能解决已观察缺陷的改动。保存数据定义与真实测量，允许改变实现、目录、字体、排版和标注。新分析若影响结论，先明确其问题、单位和所需输入。
4. 在固定目录修订；上游重处理与冻结表绘图分开。默认单 panel 仅更新目标 panel 及 receipt；运行前后检查 canonical 整版与其他 panel 未改变。--out 的所有导出写入指定目录，裁切所需整版副产品在入口说明中注明。期刊导出与通用版本隔离。
5. 检查受影响输出与科学单位、尺度、键、统计标注；检查实物尺寸和独立 panel 的 keys。只验证实际改动范围。将结果和未解决项留档。

## 六个决策模型

- **结果是一条依赖链。** Figure → panel → 源表 → 处理脚本 → accession/QC → 环境。Sandve 与 Wilson 支持可追踪结果；本项目用 manifest/索引落实。局限：记录路径不证明方法有效，也不证明有权重新分发输入。
- **干净会话是可运行性的最低证据。** R4DS 的 scripts/projects 要求持久代码与相对项目路径；本项目以 Rscript --vanilla 和共享 Python 引擎验证。局限：冻结表重绘通过不证明远端 raw-to-figure 通过。
- **显示单位先于漂亮样式。** Wickham 的观察单位整洁数据原则与 Wilke 的 faithful-data 讨论共同支持标明患者/供者/基因，不能把细胞当独立 n。局限：好的图不能弥补未经验证的恶性标签或混杂。
- **图的含义必须随图走。** Wilke 的 legend、axis、redundant coding 指导独立 panel 保留必要 keys；Nature skill 提供尺寸/字体/导出 QA。局限：通用 profile 不是 Blood Advances 或 Scientific Reports 的硬要求。
- **冻结与改进并行。** 可修复已知键错误、非收敛成功标记、裁剪和版本混乱，但保留原来源、冻结值和变化说明。局限：改变标签/模型属于重新分析，不能包装成视觉改进。
- **公开包是一项明确选择。** Wilson 的版本/许可规则和用户模板支持可审计代码发布；本项目全本地数据与公开代码分开。局限：文件小、来源公开或去标识都不自动赋予重新分发许可。

## 可操作规则

如果已有原始 data，先复用，不重新从 hulu 下载。缺输入时报告确切文件，不用虚构数据填补。
如果 QC 日志有未收敛/回退，报告状态；“程序返回”不能写成“QC通过”。
如果图中有误差线或 P 值，追溯现有结果/方法，标明单位与校正；不重算凑显著。
如果 standalone panel 有裁剪/混入其它 panel，隔离对象再导出；不要仅裁原图后宣布可复用。
如果生成期刊尺寸版本，输出到明确版别/指定目录；不要覆写 canonical panel。
如果改动产生 new/final2 文件夹，优先固定入口和 Git 历史；保留唯一当前版本。
如果发布 GitHub，检查实际提交清单与来源许可边界，排除患者表、原始/中间矩阵和凭证。

## 方法张力与边界

人可读的逐 panel 目录与少复制共享实现之间，用轻量 R/Python 入口解决。冻结可核验的 baseline 与勇于改进之间，用隔离输出、差异记录与版本控制解决。视觉简化与完整 denominators 之间，以主图直观、补图/表可审计共同解决。

这是一轮有边界的工程/图形主题提炼；没有完整调研人物生平、私有想法或最新观点。上游预处理保留了部分本地/远端依赖，不能保证陌生机器端到端运行。Harmony 未收敛、患者独立性和 DepMap 许可缺口仍需真实证据解决。配色 preset 不能自动提高生物学有效性。

## 按需阅读证据

- 数据、QC 与 clean-session：`references/research/01-reproducibility.md`。
- 独立 panel、keys 与视觉表达：`references/research/02-visualization.md`。
- 发布、版本和运行边界：`references/research/03-release.md`。

> 本主题框架依 [女娲 · Skill造人术](https://github.com/alchaincyf/nuwa-skill) 的方法生成；上游创建者：[花叔](https://x.com/AlchainHust)。本项目提炼与应用：Codex，2026-09-30。

## 缓存、期刊与独立 panel 的验收

改变 QC 阈值、细胞/样本筛选、注释、归一化、批次校正、模型或随机种子时，列出受影响的缓存。只有输入哈希、参数和相关代码版本一致才复用；否则重算相关阶段。缺缓存来源记录标为未知；缺输入则停止依赖该变更的结果更新。旧 H5AD 不能代表新参数已经运行。

期刊尺寸、字体与字母样式先查目标期刊官方指南，记录 URL、访问日期与适用文章类型。Nature skill 的通用 profile 仅作为视觉建议，不能代替 Blood Advances 或 Scientific Reports 的投稿要求。指南尚未核验时保留当前已交付期刊版并标注未重新核验，不自动套用 Nature 170 mm。

独立 panel 用 owner artists 导出；逐项核对颜色、形状、大小和线型的含义，以及 n、P、CI、相关系数与冻结表/图注一致。缺必要 key、统计文字、坐标单位或出现邻图内容即 FAIL；无关共享 keys 属于可读性 WARN，不能把机械 QA 写成语义完整认证。
