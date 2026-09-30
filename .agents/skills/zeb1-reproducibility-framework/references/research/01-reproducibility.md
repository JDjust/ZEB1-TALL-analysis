# ZEB1 可复现工作流视角：来源与项目实查

审查日期：2026-09-30。依据女娲的主题框架工作流，提取可操作方法，不模拟 Hadley Wickham 本人、不代替其发表意见。本轮只读取已有代码、数据表、日志并核验两个哈希；未重新执行统计、单细胞处理或服务器任务，未修改科学数据。

## 1. 实际取得的一手方法来源

1. Hadley Wickham、Mine Çetinkaya-Rundel、Garrett Grolemund，*R for Data Science (2e)*，第 6 章 **Workflow: scripts and projects**。实际 HTTP 取得并阅读：<https://r4ds.hadley.nz/workflow-scripts.html>。缓存：`references/sources/articles/r4ds-workflow-scripts.html`，提取正文：同名 `.txt`。
   - §6.1.1：脚本声明依赖；共享的分析脚本不直接执行 `install.packages()`。
   - §6.1.3：名称可读、可排序；不同类型文件分别放入目录。
   - §6.2.1：代码和数据文件共同构成持久记录；重启 R 后重跑，暴露隐含环境依赖。
   - §6.2.3–6.2.4：一个分析对应一个项目；脚本产生图；项目内部使用相对路径。
2. Hadley Wickham，**Tidy data**，作者官方文章页：<https://vita.had.co.nz/papers/tidy-data.html>。实际 HTTP 取得并阅读页面摘要与书目信息；没有阅读论文 PDF 全文。缓存：`references/sources/articles/wickham-tidy-data.html`，提取正文：同名 `.txt`。
   - 页面明确给出：每个变量一列、每个观察一行、每类观察单位单独一张表。适用于整理数据集登记、样本登记、QC 表及 panel 对应表。

检索失败记录：`https://r4ds.hadley.nz/workflow-projects.html` 返回 HTTP 404，不作为来源；有关项目的内容已在上面的第 6 章核实。Firecrawl 在父任务中已确认未认证，因此本轮通过 `Invoke-WebRequest` 直接取得已知官方页面。

上述来源支持组织与复现原则，不支持下文生物学结论，也不提供本项目 QC 阈值的合理性证明。

## 2. 本地证据与具体 QC 状态

以下路径均相对于整理后项目根目录；均已实际阅读。代码存在、输出存在、执行日志存在是三种不同证据，不互相替代。

### GSE227122：执行证据充分，但整合收敛存在明确限制

- 代码：`data/analysis/modules/module5/code/01_prepare_module5.py:154–267`。
- 上游实际日志：`data/analysis/total/data/validation/scrna_patient_qc/01_prepare_module5.log`。
- 登记：`data/analysis/total/data/validation/scrna_patient_qc/analysis_provenance.json`；本轮算得日志 SHA256 与登记一致。
- 已执行规则：先计算线粒体/QC 指标；保存过滤前逐细胞 QC；细胞至少表达 200 个基因；基因至少在 20 个细胞中检测；保留 `n_genes_by_counts < 8000`、`pct_counts_mt < 20`；按 sample 运行 Scrublet 并排除预测 doublet。
- 日志明确记录原始合并对象 **41,316 × 33,538**；过滤后 **31,297 × 17,675**（日志第 37 行）；Scrublet 标记 **121** 个 doublet（第 312 行）。过滤后 31,297 是 doublet 排除前数量，不能称为所有 QC 完成后的最终细胞数。
- 代码将过滤后原始 counts 存入 `layers['counts']`，之后总量归一化至 10,000、`log1p`；HVG/PCA/邻居/Harmony/Leiden 是下游表示与标签步骤。QC 阈值顺序应按代码记录，避免用一句“标准流程”代替。
- **明确未通过的条件**：日志第 343 行为 `Harmony stopped before convergence after 10 iterations`；第 358 行仍为 `omicverse harmony OK`，其原因是代码只检查返回的坐标存在。`OK` 不等于收敛。第 367 行有 igraph Leiden 完成记录。
- 代码遇 Scrublet 异常会填充 `predicted_doublet=False`，遇 HVG/Harmony/Leiden 异常会选择回退；本次日志证明 Scrublet/Leiden 正常执行，不能因此保证未来运行都会采用同一算法。
- 代码在已存在 `GSE227122_annotated.h5ad` 时直接读缓存，不核查输入/参数/代码指纹。未来参数改变可能复用旧处理对象；本轮没有运行该路径。

### GSE227122 / GSE248287：患者级深度与细胞数敏感性确有输出

- 代码：`data/analysis/total/code/scrna_patient_qc.py`、`scrna_patient_qc_sensitivity.py`。
- 实际任务输出：`data/analysis/total/data/validation/scrna_patient_qc/patient_qc_1730.out`；README 报告 Hulu job 1730 成功，本轮未重新查询 Slurm。
- 输出：同目录 `patient_cell_depth_detection.tsv`、`cell_threshold_sensitivity.tsv`；本轮确认后者共 **10** 行，即两队列 × 五阈值。
- 患者单位为 GSE227122 **10** 人、GSE248287 **15** 人；代码做 one-to-one 合并，核对恶性细胞数、summed UMI、ZEB1 counts 与冻结 pseudobulk 表一致。
- 已执行敏感性：最低恶性细胞数 **30 / 100 / 250 / 500 / 1000**；保持全队列原先基因尺度，保存所有排除 ID；每个患者集合固定 seed；10,000 次置换、2,000 次患者 bootstrap。重复阈值保留相同患者时得到相同估计，不能称五次独立验证。
- GSE227122 恶性细胞来自项目启发式标签，GSE248287 来自作者标签；pseudobulk 使用未整合 counts 不能独立消除 GSE227122 聚类标签对 Harmony 的依赖。
- 现有输出中两队列相关方向相反、区间较宽。整理或美化不能将其改成方向一致的验证。
- **尚未建立**：这些 retained diagnosis-cell 统计不重建每位患者从原始输入到所有 QC 排除后的完整流程数量；GSE248287 排序/混合样本中的细胞比例也不是原位肿瘤负荷。

### TARGET：STAR 汇总代理已执行，不是读段质控

- 代码：`data/analysis/total/code/target_star_raw_qc.py`。
- 输出与登记：`data/analysis/total/data/validation/target_star_raw_sample_qc.tsv`、`target_star_raw_qc_manifest.json`。
- 本轮读表确认 **265** 行；现场计算 TSV SHA256 与 manifest 的 `output_sha256` 一致。
- 已实施代码约束：265 个唯一 USI；合法 basename；原文件存在；STAR gene-model header；`unstranded` 整数计数；四个 `N_*` 汇总类别齐全且 assigned counts 非零。输出检测基因数、线粒体/蛋白编码比例、assigned fraction 等代理。
- 这些是 STAR count 文件汇总指标；未读取 FASTQ/BAM。本轮没有发现或执行 FastQC/read trimming/比对重跑证据，不能声称完成读段级 mapping QC。
- 生存核查文档：`data/analysis/total/data/validation/target_survival/README.md`。该记录保留 `None` 的原始事件类别，区分真正缺失，262 位 OS/EFS 纳入者及 14/27 个事件；不是新的独立验证队列。此轮未独立重算 Cox。

### GSE144035：比较定义修正和计数检查已记录

- 代码：`data/analysis/total/code/gse144035_paired_audit.py`。
- 文件：`data/analysis/total/data/validation/gse144035/README.md`、`paired_model_changes.tsv`；本轮确认表中 **12** 行（6 对模型 × 2 基因）。
- 已有审计：14 个 GSM，8 个 VTL 模型；只有 6 个模型具有 Plus/Minus Dox 对应；严格取整数计数列、检查非负整数，以供应计数表全部行计算 library sum 和 `log2(1+CPM)`。
- 作者表 `Independent -Dox` 被误作 OFF/ON 的旧比较已撤回；原行保留供审计。不能将作者汇总列视为同模型干预效应；不能由 6 对模型推断独立生物学重复数。
- 该程序尾部会修改 Figure 6E source table；复现“审计”命令可能有写入副作用，应明确写入范围或将清洗与核查拆开。此轮没有执行程序。

### GSE287751：作者 QC 条件和项目执行应分开标明

- 文件：`data/analysis/total/data/validation/gse287751/README.md`、`metadata_provenance.json`。
- 原始作者 README 被保存在 provenance：作者用 HTO 排除 multiplets/unassigned；过滤 `<2000` 基因或线粒体比例 `>6%` 的细胞。这是作者描述，不是本轮独立重做的 QC。
- 项目记录：11,661 个保留细胞的 barcode 全匹配原矩阵；八个作者条件到 HTO 唯一映射；按条件对 raw UMI 求和，再做 `log2(1+CPM)`；README 记录 job 1715 成功，此轮没有查询 Slurm。
- 四种 context 是不同条件，不能视作同一比较的独立重复；现有表只给描述性方向，不生成细胞伪重复 P 值。

### 版本与独立性：现有登记已承认未解决项

- `data/analysis/total/data/dataset_registry.tsv` 提供 accession、来源 URL、父子 series、缓存 hash 和访问限制；本轮只查阅前部记录，不声称逐行验证整个登记。
- `data/analysis/total/data/validation/cohort_overlap/README.md` 明确只在兼容 ID namespace 比较；无法映射的队列独立性是未解决，不能写为零重叠。
- `data/analysis/total/data/validation/depmap_release/README.md` 明确四输入已锁定哈希，但 Public26Q1 官方版本、发布日期和许可尚未独立核实；HTTP200 challenge 不是 release 列表。不要以相关 Figshare 的日期/许可补齐。

## 3. 三项最优先可实施改进

1. **给整理后的 data 建立唯一索引与 QC 说明。** 建 `data/DATA_PROCESSING.md` 和 tidy 的 `data/dataset_index.tsv`：每行记录 accession/观察单位/本地处理表/处理脚本/原始服务器路径/具体 QC/执行证据/已知限制/版本状态；另用 panel 索引指向同一输入。保留服务器 raw locator，默认不触发下载。把上述具体阈值、作者已做与项目已做、Harmony 未收敛、DepMap 未核实并列写清楚。验收：从任一最终 panel 能找到输入、处理代码和真实 QC 证据，不存在第二份“最终数据”。
2. **每个 panel 的 R 入口从干净会话运行并记录来源。** 路径相对项目根，D 盘 R/Python 路径由配置或环境变量声明；共享 Python 绘图实现只有一份；入口声明 panel、输入表、单位/变换、输出位置，并执行依赖预检。重绘后记录输入 SHA256、执行命令、实际 R/Python/包版本与输出 SHA256；保留绘图与上游处理两个明确层级。验收：`Rscript --vanilla` 在任意工作目录可重绘选定 panel，与固定 source table 一致；不依赖 `.RData` 或临时缓存，不静默下载/安装。
3. **未来重处理将“有结果”与“QC 成功”分开。** 将 Harmony 返回状态/收敛警告、Scrublet/HVG/聚类回退、缓存 fingerprint 明确写入处理 manifest，必要条件失败时停止或显式标记降级。加入逐阶段每样本保留/剔除数量。不能先覆盖本次冻结结果再解释差异；新算法/标签须作为敏感性分支并完整报告。验收：带有未收敛警告的试运行不会标成通过；输入/参数变化不能无提示复用旧 H5AD；现有稿件继续保留原日志揭示的限制。

## 4. 方法张力与适用边界

- 人可读目录满足定位需求；唯一源表、共享实现和 manifest 满足代码可维护需求。panel 目录可以包含 R 入口、README 与本 panel 图，避免复制整套数据或整个共享绘图模块。
- 有输入哈希的冻结图版可以被审计；全 raw-to-figure 重跑需要 Hulu 原始数据与环境。因此当前应分别称“本地最终图复现”和“有上游处理代码/执行证据”，不声称所有 raw-to-figure 在本机已重新执行。
- 勇于改进目录、入口、校验和图形规范；生物学结论的增强必须靠新证据。更漂亮的图和更整齐的目录不消除混杂、低事件数、未收敛或独立性缺口。
- 本轮两个官方方法页面足以支持这些组织原则，不构成完整人物蒸馏、对整个项目所有数据集的完整审计或统计有效性证明。
