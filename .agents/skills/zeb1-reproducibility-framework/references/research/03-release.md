# ZEB1 发布与计算复现框架：证据及项目审计

调研日期：2026-09-30。范围：读取暂存的 `clean_project` 和用户的 GitHub 结构模板；没有修改科学分析代码、数据、图或原 GitHub checkout。下述项目建议是本次审计推论，不冒充作者原话。

## 已抓取的一手来源

1. Wilson 等，**Good enough practices in scientific computing**，PLOS Computational Biology，2017，DOI：**10.1371/journal.pcbi.1005510**。页面标题已核对。原文：https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510 。已存 `../sources/wilson-2017.html` 及文字提取。原文建议逐步整理可分析数据，记录中间步骤，提供人类可读 README，明确许可证及引用方式，使用版本控制。它支持适度的工程规范，不要求每个探索项目都建复杂平台。
2. Sandve 等，**Ten Simple Rules for Reproducible Computational Research**，PLOS Computational Biology，2013，DOI：**10.1371/journal.pcbi.1003285**。页面标题已核对。原文：https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003285 。已存 `../sources/sandve-2013.html` 及文字提取。实际读取规则 1–10：追踪每个结果如何产生，避免手工数据修改，记录程序版本和脚本版本，保存中间结果与随机种子，保存图后数据，关联文字结论与底层结果，公开可共享的代码和结果。
3. **About large files on GitHub**，GitHub 官方文档：https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github 。已存 `../sources/github-large-files.html`。实际读取：超过 50 MiB 有警告，超过 100 MiB 普通 Git 推送被阻止；建议仓库理想小于 1 GB，强烈建议小于 5 GB。上述是 GitHub 技术限制，不能用于判断数据是否适合公开。

抓取说明：`npx --no-install firecrawl-cli --status` 已检查，CLI 存在但未认证。因此本次按已知权威 URL 用标准库 HTTP 抓取；没有虚称使用 Firecrawl 研究索引或完成额外文献检索。共抓取上述 3 个来源。

## 可蒸馏的主题方法

**模型 1：每张图是一条可追踪依赖链。** 从稿件的 panel 反向定位脚本、读入表、预处理程序、原始 accession 与环境。依据 Sandve 规则 1、7、9；Wilson 的 README 和项目组织建议。具体应用：panel 注释与 manifest 写明输入相对路径、列含义、统计单位、图层含义、脚本版本和已验证的运行范围。

**模型 2：把当前工作副本与公开发布包分开。** 本地保留现有数据和复现资料；公开包仅收纳明确允许共享的代码、汇总结果与元数据。依据 Wilson 的明确许可证与版本控制建议、GitHub 技术文档，并结合本项目模板明确的 raw/processed/objects 不提交要求。这是项目工程推论，不是上述作者关于本项目数据许可的判断。

**模型 3：可复现是分层事实，不是一句宣传。** A 层：冻结表重绘；B 层：现有中间对象重新计算；C 层：原始 accession 重新处理。分别记录已运行、未运行、缺输入和需服务器步骤。依据 Sandve 关于完整工作流、程序版本和中间结果的规则。仅 A 层通过时，不能标注 C 层已验证。

**模型 4：使一次有效变化留下一个清楚版本。** 根目录和用户常用文件夹只有当前有效稿件/图；探索记录放 Git 历史或外部可回滚归档。依据 Wilson 的 small changes 和版本控制建议。这并不授权丢掉仍能解释方法选择的独特分析结果。

## 现场确认的风险及可实施改进

| 风险 | 现场证据 | 可实施改进与验收 |
|---|---|---|
| 把完整本地数据直接推上 GitHub | `data/analysis/total/rebuild_2026/data/deepen_2026/a7_lineage/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz` 1481.88 MiB；`a4_yayon/thymus_scrna_tcell_subset.h5ad` 1388.37 MiB；`a7b_malignant/GSE248287_10x_RNA_matrix.mtx.gz` 821.96 MiB；`data/analysis/data/TALL_X01_vst.tsv` 530.20 MiB | 用独立 public_release 目录和明确允许清单构建公开包；禁止整个 data/analysis 递归复制。推送前检查 Git 实际 tracked/staged 文件，所有文件小于 100 MiB，且非批准矩阵、对象、压缩包均排除。保留本地文件，不再从 hulu 下载。 |
| 小文件仍包含个体/临床字段 | 实际路径包括 `data/analysis/data/polonen_syn54032669/audit_clinical_1309.tsv`，`data/analysis/total/data/clinical_mrd_persistence_20260923/paired_patients_predictions.tsv`，`.../validation/clinical_endpoint_audit/protocol_model_predictions.tsv` | 文件大小和扩展名只能辅助排除，不能代替许可/字段检查。公开 source-data 默认仅审定聚合表；个体级图数据保留本地并说明受限输入路径与 accession，未核许可不放 public release。患者编号即使去标识也不自动意味着可自由再发布。 |
| 迁移后科学预处理脚本未保证可运行 | 多个 `data/analysis/total/rebuild_2026/code/deepen_2026/*.py` 硬编码 `D:/_bioinformation/ZEB1/data/analysis...`；`a4_parse_cellxgene.py` 和 `a4_inspect_cellxgene_meta.py` 依赖旧 `C:\Users\Ls180\.cursor\...\agent-tools\...txt` | 保留科学计算步骤并标注 historical/unvalidated；按必要工作流逐步替换为项目根目录与明确输入参数。缺失旧工具缓存的脚本不能列入自动 pipeline 或承诺端到端复现。 |
| journal 导出污染 canonical panel | `script/_engine/render_figure.py` 的 `_render()` 中 plate output 使用 `DEST`，但 paneldir 固定写入 `PROJECT/script/{Figure,SFigure}/...` | journal 导出到独立输出目录，并让 panel 输出跟随导出范围；或 journal 模式不生成 canonical panel 文件。验收：journal 导出前后 canonical panel 校验值不变。 |
| 单 panel 运行有额外覆盖 | 同一 `_render()` 在 `SELECTED` 过滤 panel 前保存整张 PDF/PNG/SVG | 可保留整版重建机制，但 README 明确写出副作用；更好是 panel-only 临时绘制整版，再仅替换对应 panel 文件。与用户约定调用方式一致后验收。 |
| 输入记录未完整 | `render_figure.py` 只包装 `pd.read_csv`，记录路径无哈希；没有覆盖所有文件读取 API | 实際读入表形成 manifest，记录 path、SHA256、脚本、环境及 seed；不把 CSV 跟踪描述为完整数据 lineage。外部路径用安全标签/明确处理，避免 `relative_to(PROJECT)` 直接抛错。 |
| 字节一致与科学一致混为一谈 | receipt 有 `rendered_at`；PDF 可含创建时间；渲染随机 jitter 多处已有显式 seed | 验收数字/点数/轴范围、内容及渲染图，而非要求所有重绘 PDF 的 SHA 必然相同。manifest 区分输入哈希、科学指标和仅输出存档哈希。 |
| 公开仓库夹杂运行缓存 | staged `script/_engine/__pycache__` 有 21 个 `.pyc` | `.gitignore` 排除 `__pycache__/`、`*.pyc`、`.partial`、作业日志、工具缓存；发布构建仅复制源码。 |

## 建议的最小公开结构

遵循用户模板的功能规范，不再复制制造第二套稿件或整套图。

```text
public_release/
  README.md                  # 当前研究问题、验证范围、输入和图的映射
  LICENSE                    # 只对明确自有代码声明授权
  CITATION.cff               # 已知作者和仓库；无 DOI 则不捏造
  .gitignore
  requirements.lock.txt      # 实际运行环境的版本记录，按事实命名
  environment.yml            # 可用环境定义；有经过验证的环境才作保证
  run_all.sh                 # 仅运行实际验证过的绘图入口
  scripts/06_figures/         # 当前 plotting engine 与 panel 入口
  scripts/01_preprocess/      # 仍服务当前结果的预处理源程序及 QC 文档
  src/                       # 共享路径/风格/复现工具
  data/README.md             # accession、原始来源、许可和输入获取说明
  data/metadata/             # 审定非个体标识元数据及冻结基因定义
  results/tables/            # 审定聚合表
  results/figures/           # 当前最终图；必要格式，不复制所有预览
  docs/                      # QC、指标定义、verified/not verified、环境
```

R 入口默认通过环境变量选择 Python（如 `ZEB1_PYTHON`），本机默认值可以是已验证的 D 盘解释器；公开脚本不能假定其他人也有 D 盘。本地 Rscript 路径同理提供配置，同时记录实际执行程序及版本。R 包未重新运行或未安装时，不虚构 `renv.lock`。Docker/Snakemake 是可选升级，不应为目录整洁而生成未运行的复杂骨架。

## 诚实边界与流派张力

1. 保存所有中间结果有利追溯，但会增加磁盘占用；公开 Git 适合小型当前源码及审定表，本地/服务器适合大对象。此处不是要求删除大数据。
2. 固定版本帮助重现，但过度锁定旧环境会妨碍维护。旧分析环境与当前绘图环境分别记录；不能用当前 Python 版本冒充旧分析时版本。
3. 最简目录便于使用，但过度删除方法探索可能使筛选/敏感性决定无法解释。删除重复展示版本，不删除独特科学证据。
4. 本次未检查所有表的字段、数据条款或整个结果链；没有证明本项目能从原始数据一键重跑。
5. 来源给出方法原则，不能证明现有 ZEB1 生物学结论成立；仍须独立统计和证据审计。

## 对新主题 skill 的可执行规则

收到“整理/优化/发布/复现”时，先定位当前稿件、panel 清单、source-data 清单和成功运行记录；若缺记录，先进行冻结表重绘。发布前构造实际文件允许清单，检查大小、路径、临床字段和许可；未获审定个体表保持本地。每次输出报告列出 verified（冻结表重绘）、retained unvalidated（科学预处理代码）、external dependency（本地数据或服务器对象）。不得把目录重排称为全流程验证，也不得创建假版本、假作者或假引用。


## 本轮修订状态

此前审计发现的 --out 写入其他位置、默认单 panel 覆盖 canonical 整版已在 render_figure.py 修复。验收见 docs/OUTPUT_ISOLATION.json。原表为审计快照，上游完整流程的路径与环境限制仍然有效。
