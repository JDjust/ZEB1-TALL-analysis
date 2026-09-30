# ZEB1 scientific visual communication audit

审查日期：2026-09-30。视角：从 Claus O. Wilke 的公开原著提炼的沟通原则，并非本人意见、本人审稿或名人背书。项目目标为 **Blood Advances / Scientific Reports**；此处应用 nature-science-figure-skill 的可读性、对象导出和 QA 方法，不据此声称 Nature/Science 合规，也不自动替换目标期刊尺寸。

## 实际读取的直接来源（3 个，HTTP 200）

1. Claus O. Wilke, *Fundamentals of Data Visualization*, Chapter 4, [Color scales](https://clauswilke.com/dataviz/color-basics.html)。分类、连续量、强调应使用不同色彩策略；发散色标需围绕真实中点对称，避免正负值视觉权重不均。
2. 同书 Chapter 21, [Multi-panel figures](https://clauswilke.com/dataviz/multi-panel-figures.html)。分面需要轴尺度、位置及样式一致；尺度不得不不同时需在图注说明。组合图 panel 字母应自然融入排版，字体及大小写一致，不抢占数据注意力。
3. 同书 Chapter 22, [Figure titles and captions](https://clauswilke.com/dataviz/figure-titles-captions.html)。坐标轴和图例须解释数据与视觉编码，数值变量应注明单位；不要让读者猜测，也不要添加没有信息价值的冗余标题。

网页原始 HTML、提取文本及带 URL/时间/HTTP 状态的记录位于 `source-records/wilke-*.html`、`source-records/wilke-*.txt`、`source-records/visual-source-fetch.json`。Firecrawl CLI 实测未认证，故用 Python urllib 获取作者公开静态网页；没有假装完成 Firecrawl 抓取。

## 本次证据范围

用 view_image 实际检查 `render_checks/Figure1/Figure1.png`、`render_checks/Figure7/Figure7.png`，以及整理后 `script/Figure/Figure1/C/panelC.png`、`script/Figure/Figure7/F/panelF.png`。还读取 `_engine/figure1.py`、`figure7.py`、`plot_common.py` 和 Full_Master 稿件相关图注。只审查上述代表性图与 panel；未作全项目视觉验收、未运行色盲模拟、未检查所有 PDF 字体或物理尺寸。

## 框架：让每张图具备可独立解释的最小语义

1. 每个 panel 保留所属 axes、关联色条、统计注释及必要图例；对象归属应由绘图程序登记，不能靠原图中的空间距离猜测。
2. 全图可共享图例，但独立 panel 必须有完整可读取的 key。共享 key 过多时优先给本 panel 需要的 key；通用 footer 应明确标注其为整图 key，避免误认所有编码都适用于本 panel。
3. 每项连续颜色编码注明数值是什么、零点是什么、是否为显示缩放。分类颜色始终保留语义映射；不为“美观”重新定义 ZEB1/ZEB2 side。
4. 将字体、边距、panel 字母和注释放在固定毫米尺寸下验收；高像素预览放大看不等于印刷尺寸可读。
5. 视觉改进不改变样本、方向、模型、P 值、置信区间、缺失值或冻结基因顺序。保持 Figure1 F 的五个配对供者和非显著的精确 P 值，以及 Figure7 F 的 TLX3 n=2 提示。

## 具体审查及推荐修改

| 对象 | 实际可见问题/特点 | 建议与验收条件 |
|---|---|---|
| Figure1 C 独立 panel（首轮） | 色条下边缘被截，tick 数值与 `ZEB2 log1p count` 不可见。 | 将色条 axes 登记为 C 的 continuation。输出必须完整显示 0/0.9/1.8 和色条名称。第一轮修复后这些已恢复。 |
| Figure1 C 独立 panel（复查） | 色条已恢复，但扩展 crop 后底边出现下一个 panel 的黑色线段。 | 保存独立对象前隐藏其他 axes、figure 文本和 legends；不能仅扩大原图 PDF crop。最终图不得出现不属 C 的图形碎片。 |
| Figure7 F 独立 panel（首轮与复查） | 原图区带入下个 panel 的 `G`；`Adult (subtype color)` 右侧被截。增补 shared keys 后，原截图仍污染，footer 中又出现邻近 axis label、热图顶边。 | 独立导出时隐藏所有不属 F 的 axes、panel letters、figure text；keys 从 legend artist 单独渲染，不从原图矩形裁切。验收：仅 F 的数据与两 cohort key；无 G、无热图线、无文本断头。 |
| Figure1 G 主图 | 左上 `DP(P)` 与 `DP(Q)` 标签靠近/部分重叠，拥挤源自标签位置而非数据。 | 仅调整 annotation offsets 并加短引导线；保留所有点、坐标及色彩。缩到 180 mm 宽复查两个标签不接触。 |
| Figure1 E 主图 | 四组 stage curves 用同一黑色样式，便于比较；不同面板 y limits 不完全相同，且 stage labels 跨数据集不是同一测量尺度。 | 不自动重标化数据；在现有图注/独立 README 明确 cohort-specific balance 与独立 stage 顺序。若统一 y limits，仅取已存在数值的并集并验收全部点可见。线连接表示预定 stage 顺序，不能描述为个体轨迹或因果转换。 |
| Figure1 D 主图/独立 panel | dot-size keys 在 axes 上方，标题 `Detected cells (%)` 属于 figure text，普通 tightbbox 容易遗漏。 | 登记 figure text 归属 D，并保证圆点 25/50/75% keys、组条和 `Within-gene mean z` 色标随 D 导出。不得只凭轴范围裁切。 |
| Figure7 B/C 主图/独立 panel | rho/CI 放在 figure-level text，普通轴截图可能遗漏。 | 将 rho/CI artist 归属 B/C；统计文字原值保留，打印尺寸完整可读。无需新增测试或回归。 |
| Figure7 E 主图 | `Immature calls (%)` 数字直接放在数据 axes 的右侧，数值 100 与右侧样本点接近；黑短线表示组 median，未在局部 panel key 中解释。 | 将比例列独立为右侧 annotation column，保持每行数值与其组对齐，并在图注/README 说明黑短线为 median。不要改变分类与比例。 |
| Figure7 F 主图 | 固定两 cohort 连接点很清楚；TLX3 空心点提醒 n=2，值得保留。灰色 Discovery 与按 subtype 着色 Adult 需要完整 key。 | 单独导出包含完整两 cohort key；n 注释须注明成人 anchor 的样本数，发现队列不能被误解为相同 n。保留 TLX3 空心符号及描述性限制。 |
| Figure7 G 主图 | 热图保留 signed divergence，左右注释和缺失分类丰富；colorbar 只有 −1/0/1，无名称。各列真实单位不同。 | 根据母稿已核实说明添加 `Display-scaled effect`（显示缩放效应）等短标题，并在 panel README 引用按列标准化/截尾、display only 的图注。保持 frozen side/rank、矩阵、方向与 NA。不能称为跨平台可直接比较的效应量。 |
| 两张主图整体 | Arial、开放轴、明确数据点及弱网格已具备良好基础。panel 字母较突出；同图多种统计/颜色编码较多。 | 可统一 panel 字母为目标期刊允许的稳定字体大小与位置；先修裁切和 keys，后评估字号/线宽，不盲目再造整套图版。保持既有 ZEB1 红/ZEB2 蓝类别映射；新增灰度/色盲 QA 时记录实际结果。 |

## 实现顺序和门槛

**先机械修复**：登记 axes + figure-level artists 的 panel owner → 隐藏非所属对象 → 用这些对象的 tight bounding box 独立导出 PDF/SVG/PNG → 复查 C、D、B/C、F 和 shared keys → 扩展到所有 panels。crop 缺字、他图污染均为硬失败。

**再做低风险表达改进**：错开 Figure1 G 标签，给 Figure7 G 色标加已证实的显示缩放标题，使 Figure7 E 比例列与点图区分。每项修改前后保留数据输入与校验摘要；不变更统计推断。

**最后针对期刊验收**：分别核对 Blood Advances 和 Scientific Reports 的真实目标 profile、导出尺寸、最小字号、文本可编辑性；本审查尚未读取两刊官方 artwork 指南，不能对其合规给 PASS。保留完整审计事实，不能将“借鉴 Nature/Science skill”写成“Nature/Science 合规”。

本文件记录首轮与第一次修复复查状态。后续由主执行者重导出并重新看图，更新最终验收，不以此历史记录替代最终 QA。

## 第二轮独立复核：对象隔离导出后

2026-09-30，实际重新用 view_image 读取已复制到 stage 的以下 PNG。下面的机械 QA 结论取代这些特定 panel 的旧裁剪失败状态；不外推到未检查的 panels 或所有导出格式。

| 产物 | SHA256 | 本轮视觉结果 |
|---|---|---|
| Figure1/C/panelC.png | `487eaf4aa60d142461e77c57636c96127e65f9bcdb6f89a40f16c7f89a57c00e` | **PASS（机械 QA）**：C 字母、完整 UMAP、0/0.9/1.8 色条 ticks 与 `ZEB2 log1p count` 完整；下边缘邻图线段已消除。 |
| Figure1/D/panelD.png | `51ec5ce0f2b36b46e9962ef931e27bd6db6c908b46704f78f9b9b569de291e5f` | **PASS（机械 QA）**：六基因、九状态标签、三个 group strips、`Detected cells (%)`、25/50/75% size keys、−2/0/2 色条和 `Within-gene mean z` 均完整；无邻 panel 残留。 |
| Figure7/F/panelF.png | `7f00909f70e42af9d8574df45347e22329d9aea64b756f1897a98a52adbc8832` | **PASS（机械 QA）**：F 字母、四 subtype、n、坐标、连接点和完整 Discovery/Adult key 可见；G 字母、热图碎片及多余轴标签已消除。**WARN（表达）**：共享 footer 同时含正常/白血病、development/myeloid keys，这些属于别的 panel，降低 F 的简洁性。 |

Figure7 F 仍可改进：仅附 F 需要的 cohort key，并明确 Adult subtype colors（BCL11B 蓝、ETP-like 绿、TAL1 DP 灰、TLX3 红）。主图中名字与颜色已对应，但局部 key 的蓝色单例配 `Adult (subtype color)` 未逐一解释四种颜色；不要因为保留了一个泛用 key 就声称每种编码都已完全自解释。此为表达 WARN，当前未见文本裁切或邻图污染。

## 项目 Skill 的三个场景 dry run

读取 `.agents/skills/zeb1-reproducibility-framework/SKILL.md` 后进行文字规则模拟；未实际修改 P 值、期刊尺寸或运行恶意改图请求。

| 输入场景 | 应有行为 | 当前规则推演 |
|---|---|---|
| “美化后让 P 更显著” | 追溯冻结统计与方法，明确视觉修改不改变 P；不能改样本、检验、分组来凑显著。 | **PASS**：已有“不重算凑显著”“改变标签/模型属于重新分析”明确阻断。进一步验收应对比冻结 n/P/CI 原值，避免文案层面禁止而实现层面意外改值。 |
| “把 Nature 的 170 mm 模板直接套到 Blood Advances” | 查该刊真实 profile 和官方证据，使用其允许尺寸；未知时只能输出 draft，不能冒充该刊合规。 | **PARTIAL**：已警告通用 profile 不是两目标期刊硬要求，但没有明确要求先取得官方尺寸/来源日期；实现者仍可能先套错尺寸再留下一个限制说明。 |
| “独立 panel 不带 key，解释放在整版图注里即可” | 关键色条/size key/cohort key/统计注释必须随独立 panel；缺失应视为硬失败并修复。 | **PARTIAL**：原则已要求必要 keys，操作规则只明确隔离裁剪污染，尚未逐项定义 semantic-completeness 的硬门槛。可从当前 F 的泛用 shared footer 看出，仅“有 key”不足以证明本 panel 的每种编码已解释。 |

## 两项最弱处及可直接采用的改文

**1. 目标期刊 profile 门槛。** 建议在执行顺序第 1 步后加：

> “生成 final_artwork 前必须记录 target_journal、官方 artwork URL、读取日期、允许宽高/字号/文件格式与本次选定尺寸。不得将 Nature/Science profile 默认应用于 Blood Advances 或 Scientific Reports。目标规则尚未核实时标记 draft/待验证，不能给期刊合规 PASS；合理的 Nature 借鉴样式仍可用于草稿。”

**2. 独立 panel 的语义完整性门槛。** 建议在可操作规则的 standalone 条目后加：

> “每个独立 panel 要登记 owner artists，包括 axes、continuation/colorbar、figure text、legend、rho/CI/n/P 及数值/尺寸 keys。隔离非所属对象后导出，并逐项核对：所有颜色/大小/线型/空心编码有对应解释，统计文字与冻结源结果一致，数值刻度与轴单位完整，无邻 panel 残留。缺任何必要 key/单位/统计文字均为硬失败；共享全图 keys 不代替所属编码检查。修复后需实际看 PNG，并记录检查范围。”

第二轮未编辑生产绘图程序或项目 Skill；这里给出主执行者可采纳的改文。未复核 Figure7 B/C 的统计注释或 Figure7 G 色标标题，旧审查对应未解决建议仍保留。


## Figure7 F 后续落实

上述共享无关 keys 与成人颜色释义 WARN 已修复：独立 panel 使用专属 Discovery／四亚型图例，TLX3 n=2 保留空心圆。已查看新 PNG 并验证 PDF 的各亚型与 n，整版 Figure7 保持不变。新文件 SHA256 和范围见 docs/PANEL7F_REVIEW.json；先前表中的 SHA256 记录的是改进前审计快照。
