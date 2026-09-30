# 本轮技能应用与已完成改进

女娲主题流程使用三组并行视角，实际读取 R4DS/Wickham、Wilke 原著，以及 Wilson/Sandve/GitHub 官方文档，共8个一手来源。以主题框架而非人物模仿提炼六个决策模型。每组报告、原始页面/正文、访问记录在项目 `.agents/skills/zeb1-reproducibility-framework/references/`，方法分歧和局限明确保留。

nature-science-figure-skill 已用于 final-artwork 的字体/尺度/vector/standalone QA。实际目标是 Blood Advances / Scientific Reports，因此不直接套 Nature 的宽高上限，也不宣称本项目通过 Nature 正式合规。当前冻结语义色采用明确自定义 `zeb_semantic_locked`，保留已有 ZEB1/ZEB2/ETP 语义和图注；不把近似 article-inspired 色板描述成期刊官方推荐。

完成的具体改进：

- Article、Image、data、script 固定入口；旧版本移到项目外可恢复备份。
- 最新版20张图/92个实际 panel 的完整入口，R 调用 Python，并保留PDF/PNG/SVG。
- S4D恢复基于gene key的修订数据；S6恢复定向修正的 GSE146901/恶性效应源；补图标题/字母按实际已发布母版还原。
- panel真正隔离所属axes/text再导出，色条完整，shared figure keys单独无其它图形渲染；不带下一 panel 字母或图形。
- 单panel运行不会改 Image整图；--out期刊导出不会改通用panel。
- 本地路径自动寻找项目根；R/Python可配置；实际软件记录和输入SHA。无隐式下载/安装。
- 数据集来源、源码QC行、处理入口与已知限制可查询；不掩盖Harmony未收敛/版本许可未核实。
- GitHub按用户模板筛选代码与非识别元数据；矩阵/患者表/旧稿/缓存排除；README明确可运行范围。

整图在整理前后逐像素一致，是没有误改现成数据图的验收；独立panel导出与复现工程则做了实质优化。新的视觉/统计改进仍可在这套固定入口迭代，无需再创建一份文件夹。
