# 已核对的 QC 与复现边界

这里记录对旧代码/日志的观察；本轮没有重跑上游处理。完整证据见 `.agents/skills/zeb1-reproducibility-framework/references/research/01-reproducibility.md`。

| 数据/问题 | 现有执行证据 | 仍需保留的限制 |
|---|---|---|
| GSE227122 | module5 代码与原日志：41,316 初始细胞；31,297 是 doublet 排除前过滤后数量；121 predicted doublets。阈值/顺序见原脚本与日志 | Harmony 在 10 次迭代后未收敛；原 “harmony OK” 只是返回坐标。恶性标签是项目启发式，不能因为 pseudobulk 用未整合 counts 就宣称标签不受整合影响 |
| GSE227122 / GSE248287 | 两队列 × 五细胞数阈值，共10行现成敏感性结果；10人与15人以患者为单位 | 阈值保留同一患者时不构成独立验证；区间与方向限制保留；不重新推断肿瘤比例 |
| TARGET STAR | 265 行代理 QC 表，文件 SHA 与原 manifest 相符 | 汇总 count QC 不是 FASTQ/BAM 比对/读段 QC；未独立重跑 Cox |
| GSE144035 | 6 对模型×2基因的12行配对核查已存在 | VTL 条件/模型不是自动的独立生物学重复；旧错误比较只能作为撤回来源保留 |
| GSE287751 | barcode 匹配和条件 pseudobulk 有既有执行记录 | 作者 HTO/基因数/mito QC 与本项目独立执行应分开；不同 context 不当作同一比较重复 |
| S4D 键对齐 | 已有 KEYED 表24356行，unique gene_id/symbol；使用修订图源 | 旧位置配对源仍作为修复输入/审计证据；绘图必须用 KEYED 表。没有重拟合 PC 模型 |
| DepMap Public26Q1 | 既有四个输入 hash 已锁定 | 官方 release 身份、发布日期和许可尚未独立核实；不能用其它 Figshare 记录代替 |
| 队列重叠 | 只在可兼容 ID namespace 做过审计 | 未能映射的队列不能声称零重叠 |
| 复现范围 | 20张当前图从冻结输入成功重绘 | 上游大数据/环境不在本轮重新运行范围；保留代码不等于全链验证 |

未来重新处理时必须记录 convergence/回退/缓存 fingerprint。当前冻结分析不因目录整理而自动被覆写。研究结论增强需要真实新证据，不能由更美观图版或更整齐目录代替。
