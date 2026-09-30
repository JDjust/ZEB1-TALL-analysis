# 从已有数据到图：操作路线

1. 在 `dataset_index.tsv` 选择 accession，查看原 `dataset_registry.tsv` 来源与许可状态；关联是现有源码/路径证据，不等于所有条目已完成独立验证。
2. 在 `metadata/code_dataset_index.tsv` 找处理脚本；在 `metadata/qc_code_evidence.tsv` 找过滤/转换/样本合并的实际代码行。`analysis/modules/module*/report.md` 和各 validation README/manifest 记录原执行结果。
3. 将 counts、TPM/FPKM、logCPM 及显示 z 分开；检查 gene_id/version、患者或供者键、重复样本、排除名单、pairing、批次与亚型混杂。依据本数据集阈值及分布，不复制统一“标准QC”。
4. 单细胞步骤还要检查 QC 是作者已做还是本项目重做、doublet/ambient 是否有证据、Harmony 收敛/降级、恶性标签来源、counts pseudobulk 和供者层推断。当前已证实的问题见 QC_LIMITATIONS.md。
5. 若只需画图：直接运行 `script/`。64个实际绘图读表输入均有 SHA256；R入口不安装包/下载数据/重拟合模型。
6. 若要重新处理：在独立任务明确输入、服务器路径、软件环境与输出范围，先核对现有冻结结果，再决定替换。保留可回滚版本和原始输入；记录每阶段纳入/排除数、版本与随机种子。

本轮没有把 Hulu 原始矩阵拉回。data 中的大对象是整理前已有本地文件，保留其存在并排除公开仓库。预处理源中的服务器路径仅作 locator；没有将“存在路径”写成“服务器已核验”。
