# 环境记录范围

requirements.txt 记录已实际用于本轮重绘的 Python 包。renv.lock 记录已安装的八个 R 包版本，未由 renv 完整解析递归依赖，未测试 renv::restore；不是已经验证的上游计算环境锁。R 绘图入口只用 base R 调用 Python，已用 Rscript --vanilla 验证。environment.yml 未在干净机器创建。单细胞/空间处理所需额外包应按对应原脚本逐数据集核对，不能从当前重绘成功推定安装完整。
