import scanpy.tools._score_genes as m
print(m.__file__)
names = [n for n in dir(m) if "gene" in n.lower() or n in {"s_genes", "g2m_genes"}]
print(names)
