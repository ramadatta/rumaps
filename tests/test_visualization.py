import scanpy as sc
from rumaps.visualization import static_umap, interactive_umap

def test_static_umap():
    adata = sc.datasets.pbmc3k()
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    static_umap(adata, labels_key="louvain")

def test_interactive_umap():
    adata = sc.datasets.pbmc3k()
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    plot = interactive_umap(adata, labels_key="louvain")
    assert plot is not None
