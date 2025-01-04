import os
import scanpy as sc
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from rumaps.visualization import umapviz

# Define the path to test AnnData file and logo file
TEST_ADATA_PATH = "test_adata.h5ad"
LOGO_PATH = "test_logo.png"

# Load a test AnnData object
def load_test_adata():
    return sc.read_h5ad(TEST_ADATA_PATH)

def test_basic_umap_plot():
    """Test a basic UMAP plot."""
    adata = load_test_adata()
    umapviz(adata, labels_key="cluster_key", title="Basic UMAP Plot")

def test_save_plot():
    """Test saving a UMAP plot to a file."""
    adata = load_test_adata()
    output_path = "test_plot.png"
    umapviz(adata, labels_key="cluster_key", title="Save Test", save_path=output_path)
    assert os.path.exists(output_path), "Plot was not saved"
    os.remove(output_path)

def test_highlight_clusters():
    """Test highlighting specific clusters in the UMAP plot."""
    adata = load_test_adata()
    umapviz(
        adata,
        labels_key="cluster_key",
        title="Highlight Clusters",
        highlight_labels=["Cluster1", "Cluster2"],
        highlight_label_keywords={"fontsize": 12, "fontweight": "bold"}
    )

def test_custom_colormap():
    """Test a UMAP plot with a custom colormap."""
    adata = load_test_adata()
    custom_cmap = ListedColormap(cm.gist_earth(np.linspace(0.3, 0.9, 10)), name="CustomCMap")
    umapviz(
        adata,
        labels_key="cluster_key",
        title="Custom Colormap",
        custom_cmap=custom_cmap
    )

def test_dark_mode():
    """Test a UMAP plot in dark mode."""
    adata = load_test_adata()
    umapviz(
        adata,
        labels_key="cluster_key",
        title="Dark Mode Test",
        darkmode=True
    )

def test_add_logo():
    """Test adding a logo to the UMAP plot."""
    if not os.path.exists(LOGO_PATH):
        print(f"Logo file {LOGO_PATH} not found, skipping test.")
        return

    adata = load_test_adata()
    umapviz(
        adata,
        labels_key="cluster_key",
        title="Logo Test",
        logo_path=LOGO_PATH
    )

def test_dynamic_labels():
    """Test dynamic label sizing."""
    adata = load_test_adata()
    umapviz(
        adata,
        labels_key="cluster_key",
        title="Dynamic Labels",
        dynamic_label_size=True
    )

def test_adjust_point_size():
    """Test adjusting point size and transparency."""
    adata = load_test_adata()
    umapviz(
        adata,
        labels_key="cluster_key",
        title="Adjusted Points",
        point_size=10,
        alpha=0.5
    )

def test_prune_small_clusters():
    """Test pruning small clusters from the UMAP plot."""
    adata = load_test_adata()
    label_sizes = adata.obs["cluster_key"].value_counts()
    pruned_labels = adata.obs["cluster_key"].copy()
    pruned_labels[np.in1d(pruned_labels, label_sizes[label_sizes < 300].index)] = "Small Cluster"
    adata.obs["pruned_cluster_key"] = pruned_labels

    umapviz(
        adata,
        labels_key="pruned_cluster_key",
        title="Pruned Clusters"
    )

def test_medoid_labels():
    """Test using medoids for label placement."""
    adata = load_test_adata()
    umapviz(
        adata,
        labels_key="cluster_key",
        title="Medoid Labels",
        use_medoids=True
    )
