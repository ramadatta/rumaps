import datamapplot as dmp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from PIL import Image
import pandas as pd

def umapviz(
    adata,
    labels_key,
    title="UMAP Plot",
    sub_title=None,
    title_keywords=None,
    sub_title_keywords=None,
    point_size=15,
    marker_type="o",
    alpha=0.8,
    cmap=None,
    darkmode=False,
    dynamic_label_size=False,
    label_over_points=False,
    highlight_labels=None,
    highlight_label_keywords=None,
    label_margin_factor=1.0,
    font_size=12,
    font_family="Arial",
    font_weight=None,
    label_linespacing=None,
    label_wrap_width=20,
    arrowprops=None,
    use_medoids=False,
    add_glow=True,
    noise_label=None,
    noise_color="gray",  # Default color for noise clusters
    custom_cmap=None,
    logo_path=None,
    save_path=None,
    prune_clusters=False,
    min_cluster_size=0,
    small_cluster_label="Small Cluster",
    figsize=(10, 8),
    save_format="pdf",
):
    """
    Generate a highly customizable UMAP plot.
    """
    umap_data = adata.obsm["X_umap"]
    labels = adata.obs[labels_key].copy()

    # Ensure labels are categorical
    if not isinstance(labels, pd.Categorical):
        labels = pd.Categorical(labels)

    # Ensure "Small Cluster" is a valid category
    if small_cluster_label not in labels.categories:
        labels = labels.add_categories([small_cluster_label])

    # Handle colors from adata.uns
    label_color_map = None
    if custom_cmap:
        print("Using custom colormap provided by the user.")
        cmap = custom_cmap
    elif f"{labels_key}_colors" in adata.uns:
        # Create label-to-color mapping
        label_categories = adata.obs[labels_key].cat.categories
        colors = adata.uns[f"{labels_key}_colors"]
        label_color_map = dict(zip(label_categories, colors))
        print(f"Using colors from adata.uns['{labels_key}_colors']")
    else:
        cmap = plt.cm.tab20
        print("Using default colormap (tab20).")

    # Handle logo
    logo = None
    if logo_path:
        logo = np.asarray(Image.open(logo_path))

    # Ensure arrowprops is valid
    if arrowprops is None:
        arrowprops = {}  # Initialize as an empty dictionary

    # Prune small clusters
    if prune_clusters:
        label_sizes = pd.Series(labels).value_counts()
        small_clusters = label_sizes[label_sizes < min_cluster_size].index
        labels[np.in1d(labels, small_clusters)] = small_cluster_label

    # Grey out clusters specified in noise_label
    if isinstance(noise_label, list):
        for cluster in noise_label:
            labels[labels == cluster] = small_cluster_label
    elif isinstance(noise_label, str):
        labels[labels == noise_label] = small_cluster_label

    # Ensure highlight_label_keywords is valid
    if highlight_label_keywords is None:
        highlight_label_keywords = {}  # Initialize as an empty dictionary

    # Dynamic label size setup
    if dynamic_label_size:
        font_size = None  # Let `datamapplot` handle dynamic sizing internally

    # Plot the UMAP with the specified parameters
    fig, ax = dmp.create_plot(
        umap_data,
        labels,
        title=title,
        sub_title=sub_title,
        title_keywords=title_keywords,
        sub_title_keywords=sub_title_keywords,
        point_size=point_size,
        marker_type=marker_type,
        alpha=alpha,
        cmap=cmap,
        darkmode=darkmode,
        dynamic_label_size=dynamic_label_size,
        label_over_points=label_over_points,
        label_font_size=font_size if not dynamic_label_size else None,
        highlight_labels=highlight_labels,
        highlight_label_keywords=highlight_label_keywords,
        label_margin_factor=label_margin_factor,
        font_family=font_family,
        font_weight=font_weight,
        label_linespacing=label_linespacing,
        label_wrap_width=label_wrap_width,
        arrowprops=arrowprops,
        use_medoids=use_medoids,
        add_glow=add_glow,
        noise_label=small_cluster_label,
        noise_color=noise_color,
        label_color_map=label_color_map,  # Pass the label-to-color mapping
        logo=logo,
        figsize=figsize,
    )

    # Save the plot if requested
    if save_path:
        fig.savefig(save_path, format=save_format, bbox_inches="tight")
        print(f"Plot saved to {save_path}")

    plt.show()
    return fig, ax
