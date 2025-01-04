
import datamapplot as dmp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns

def umap_iviz(
    adata,
    labels_key,
    label_layers=None,
    title="Interactive UMAP",
    sub_title=None,
    hover_columns=None,
    enable_search=True,
    save_path=None,
    custom_cmap=None,
    colormaps=None,  # Simpler dictionary for dynamic coloring
    colormap_rawdata=None,  # List of raw data arrays for advanced customization
    colormap_metadata=None,  # List of metadata dictionaries for colormap control
):
    """
    Generate an interactive UMAP plot for AnnData objects with hover text and advanced colormap support.
    
    Parameters:
    - adata: AnnData object containing the data.
    - labels_key: Primary categorical field for coloring.
    - label_layers: Additional label layers for secondary annotations.
    - title: Title of the interactive plot.
    - sub_title: Subtitle of the interactive plot.
    - hover_columns: Columns to include in the hover text.
    - enable_search: Enable search functionality in the interactive plot.
    - save_path: Path to save the plot.
    - custom_cmap: Custom colormap for primary coloring.
    - colormaps: Dictionary of additional fields for dynamic colormap support.
    - colormap_rawdata: List of raw data arrays for advanced colormap customization.
    - colormap_metadata: List of metadata dictionaries for advanced colormap customization.
    """
    # UMAP coordinates and labels
    umap_data = adata.obsm["X_umap"]
    labels = adata.obs[labels_key]

    # Check for colors in adata.uns
    label_color_map = None
    if f"{labels_key}_colors" in adata.uns:
        label_categories = adata.obs[labels_key].cat.categories
        colors = adata.uns[f"{labels_key}_colors"]
        label_color_map = dict(zip(label_categories, colors))
        print(f"Using colors from adata.uns['{labels_key}_colors']")
    else:
        print(f"No predefined colors found for '{labels_key}', using default colormap.")
        cmap = plt.cm.tab20 if custom_cmap is None else custom_cmap
        unique_labels = labels.unique()
        label_color_map = {label: cmap(i / len(unique_labels)) for i, label in enumerate(unique_labels)}

    # Process additional label layers
    processed_label_layers = []
    if label_layers:
        label_layers = [label_layers] if isinstance(label_layers, str) else label_layers
        for layer in label_layers:
            if layer in adata.obs:
                processed_label_layers.append(adata.obs[layer])
            else:
                raise KeyError(f"[Error] Layer '{layer}' not found in adata.obs")

    # Generate hover text
    hover_data = (
        adata.obs[hover_columns].apply(lambda row: "\n".join([f"{col}: {row[col]}" for col in hover_columns]), axis=1)
        if hover_columns else
        None  # No hover text if hover_columns is not provided
    )

    # Validate colormap_rawdata and colormap_metadata
    if colormap_rawdata and colormap_metadata:
        if len(colormap_rawdata) != len(colormap_metadata):
            raise ValueError("[Error] 'colormap_rawdata' and 'colormap_metadata' must have the same length.")
        print(f"Using advanced colormap control with {len(colormap_rawdata)} fields.")

    # Create interactive plot
    plot = dmp.create_interactive_plot(
        umap_data,
        labels,
        *processed_label_layers,
        hover_text=hover_data,
        title=title,
        sub_title=sub_title,
        cmap=custom_cmap,  # Use user-provided colormap
        enable_search=enable_search,
        label_color_map=label_color_map,  # Apply predefined or fallback colors
        colormaps=colormaps,  # Simpler dynamic colormaps
        colormap_rawdata=colormap_rawdata,  # Advanced colormap customization
        colormap_metadata=colormap_metadata,  # Metadata for colormap customization
    )

    # Save plot if requested
    if save_path:
        plot.save(save_path)

    return plot



# def umap_iviz(
#     adata,
#     labels_key,
#     label_layers=None,
#     title="Interactive UMAP",
#     sub_title=None,
#     hover_columns=None,
#     enable_search=True,
#     save_path=None,
#     custom_cmap=None,
#     colormaps=None,  # Directly pass colormaps dictionary
# ):
#     """
#     Generate an interactive UMAP plot for AnnData objects with hover text and colormap support.
    
#     Parameters:
#     - adata: AnnData object containing the data.
#     - labels_key: Primary categorical field for coloring.
#     - label_layers: Additional label layers for secondary annotations.
#     - title: Title of the interactive plot.
#     - sub_title: Subtitle of the interactive plot.
#     - hover_columns: Columns to include in the hover text.
#     - enable_search: Enable search functionality in the interactive plot.
#     - save_path: Path to save the plot.
#     - custom_cmap: Custom colormap for primary coloring.
#     - colormaps: Dictionary of additional fields for dynamic colormap support.
#     """
#     # UMAP coordinates and labels
#     umap_data = adata.obsm["X_umap"]
#     labels = adata.obs[labels_key]

#     # Check for colors in adata.uns
#     label_color_map = None
#     if f"{labels_key}_colors" in adata.uns:
#         label_categories = adata.obs[labels_key].cat.categories
#         colors = adata.uns[f"{labels_key}_colors"]
#         label_color_map = dict(zip(label_categories, colors))
#         print(f"Using colors from adata.uns['{labels_key}_colors']")
#     else:
#         # Fallback to default colormap if no colors are found
#         print(f"No predefined colors found for '{labels_key}', using default colormap.")
#         cmap = plt.cm.tab20 if custom_cmap is None else custom_cmap
#         unique_labels = labels.unique()
#         label_color_map = {label: cmap(i / len(unique_labels)) for i, label in enumerate(unique_labels)}

#     # Process additional label layers
#     processed_label_layers = []
#     if label_layers:
#         label_layers = [label_layers] if isinstance(label_layers, str) else label_layers
#         for layer in label_layers:
#             if layer in adata.obs:
#                 processed_label_layers.append(adata.obs[layer])
#             else:
#                 raise KeyError(f"[Error] Layer '{layer}' not found in adata.obs")

#     # Generate hover text
#     hover_data = (
#         adata.obs[hover_columns].apply(lambda row: "\n".join([f"{col}: {row[col]}" for col in hover_columns]), axis=1)
#         if hover_columns else
#         None  # No hover text if hover_columns is not provided
#     )

#     # Create interactive plot
#     plot = dmp.create_interactive_plot(
#         umap_data,
#         labels,
#         *processed_label_layers,  # Include processed label layers
#         hover_text=hover_data,
#         title=title,
#         sub_title=sub_title,
#         cmap=custom_cmap,  # Use user-provided colormap
#         enable_search=enable_search,
#         label_color_map=label_color_map,  # Apply predefined or fallback colors
#         colormaps=colormaps,  # Pass additional colormaps
#     )

#     # Save plot if requested
#     if save_path:
#         plot.save(save_path)

#     return plot
