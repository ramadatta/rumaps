import datamapplot as dmp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def interactive_umap(
    adata,
    labels_key,
    label_layers=None,
    title="Interactive UMAP",
    sub_title=None,
    hover_columns=None,
    enable_search=True,
    save_path=None,
    custom_cmap=None,
):
    """
    Generate an interactive UMAP plot for AnnData objects with hover text support.

    Parameters:
        adata: AnnData object containing UMAP data and labels.
        labels_key: Key in `adata.obs` for cluster labels.
        label_layers: Additional label layers to include (list of strings or single string).
        title: Title for the plot.
        hover_columns: List of columns to include in hover text.
        enable_search: Enable search functionality.
        save_path: Path to save the HTML output.
    """
    print("Debug: Start of interactive_umap function")
    print(f"UMAP data shape: {adata.obsm['X_umap'].shape}")
    print(f"Labels key: {labels_key}")

    umap_data = adata.obsm["X_umap"]
    labels = adata.obs[labels_key]

    # Generate colormap
    cmap = plt.cm.tab20 if custom_cmap is None else custom_cmap
    unique_labels = labels.unique()
    label_color_map = {label: cmap(i / len(unique_labels)) for i, label in enumerate(unique_labels)}

    print(f"Label color map: {label_color_map}")

    # Use existing colors for the primary labels if available
    if f"{labels_key}_colors" in adata.uns:
        colors = adata.uns[f"{labels_key}_colors"]
        label_color_map = dict(zip(unique_labels, colors))
    print(f"Primary label color map: {label_color_map}")

    # Process additional layers and assign colors
    processed_label_layers = []
    if label_layers:
        for layer in label_layers:
            if layer in adata.obs:
                layer_data = adata.obs[layer]
                processed_label_layers.append(layer_data)
                
                # Check for colors in adata.uns
                if f"{layer}_colors" in adata.uns:
                    layer_colors = adata.uns[f"{layer}_colors"]
                    unique_layer_labels = layer_data.unique()
                    layer_color_map = dict(zip(unique_layer_labels, layer_colors))
                    
                    # Add these colors to the main color map
                    for label, color in layer_color_map.items():
                        if label not in label_color_map:
                            label_color_map[label] = color
                else:
                    # Assign new colors if none exist in adata.uns
                    unique_layer_labels = layer_data.unique()
                    for i, label in enumerate(unique_layer_labels):
                        if label not in label_color_map:
                            label_color_map[label] = cmap((len(label_color_map) + i) / len(unique_labels))
            else:
                raise KeyError(f"Layer '{layer}' not found in adata.obs")
    
    print(f"Final label color map: {label_color_map}")


    # Generate hover text
    if hover_columns:
        print(f"Generating hover_data from columns: {hover_columns}")
        hover_data = adata.obs[hover_columns].apply(
            lambda row: "\n".join([f"{col}: {row[col]}" for col in hover_columns]), axis=1
        )
    else:
        print("No hover_columns provided, hover text will be empty.")
        hover_data = pd.Series([""] * len(labels), index=adata.obs.index)

    # Create interactive plot
    plot = dmp.create_interactive_plot(
        umap_data,
        labels,
        *processed_label_layers,
        # *label_layers if label_layers else [],
        hover_text=hover_data,
        title=title,
        sub_title=sub_title,
        cmap=cmap,
        label_color_map=label_color_map,
        enable_search=enable_search,
    )

    if save_path:
        plot.save(save_path)
        print(f"Interactive plot saved to {save_path}")

    print("Debug: Plot creation complete")
    return plot
