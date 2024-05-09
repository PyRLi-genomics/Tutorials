import pandas as pd
import matplotlib.pyplot as plt

def plot_volcano(file_path, id_col='gene_name', x_col='log2FoldChange', y_col='padj',
                 up_color='red', down_color='green', neutral_color='grey', label_genes=False,
                 genes_to_label=None, plot_title='Volcano plot', point_size=50, alpha=0.6,
                 x_lim=None, x_break=None, y_lim=None, y_break=None,
                 x_label_size=12, x_title_size=14, y_label_size=12, y_title_size=14,
                 label_size=10, log2FC_threshold=(1, -1), padj_threshold=0.05):
    """
    Generate a volcano plot for gene expression data.

    Parameters:
    - file_path: Path to the CSV file containing the data.
    - id_col: Column name for gene identifiers.
    - x_col: Column name for log2 fold change values.
    - y_col: Column name for adjusted p-values.
    - up_color, down_color, neutral_color: Colors for upregulated, downregulated, and not significant points.
    - label_genes: Boolean to indicate if genes should be labeled.
    - genes_to_label: List of gene names to label or number to label the top changing genes.
    - plot_title: Title of the plot.
    - point_size: Size of the points in the plot.
    - alpha: Transparency of the points.
    - x_lim, y_lim: Tuple for x and y axis limits (min, max).
    - x_break, y_break: Break size on the x and y axis.
    - x_label_size, y_label_size: Font size for axis labels.
    - x_title_size, y_title_size: Font size for axis titles.
    - label_size: Font size for data point labels.
    - log2FC_threshold: Thresholds for significant log2 fold change.
    - padj_threshold: Threshold for significant adjusted p-value.
    """

    # Load data
    data = pd.read_csv(file_path)

    # Determine colors for each data point based on significance thresholds
    conditions = [
        (data[x_col] >= log2FC_threshold[0]) & (data[y_col] <= padj_threshold),
        (data[x_col] <= log2FC_threshold[1]) & (data[y_col] <= padj_threshold)
    ]
    choices = [up_color, down_color]
    data['color'] = pd.Series(neutral_color, index=data.index)
    data['color'] = pd.Series(np.select(conditions, choices, default=neutral_color), index=data.index)

    # Plotting
    plt.figure(figsize=(10, 8))
    plt.scatter(data[x_col], -np.log10(data[y_col]), c=data['color'], alpha=alpha, s=point_size)

    # Label genes if required
    if label_genes:
        if genes_to_label is None or type(genes_to_label) == int:
            top_genes = data.nlargest(genes_to_label or 10, key=lambda x: abs(x[x_col]))
        else:
            top_genes = data[data[id_col].isin(genes_to_label)]

        for _, row in top_genes.iterrows():
            plt.text(row[x_col], -np.log10(row[y_col]), row[id_col], fontsize=label_size)

    # Draw significance threshold lines
    plt.axhline(-np.log10(padj_threshold), color='grey', linestyle='dashed')
    plt.axvline(log2FC_threshold[0], color='grey', linestyle='dashed')
    plt.axvline(log2FC_threshold[1], color='grey', linestyle='dashed')

    # Set limits and labels
    plt.xlim(x_lim or (data[x_col].min() - x_break if x_break else None, data[x_col].max() + x_break if x_break else None))
    plt.ylim(y_lim or (-np.log10(data[y_col].max()) - y_break if y_break else None, 0 + y_break if y_break else None))
    plt.xlabel('Log2 Fold Change', fontsize=x_title_size)
    plt.ylabel('-Log10 adjusted P-value', fontsize=y_title_size)
    plt.title(plot_title)

    # Show plot
    plt.show()

# Example usage
plot_volcano('path_to_your_data.csv')
