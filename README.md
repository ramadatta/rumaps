
# RUMAPS

RUMAPS is a Python package designed to create refined UMAP visualizations with clean labeling and interactive exploration. Specifically, this package is created to fit in single cell analysis packages (scanpy and anndata).  RUMAPS provides intuitive tools to enhance your data visualization workflow.



# Static plot using RUMAPS
![static](https://github.com/ramadatta/rumaps/blob/main/images/static_plot.png)

# Interactive plot using RUMAPS
![interactive](images/interactive.gif)

---

## Key Features
- **Polished UMAP Visualizations**: Generate clean and well-structured UMAP plots with custom labeling.
- **Interactive Exploration**: Add hover text, search functionality, and custom layers to dive deeper into your data.
- **Streamlined Integration**: Seamlessly integrate RUMAPS into your existing data analysis pipelines with scanpy and anndata object.
- **Interactive UMAP HTML:**: Interactive UMAP visualizations can be saved and shared as HTML files, enabling ease of sharing and exploration.

---

## Installation

<!-- Install RUMAPS via pip: -->

<!-- ```bash -->
<!-- pip install -i https://test.pypi.org/simple/ rumaps==0.1.2 -->
<!-- ``` -->

Install RUMAPS directly with pip
```
pip install git+https://github.com/ramadatta/rumaps.git
```

Install RUMAPS in a conda environment:
```
conda create -n rumaps
conda activate rumaps
pip install "matplotlib>=3.8"
pip install datashader
pip install "colorspacious>=1.1"
pip install "scikit-image>=0.22"
pip install dask
python -m pip install "dask[dataframe]" --upgrade
pip install git+https://github.com/TutteInstitute/datamapplot.git
```

Note: Some functionalities are not available in pip version package of ```datamapplot```. So, one may need to install from the main branch on github to have ```rumaps``` working.
```
pip install git+https://github.com/TutteInstitute/datamapplot.git
```

Install RUMAPS from Source Code:
```bash
git clone https://github.com/ramadatta/rumaps.git
cd rumaps
pip install -e .
```


---

## Usage

Documentation
For detailed usage examples, please check the tutorial notebook:

[RUMAPS Tutorial Notebook](https://colab.research.google.com/drive/18SynVkqi3sw7ZSXyTUu_PIarvZr-mlQV?usp=sharing)


---

## Dependencies
RUMAPS requires the following Python packages:
- `numpy`
- `pandas`
- `matplotlib`
- `datamapplot`
- `scanpy`
- `seaborn`
- `Pillow>=9.0.0`


Install these dependencies using:

```bash
pip install -r requirements.txt
```

---

## Contributing

We welcome contributions from the community! If you have suggestions, bug fixes, or enhancements, feel free to open an issue or submit a pull request.

---

## License

RUMAPS is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

---

## Acknowledgments

- RUMAPS was developed using the datamapplot library. Special thanks to the team behind DataMapPlot for their amazing work.
- Publicly available [Tabula Muris h5ad](https://figshare.com/articles/dataset/Tabula_Muris_Scanpy/13363628?file=25753577) file was used to demonstration purposes.

---

## Contact

For any questions or issues, please contact us at [srsridatta.prakki@helmholtz-munich.de].
