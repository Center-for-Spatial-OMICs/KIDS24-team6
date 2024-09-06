### Part I : GPU Grid-Based Pixel Classification


- Input: 
  - High-resolution H&E tif image [Download from GitHUb](https://github.com/stjude-biohackathon/KIDS24-team6/blob/main/Data/tissue_hires_image.tiff)
  - Grid size parameter for dividing the image
  - Number of clusters for KMeans classification
- Output: 
  - Classified grids with coordinates marked on the original H&E tif image
  - Visualization of clustered grids overlaid on the H&E image

### Part II: Match Barcodes sequencing information to H&E

- Input:
  - Visium HD spatial transcriptomics data (barcode file, tissue positions, gene expression)
  - High-resolution H&E pathology image
  - scalefactors_json.json file for scaling factor adjustments

- Output:
  - Annotated H&E image with overlaid barcode locations and gene expression data.

## Part III: Interactive RShiny App 

[Copy it to Run](https://github.com/stjude-biohackathon/KIDS24-team6/blob/main/Code/Shiny_app.R)


- **Upload Files:**  
  - Upload a CSV file and a PNG image to initiate the analysis.

- **Expression Cluster Selection:**  
  - After uploading files, select a specific expression cluster to view on spatial coordinates.  
  - Option to display all clusters together in a graph.

- **Scale and Offset Adjustments:**  
  - Adjust the size of the dots overlaid on the image using X/Y Scale Factor tabs.  
  - Shift the position of the dots with X/Y Offset controls.

- **Rotation and Flipping:**  
  - Rotate dots 90° counterclockwise.  
  - Flip the dots horizontally or vertically.

- **Pixel Cluster Selection:**  
  - Choose a specific pixel cluster to view its distribution.  
  - Option to plot all pixel clusters together.

- **Seurat Section:**  
  - View expression and pixel clusters overlaid on the HE image.  
  - Option to show all clusters together.

- **Gene Expression Visualization:**  
  - Enter a gene name to view its spatial expression on the image.

- **Reset Adjustments:**  
  - Reset all markers, buttons, and inputs back to their default state.