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