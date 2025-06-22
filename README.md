# Fura2 Imaging & Analysis Pipeline

A multi‐step MATLAB workflow:

## 1. Split 340/380 into Cellpose stacks
- **Inputs**: `Input.mat` (positions), `imageDir`
- **Process**:  
  - Find `_340` & `_380` multi‐page TIFFs per position  
  - Combine corresponding planes (380→R, 340→G), median‐filter, rescale  
  - Save `stackImages/{position}_imPlane{n}_cp.tif`
- **Output**: Cellpose‐ready stacks

## 2. Merge Cellpose masks
- **Inputs**: `cellMaskDir` masks, `imageDir` stacks, `input` list
- **Process**:  
  - For each sample, collect `_imPlane{n}_cp_cp_masks.tif`  
  - Append into `{sample}_cp_masks.tif`
- **Output**: Stacked mask TIFFs

## 3. Timelapse analysis
- **Inputs**: `cpDir` masks, `imDir` raw images, `input` prefixes  
- **Steps**:  
  1. **Normalize** object IDs across planes → `cpObjectNums.mat`  
  2. **Extract** f340/f380/AMT/GFP intensities per object & frame (background-corrected)  
  3. **Save** `singleOutput` with per-object time series
- **Outputs**: `backCoord.mat` (if manual), `singleOutput.mat`

## 4. Interactive gating & metrics
- **Inputs**: `singleOutput`, `input` layout  
- **Process**:  
  - Histograms & polygon ROI for:  
    - Object area  
    - Total Fura2 intensity  
    - Ratiometric trend  
    - NGFR vs AMT  
  - Compute average ratio per sample  
  - Export gating plots as PDF  
- **Outputs**: `avgOutput`, gated subsets (`negData*`), PDF figures

## 5. Validate mask alignment
- **Inputs**: `cpObjectNums.mat`, raw & mask dirs  
- **Process**:  
  - Re‐create merged intensity stacks  
  - Mask out non-aligned objects  
  - Save `{sample}_intensityImage.tif`
- **Output**: Per-sample intensity TIFFs

## 6. Export to Excel
- **Inputs**: `singleOutput`, `input`  
- **Steps**:  
  - Flatten single-cell values → `singleCellOutput.xlsx`  
  - Pivot for GraphPad → `output.xlsx`
- **Outputs**: Excel files

## 7. Batch rename & copy
- **Input**: `tmp` folder of TIFFs  
- **Process**: copy to `newfolder`, insert `Iono` in filenames  
- **Output**: Renamed files in `newfolder`
