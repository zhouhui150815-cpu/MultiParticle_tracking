# MultiParticle_tracking
Java implementation for automated video detachment

### Overview
MultiParticle_Tracking is an ImageJ/Fiji plugin designed for tracking cell/particle detachment events in time-lapse image sequences. It automatically detects cells within a user-defined pentagonal probe region, tracks their presence across frames, calculates detachment-related parameters, generates data tables, and creates pseudo-color heatmaps based on detachment time differences.

### Features
- **Integrated UI**: All parameters in one window with clear workflow
- **Flexible Probe Region**: Draw pentagon probes or load existing ROIs
- **Automatic Cell Detection**: Otsu thresholding with morphological operations
- **Detachment Time Analysis**: Track cells from first frame across subsequent frames
- **Pseudo-color Heatmap**: Visualize detachment time distribution with color mapping
- **Comprehensive Statistics**: Generate tables with area, perimeter, circularity, detachment time, etc.

### Installation
1. Copy the `MultiParticle_Tracking.java` file to your ImageJ/Fiji `plugins` folder
2. Restart ImageJ/Fiji
3. The plugin will appear in the `Plugins` menu as "MultiParticle_Tracking"

### Usage Guide
1. **Open Video/Sequence**
   - Click "Open Video/Sequence" button
   - Select your time-series image (TIFF stacks, AVI, etc.)

2. **Set Scale**
   - Click "Set" in the Scale Setting panel
   - Draw a line of known length using the line tool
   - Click "Confirm" and enter the actual length in micrometers

3. **Define Probe Region**
   - Click "Draw Pentagon Probe" to draw a pentagon region
   - Or click "Load Existing ROI" to load a pentagon ROI from ROI Manager
   - The probe will be displayed as a red outline (not added to ROI Manager)

4. **Configure Parameters**
   - **Min/Max Area (μm²)**: Filter cells by size
   - **Min/Max Circularity**: Filter cells by shape (0-1, 1=perfect circle)
   - **Cell Distance Deviation (%)**: Allowed distance from probe (0-10000%)

5. **Preview Detection**
   - Click "Preview First Frame" button
   - Cells near probe appear in purple (eligible), distant cells in yellow (ineligible)
   - Adjust parameters based on preview results

6. **Start Analysis**
   - Click "Start Tracking Analysis" button
   - The plugin will track cells across all frames
   - Results table and pseudo-color settings will appear automatically

7. **View Results**
   - **Results Table**: Shows all eligible cells with parameters
   - **Pseudo-color Settings**: Adjust LUT type and thresholds
   - **Analysis Summary**: Provides statistical overview

8. (optional)**Pseudo-color Settings**
- **Select LUT Type**: Jet, Spectrum, Red, Green, Blue
- **Set Detachment Constant C**: Default 10.78, scales color mapping
- **Adjust Threshold Range**: Set min/max thresholds (seconds)
- **Click "Confirm" to regenerate the pseudo-color map.

### Output Description

#### Results Table Columns
| Column | Description |
|--------|-------------|
| Cell ID | Cell identifier (format: cell1, cell2, ...) |
| Initial Area (μm²) | Initial area |
| Initial Perimeter (μm) | Initial perimeter |
| Initial Circularity | Initial circularity |
| Distance Deviation (%) | Distance deviation from probe |
| Detachment Time (s) | Detachment time ("Not detached" if never detached) |
| Detachment Frame | Detachment frame number ("Not detached" if never detached) |

#### Pseudo-color Map
- White background
- Colored regions: Color-mapped based on detachment time
- Black outline: Probe region

#### Analysis Summary
- Analysis parameters (scale, probe region, thresholds, etc.)
- Statistics (total cells, detached cells, detachment ratio, average/max detachment time)

### Parameter Details

| Parameter | Description | Recommended Range |
|-----------|-------------|-------------------|
| Min Area | Minimum cell area (exclude noise) | Adjust based on cell size |
| Max Area | Maximum cell area (exclude large artifacts) | Adjust based on cell size |
| Min Circularity | Minimum circularity (exclude irregular shapes) | 0-1 |
| Max Circularity | Maximum circularity | 0-1 |
| Distance Deviation | Allowed distance from probe | 0-1000% |
| Frame Interval | Time between frames | Set according to acquisition |
| Detachment Constant C | Color mapping scale factor | 10.78 (default) |

### Troubleshooting

#### 1. No cells detected
- Check area thresholds (may be too large or too small)
- Ensure image is 8-bit grayscale (plugin auto-converts)
- Adjust circularity thresholds
- Verify probe region is drawn in cell-dense area

#### 2. Too many cells detected (including noise)
- Increase minimum area threshold
- Increase minimum circularity requirement
- Adjust distance deviation threshold

#### 3. No table appears after analysis
- Check console for error messages
- Ensure at least one eligible cell was detected
- Try restarting ImageJ/Fiji

#### 4. Pseudo-color map not showing or abnormal
- Ensure at least one cell was successfully tracked
- Check ROI Manager for correctly named ROIs (format: 1f-cellX)
- Adjust threshold ranges

### Development

#### Project Structure
```
MultiParticle_Tracking/
├── src/
│   └── main/
│       └── java/
│           └── MultiParticle_Tracking.java
├── target/
│   ├── classes/
│   ├── MultiParticle_Tracking-1.0.0.jar
│   └── ...
├── pom.xml
└── README.md
```

#### Dependencies
- ImageJ
- Maven (development only)

### Author
**Liu Ruye, Sichuan University**
