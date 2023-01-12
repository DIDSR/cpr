# Color Performance Review Tool for Endoscopy Devices

## 1. Introduction
The Color Performance Review (CPR) Tool for Endoscopy Devices is a software program written in Matlab for analyzing color performance testing data in endoscopy device. The tool accepts the testing data and then generates quantitative analyses for the user to compare color performance between devices. The input testing data include the measurement data of a standard test target from the device output as well as the ground truth of the test target. 

The tool provides the following analyses: 
*	Visual simulation of the test target and sample scenes
*	Absolute color errors in comparison with the ground truth
*	Preservation of the patch order in lightness, hue, and chroma
*	Preservation of color contrast between patches

### 1.1 Software Requirements
*	Matlab Runtime 9.12 (Version 9.12 tested)
*	Matlab (Version 9.12 tested)
*	Image Processing Toolbox (Version 11.5 tested)
*	Computer Vision Toolbox (Version 9.12 tested)
*	Statistics and Machine Learning Toolbox (Version 12.3 tested)

## 2. Data Processing
### 2.0 Data Preparation
The input data should include the measured CIEXYZ values of the 24-patch ColorChecker for the test target and the endoscopy device. The input files are text files in the comma-separated value (CSV) format. Each input file contains 25 rows of the CIE X, Y, and Z values separated by commas. Lines #1~#24 describe the 24 color patches, and Line #25 describes the reference white. Follow the patch order defined by the ColorChecker.

The following is a sample input file:
```
19.1,16.2,8.9
76.6,72.2,58.4
40.6,41.7,80.7
20.1,26.4,9.8
61.3,58.6,99.5
80.7,95,106.7
66.8,58.4,12.1
24.5,20.3,77.6
47.3,29.1,24.3
5.5,3.3,12.1
66.5,85.1,19.3
74.1,74.9,11.2
10.8,6.8,46
29.3,46.6,15.1
29.2,16.4,5.3
95,107.4,14.7
56.4,38,69.2
44.2,51.1,100.4
167.3,174.2,193.4
124.8,130.3,145.4
88.1,92.2,103.4
45.2,47.3,52.6
13.6,14.2,17.4
0.6,0.6,0.8
167.3,174.2,193.4
```
*Figure 2.0: CIEXYZ input file. A sample input file of the CIEXYZ data measured from the device.*

### 2.1 Import Data
Provide filenames for loading data files that contain the measured CIEXYZ values of the 24-patch ColorChecker for the test target (e.g., `XYZ_Reference.csv`) and the subject device (e.g., `XYZ_Subject.csv`).

### 2.2 CIEXYZ Data
The input data for the ground truth (`Ref_X`, `Ref_Y`, and `Ref_Z`) and endoscopy device (`Sub_X`, `Sub_Y`, and `Sub_Z`) are combined into `CIEXYZ_Data` as a 25x6 table for inspection:

![CIEXYZ_Data](https://user-images.githubusercontent.com/45103074/212110144-2e723b08-7329-4885-82e9-8328bd7cb4b7.png)

*Figure 2.2: Verification of the input data. The CIEXYZ data of the ground truth (first three columns) and the device output (last three columns).*

### 2.3 CIELAB Data
The CIEXYZ data are converted into CIELAB data. Show the CIELAB data (`CIELAB_Data`) for the user to examine:

![CIELAB_Data](https://user-images.githubusercontent.com/45103074/212108451-dcf69076-0c6f-4f26-a860-4a5050392dfe.png)

*Figure 2.3: Verification of the color conversion. The converted CIELAB data of the ground truth (first three columns) and the device output (last three columns).*

## 3. Data Analysis

### 3.1 Visual simulation of the test target and sample scenes

#### 3.1.1 Visualize CIEXYZ data
The following charts show the simulated visual results when using D65 as the reference white. Use these charts to check excessive color shift caused by the light source and/or the device.

![xyzx2](https://user-images.githubusercontent.com/45103074/212109243-153568a8-aa52-4b69-877f-a05d6572e894.png)
*Figure 3.1.1: Visual verification of the CIEXYZ data. The charts show the simulated visual results when using D65 as the reference white. Use these charts to check excessive color shift caused by the light source and/or the device.*

#### 3.1.2 Visualize CIELAB data
The following charts show the simulated visual results when using the provided reference white. Use these charts to assess how the device would reproduce the ColorChecker.
 
![labx2](https://user-images.githubusercontent.com/45103074/212101296-06020a3f-6dc1-4144-a121-e5dfb2629699.png)
*Figure 3.1.2: Visual verification of the CIELAB data. The charts show the simulated visual results when using the provided reference white. Use these charts to assess how the device would reproduce the test target.*

#### 3.1.3 Visualize endoscopic scene
The following charts show the simulated visual results of an endoscopic scence. Use the default polyp sample or provide a different image according to the intended use.

![samplex2](https://user-images.githubusercontent.com/45103074/212101187-8bdd593a-ded8-42d8-a28a-ab9ab832f60b.png)
*Figure 3.1.3: Visualize endoscopic scenes. The charts show the simulated visual results of an endoscopic scene. Use the default polyp sample or provide a different image according to the intended use.*

### 3.2 Absolute color errors in comparison with the ground truth
The left chart shows the per-patch color difference between the endoscopy device and the ground truth. The right chart shows the boxplot. Statistics (mean, std, min, median, and max) are provided in the titles.
 
![bar_de](https://user-images.githubusercontent.com/45103074/212133919-e8808080-d472-4244-8bc5-3425f2d7f881.png)

*Figure 3.2: Absolute color errors in comparison with the ground truth. The left chart shows the per-patch color difference between the subject device and the ground truth. The right chart shows the box plot. The statistics (mean, standard deviation, minimum, median, and maximum) are provided in the titles.*

### 3.3 Preservation of the patch order in lightness, hue, and chroma

#### 3.3.1 Order in lightness, chroma, and hue - 1D view
The following charts show the patch order in lightness, chroma, and hue (1D view). The top row is the reference, the bottom row is the device output, and each line connects the same patch. Use these charts to identify any out-of-order patches and assess concordance and monotonicity.

![lch_1d](https://user-images.githubusercontent.com/45103074/212101012-e9750cec-375a-4b4e-a665-82414df1f4fa.png)
*Figure 3.3.1: Preservation of the patch order in lightness/hue/chroma -- monotonicity. The charts show the patch order in lightness, chroma, and hue. In each chart, the top row is the reference, the bottom row is the device output, and each line connects the same patch. The top chart shows the patch order in lightness for gray patches (#19-#24). The remaining three charts show the chromatic patches in the lightness, hue, and chroma order. Use these charts to identify any out-of-order patches. The Kendall Tau-a rank correlation coefficients are included in each plot.*

#### 3.3.2 Order in lightness, chroma, and hue - 2D view
The following charts show the patch order in lightness, chroma, and hue (2D view). Use these charts to identify any out-of-order patches and assess linearity, concordance, and monotonicity.
 
![lch_2d](https://user-images.githubusercontent.com/45103074/212100918-f8b7dc53-a073-4d27-b737-a16697fb3381.png)
*Figure 3.3.2: Preservation of the patch order in lightness/hue/chroma -- linearity. The top two charts show the patch order in lightness for gray patches (#19-#24) and chromatic patches (#1-#18). The lower left chart shows the chromatic patches in the hue order, and the lower right in the chroma order. Linear regression coefficients are included in each plot.*

#### 3.3.3 Three-dimensional color transfer
The following charts show the color transfer from the ground truth (spheres) to the device output (crosses) of all patches in the CIELAB color space. Rotate the 3D plots in Matlab to observe the spatial relationship. 

![labx6](https://user-images.githubusercontent.com/45103074/212100782-49a22622-ca5a-44ee-8185-e7ab3ed1088c.png)
*Figure 3.3.3: Visualization of color transfer. The top four charts show the color transfer of the chromatic patches observed from different angles. The bottom two charts show the color transfer of the gray patches to observe their lightness and chromaticity.*

### 3.4 Preservation of Color Contrast between Patches
The following charts show all datapoints plotted according to their ground truth and device output. The color contrast enhancement (CCE) is defined as follows:

![CCE_eq](https://user-images.githubusercontent.com/45103074/212115014-c9581a94-ea4c-493a-8397-644011e11c07.png)

where *i* and *j* are patch numbers. *ref(i)* and *sub(i)* are reference and subject device output for patch *#i*, respectively. *ΔE* is a function calculating the CIE color difference between two colors based on either the *ΔE<sub>00</sub>*, *ΔE<sub>94</sub>*, or *ΔE<sub>76</sub>* formulas. *CCE=1* (represented by the dotted red line) means that the device reproduces the color contrast perfectly. *CCE>1* means that the device enhances the color contrast, while *CCE<1* means otherwise.

![ccex3](https://user-images.githubusercontent.com/45103074/212100537-126d011c-7825-42b7-9a19-ed7474b0c606.png)
*Figure 3.4: Preservation of color contrast between patches. The CCE values calculated based on the ΔE<sub>00</sub>, ΔE<sub>94</sub> , or ΔE<sub>76</sub> formulas. Each colored cross represents a patch-pair where the horizontal and vertical bars are colored separately according to the patch-pair. The percentage indicates patch-pairs that have CCE>1.*
