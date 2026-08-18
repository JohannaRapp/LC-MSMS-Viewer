# LC-MSMS-Viewer
Can be used to load targeted LC-MS/MS data from Agilent HR-QTOF. The Matlab GUI imports raw data, extracts EICs (extracted ion chromatograms) for targeted metabolites, extracts MS2 fragments and compare them with predicted spectra from CFM-ID 4.0.

# Example Data
You can find example data (including a video showing how to import data) here: ADD ZENODOO link here

# How to use this GUI
1.	Convert raw files (.d) in .mzXML files using MSConvert (Settings: 64-bit, Write index, TPP compatibility, Figure 4).
2.	Install MATLAB (tested for Version R2023b). Install the following toolboxes: “Bioinformatic Toolbox”, “Parallel Computing Toolbox”, “Signal Processing Toolbox”.
3.	Download MATLAB App from GitHub. Unzip downloaded folder.
4.	Download example data from Zenodoo. Unzip downloaded folder.
5.	Move downloaded App folder and example data folder in a working folder (e.g. in the MATLAB folder).
Note: Steps 6-11 are visualized in a video which can be found in the example dataset from Zenodoo (LC_MSMS_Viewer.mp4).
6.	Open MATLAB and select the working folder. Add folders and subfolders to the MATLAB path (right click: “Add To Path” -> “Selected Folders and Subfolders”. Right click on “App_LC_MSMS_Viewer.mlapp” -> Run. Graphical user interface (GUI) will open. 
7.	Press “Load” -> “Load RawData”. Select all .mzXML files from Example Data. Select “hitlist” (you can find an example HitList in the Example data). Enter a filename and destination for saving the imported data.
8.	In the command line you can see the progress of data import. 
9.	Imported data are saved in a Matlab structure with the following name: YYYYMMDD_LC_MSMS_data_SelectedFilename.mat
Data are also exported into two Excel sheets:
YYYYMMDD_LC_MSMS_data_SelectedFilename.xlsx (MS2 fragments, fold-changes, predicted fragments, etc.)
YYYYMMDD_LC_MSMS_data_SelectedFilename_QC.xlsx (Quality Control results)
10.	To visualize EICs and MS1 and MS2 spectra load Matlab structure in the App (Figure 6). Press “Load” -> “Load LC_MSMS_data” -> Select previously saved Matlab structure (YYYYMMDD_LC_MSMS_data_SelectedFilename.mat). Select metabolite-strain pair for visualization in the Table on the left side. To change polarity use slider.
Critical: If GUI is not working ensure that the folder “App_LC_MSMS_Viewer” AND the subfolder “Functions” is added to the matlab path.
11.	Use exported Excel-Sheet with MS2 fragments for further analysis.


Note: The Matlab GUI imports the mzXML files and extracts MS1 and MS2 scans via mzXMLread.m function. Extracted ion chromatograms (EICs) for all metabolites targeted within one polarity are created comparing the m/z-vector of each MS1 scan with the exact m/z-value of the precursor ion. m/z-values and corresponding intensities, which differed less than Δ=0.003 Da from the theoretical m/z-value are used for the EIC. If several matching m/z-values are found, the one with the highest intensity is used for the EIC. The three MS2 scans that are in between the MS1 scans with the highest and second highest intensity for the precursor ion are used for MS2 spectra. Peak Picking in the MS2 spectra is performed with the findpeaks.m function of MATLAB (Height Filter: 400, Prominence: 400). MS1 scan with the highest intensity of the precursor ion is used to analyze the precursor ion purity, mass accuracy (0.003 Da) and intensity (intensity threshold 8000). MS2 fragments matching the precursor m/z and or larger than the precursor m/z are excluded. MS2 fragments with a normalized intensity (normalized on highest peak in the spectrum) smaller 0.05 are excluded.
