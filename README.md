<h1>Air-water flow data analysis software for double-tip phase-detection intrusive probes (https://doi.org/10.5281/zenodo.2448250)</h1>

<p>The software is able to calculate a wide range of air-water flow properties in high-velocity air-water flows typical for hydraulic engineering applications. The software is suitable for two simultaneously sampled signals of a phase-detection intrusive probe including one double-tip probe or two single tip probes.</p>
<p>The software was developed by Dr Stefan Felder (UNSW Sydney) during his PhD project (Felder 2013) and has been used in many scientific publications for post-processing of conductivity probe signals. If you decide to use the software for post-processing of your data, please credit the creator of the work as:</p>
<p>Felder, S. (2013). “Air-Water Flow Properties on Stepped Spillways for Embankment Dams: Aeration, Energy Dissipation and Turbulence on Uniform, Non-Uniform and Pooled Stepped Chutes.” PhD thesis, The University of Queensland.</p>
and provide details on the release version of the software as:
<p>Felder, S. (2018). “StefanFelder/Air-water-flow-data-analysis-software-for-double-tip-phase-detection-intrusive-probes (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.2448251</p>

<p>The software has been already made public as a digital appendix in the PhD thesis of Felder (2013) and is still available for download here: https://espace.library.uq.edu.au/view/UQ:301329. Note that the readme file of the digital appendix is also provided as a PDF document for further documentation on the software (Note that the digital appendix refers to two further data processing tools including for two side-by-side single tip probes and the triple decomposition technique for instationary air-water flows; the present software can be used also for two side-by-side probes and the triple decomposition software will be made available in a separate repository).</p> 
<p>The software is able to analyse two simultaneously sampled raw Voltage signals of double-tip phase-detection intrusive probes (or two single tip probes) including conductivity and fiber-optical probes:</p>

- For double tip conductivity probe signals use “Data_analysis_software_conductivity_probe”
- For double tip fiber-optical probe signals use “Data_analysis_software_fiber-optical_probe”

The software is able to provide information about the following air-water flow properties:
- Void fraction C (-) (both leading and trailing tips)
- Bubble count rate F (Hz) (both leading and trailing tips)
- Interfacial velocity V (m/s) (Not meaningful for two side-by-side probes)
- Turbulence intensity Tu (-) (Not meaningful for two side-by-side probes)
- Bubble/droplet chord time tch (s)
- Bubble/droplet chord length ch (m) (Not meaningful for two side-by-side probes)
- Auto-correlation function Rxx (-)
- Cross-correlation function Rxy (-)
- Maximum cross-correlation coefficient (Rxy)max (-)
- Auto-correlation time scale Txx (s) 
- Cross-correlation time-scale Txy (s) 
- Cluster properties based upon near-wake, constant and percentage criteria
- Inter-particle arrival time (s)


<b>1 Content:</b>
<p>The software is written in Fortran Compaq. It has been successfully compiled with any Intel Fortran compiler including the latest version of Intel Fortran using Parallel Studio XE 2018 updates 2 or 3. The repository contains:</p>

- Full source code for double-tip conductivity probe raw Voltage data
- Full source code for double-tip fiber-optical probe raw Voltage data
- Digital Appendix of Felder (2013) with additional documentation
- “Supplementary input data files for software”, which are needed to run the software comprising three input text files:
	- parameter.txt	: This file contains important input parameters (see Table DA-4 in the Digital Appendix of Felder (2013))
	- binary.txt :  The recorded input files are in binary format and must be of 8 character lengths. (Details regarding the required binary file format are available in the source code: lines 359 to 384.)
	- positions.txt :  The measurement locations must be of 4 character lengths and will appear in the summary result files.
- The folder “example raw data” containing some example binary raw data files of the thesis of Felder (2013). (The three input text files match the raw data). Note further that the raw data were acquired with a LabVIEW data acquisition software, which was also developed by Felder (2013). (Further details on this acquisition software can be found in Felder (2013) and the Digital Appendix file.)

<b>2 How to run the code:</b>
- Compile the source code with a Fortran compiler (it works with any form of Intel or Compaq Fortran compilers) to create an executable file.
- Copy this file into the same folder as your binary raw data together with the three supplementary input text files “parameter”, “binary” and “positions”.
- Run the executable file and follow the prompts in the software. The software will analyse all raw data files and will produce a range of result text files containing the wide range of air-water flow properties.

<b>3 Contact:</b>
<p>For feedback, questions and recommendations, please use the issue-section or contact the author via Email:
Stefan Felder, Senior Lecturer, Water Research Laboratory, School of Civil and Environmental Engineering, UNSW Sydney, Australia; Email: s.felder@unsw.edu.au </p>

<b>4. Selected References:</b>
<p>The software has been used in a number of studies and research projects at The University of Queensland and UNSW Sydney. Some selected references are:</p>

Chachereau, Y. and Chanson, H. (2011). "Bubbly Flow Measurements in Hydraulic Jumps with Small Inflow Froude Numbers." International Journal of Multiphase Flow, Vol. 37, No. 6, pp. 555-564.

Chanson, H. and Chachereau, Y. (2013). "Scale Effects Affecting Two-Phase Flow Properties in Hydraulic Jump with Small Inflow Froude Number." Experimental Thermal and Fluid Science, Vol. 45, pp. 234-242. 

Felder, S. (2013). “Air-Water Flow Properties on Stepped Spillways for Embankment Dams: Aeration, Energy Dissipation and Turbulence on Uniform, Non-Uniform and Pooled Stepped Chutes.” PhD thesis, The University of Queensland.

Felder, S. and Chanson, H. (2011). "Energy dissipation down a stepped spillway with nonuniform step heights', Journal of Hydraulic Engineering, ASCE, Vol. 137, pp. 1543–1548.

Felder, S. and Chanson, H. (2011). "Air-water flow properties in step cavity down a stepped chute", International Journal of Multiphase Flow, Vol. 37, pp. 732–745.

Felder, S. and Chanson, H. (2013). "Aeration, flow instabilities, and residual energy on pooled stepped spillways of embankment dams', Journal of Irrigation and Drainage Engineering, Vol. 139, pp. 880–887.

Felder, S. and Chanson, H. (2015). "Phase-detection probe measurements in high-velocity free-surface flows including a discussion of key sampling parameters", Experimental Thermal and Fluid Science, vol. 61, pp. 66–78.

Felder, S. and Chanson, H. (2016). "Air–water flow characteristics in high-velocity free-surface flows with 50% void fraction", International Journal of Multiphase Flow, vol. 85, pp. 186–195.

Felder, S. and Chanson, H. (2017). "Scale effects in microscopic air-water flow properties in high-velocity free-surface flows", Experimental Thermal and Fluid Science, vol. 83, pp. 19–36.

Felder, S. and Pfister, M. (2017). "Comparative analyses of phase-detective intrusive probes in high-velocity air–water flows" International Journal of Multiphase Flow, Vol. 90, pp. 88–101.

Felder, S. and Chanson, H. (2018). "Air-water flow patterns of hydraulic jumps on uniform beds macroroughness", Journal of Hydraulic Engineering, ASCE, Vol. 144, Paper 0001402.

Guenther, P., Felder, S. and Chanson, H. (2013). "Flow aeration, cavity processes and energy dissipation on flat and pooled stepped spillways for embankments", Environmental Fluid Mechanics, Vol. 13, pp. 503–525.

Severi, A. (2018). "Aeration Performance and Flow Resistance in High-Velocity Flows over Moderately Sloped Spillways with micro-rough bed", PhD thesis, UNSW Sydney, Australia.

Wang, H. (2014). "Turbulence and Air Entrainment in Hydraulic Jumps." Ph.D. thesis, School of Civil Engineering, The University of Queensland, Brisbane, Australia.

Zhang, G., Wang, H. and CHANSON, H. (2013). "Turbulence and Aeration in Hydraulic Jumps: Free-Surface Fluctuation and Integral Turbulent Scale Measurements." Environmental Fluid Mechanics, Vol. 13, No. 2, pp. 189–204.

