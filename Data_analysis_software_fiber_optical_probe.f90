!  dataanalysis.f90 

!****************************************************************************
!
!  PROGRAM: dataanalysis
!
!  Version 2
!  by Stefan Felder - August 2009 - 
!
!****************************************************************************

program dataanalysis

implicit none

	!Parameters

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real(kind=DBL) 
															!variables
	integer, parameter :: short = SELECTED_INT_KIND(3)

	!General variables
	!parameters to be loaded from file
	REAL :: delta_x
	INTEGER :: frequency
	INTEGER :: Correlation_steps
	INTEGER :: sample_duration
	INTEGER :: Number_devices
	INTEGER :: segment_numbers
	REAL :: PDF_V_segment_size
	INTEGER :: max_no_bubbles
	REAL :: threshold_C
	REAL :: PDF_Court_length_segment_size
	REAL :: PDF_Court_time_segment_size
	REAL :: PDF_Court_min_length
	REAL :: PDF_Court_max_length
	REAL :: PDF_Court_min_time
	REAL :: PDF_Court_max_time
	INTEGER :: cluster_or_not
	REAL :: cluster_factor
	REAL :: cluster_factor_1

	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Data_array
	CHARACTER(len=9) :: name_parameters			!name of file with parameters to be opened
	CHARACTER(len=9) :: name					!name of file with positions to be opened
	CHARACTER(len=6) :: name_binary
	CHARACTER(len=8) :: location_binary
	CHARACTER(len=4) :: position
	CHARACTER(len=8) :: position_binary				
	CHARACTER(len=4) :: ending = '.txt'
	CHARACTER(len=4) :: ending_binary = '.dat'
	CHARACTER(len=1) :: sample = '_'
	INTEGER :: number_rows = 0					!number of rows to read
	INTEGER :: status = 0							!I/O status	
	REAL(kind=DBL) :: value	= 0
	INTEGER :: number_positions	= 0
	CHARACTER(len=8), ALLOCATABLE, DIMENSION(:) :: Position_array_binary
	CHARACTER(len=4), ALLOCATABLE, DIMENSION(:) :: Position_array
	INTEGER :: number_positions_binary  =0 
	CHARACTER(len=4) :: location					
	INTEGER :: number_columns = 0					
	INTEGER :: i,j,k
	INTEGER (KIND=short) :: x, y
	character(len=8) :: header						!header gives info about size of binary file in labview code

	!data arrays for annalyses

	! Variables for correlation-analyses
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Autocorrelation			!results_file
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: Crosscorrelation			!results_file

	!Variables for Propability density distribution of Voltage signal for both probe tips
	REAL(kind=DBL) :: PDF_segments	= 0							!Number of segments for PDF
	REAL(kind=DBL), ALLOCATABLE, DIMENSION (: , :) :: PDF_V
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:) :: threshold			!results_file


	!Variables for Correlation-scale analyses
	REAL(kind=DBL), DIMENSION(14) :: CorScale_results				!summary file with Correl-
																	!ation properties 
	REAL(kind=DBL) :: velocity = 3									!default value

	!Variables for basic property analyses (air-water-interfaces, void fraction, frequency)
	INTEGER, ALLOCATABLE, DIMENSION(:,:) :: interfaces !interfaces
	REAL(kind=DBL), ALLOCATABLE, DIMENSION(:,:) :: basic_results	!C and F

	!Variables for chord time and size distributions
	REAL(kind=DBL) :: PDF_number_chord = 0
	REAL(kind=DBL) :: PDF_number_time = 0
	REAL(kind=DBL), ALLOCATABLE, DIMENSION (: , :) :: PDF_chord_size
	REAL(kind=DBL), ALLOCATABLE, DIMENSION (: , :) :: PDF_chord_time


	!Variables for cluster analyses
	REAL(kind=DBL) :: PDF_cluster_distr (20,2)
	REAL(kind=DBL) :: cluster_results (9, 1)

	!Variables for cluster analyses constant cluster_criterion
	REAL(kind=DBL) :: PDF_cluster_distr_1 (20,2)
	REAL(kind=DBL) :: cluster_results_1 (9, 1)

	!Variables for cluster analyses with wake criterion
	REAL(kind=DBL) :: PDF_cluster_distr_2 (20,2)
	REAL(kind=DBL) :: cluster_results_2 (9, 1)

	! variables for particle grouping analyses
	REAL(kind=DBL) :: PDF_particle_group (101, 11)
	REAL(kind=DBL) :: PDF_particle_group_time (101, 11)

	!Open file to get parameters from parameter file

	!Get the filename of the list of parameters for the analyses
	WRITE (*,*) 'Please input file name with parameters:'
	READ (*,*) name_parameters
	WRITE (*,501) name_parameters
	501 format (' ', 'The filename is: ', A '.txt')

	OPEN (UNIT = 21, FILE = name_parameters // ending, STATUS = 'OLD', ACTION ='READ', IOSTAT=status)
IF (status == 0) THEN

	READ (21, 100) delta_x
	100 FORMAT (//, 1F5.2)

	READ(21, 101) frequency
	101 FORMAT (/, 1I7)

	READ(21, 102) Correlation_steps
	102 FORMAT (/, 1I5)

	READ(21, 103) sample_duration
	103 FORMAT (//, 1I5)

	READ(21, 104) Number_devices
	104 FORMAT (//, 1I1)

	READ(21, 105) segment_numbers
	105 FORMAT (/, 1I4)

	READ(21, 106) PDF_V_segment_size
	106 FORMAT (//, 1F5.2)
			
	READ(21, 107) max_no_bubbles
	107 FORMAT (/, 1I6)

	READ(21, 108) threshold_C
	108 FORMAT (//, 1F5.2)

	READ(21, 109) PDF_Court_length_segment_size
	109 FORMAT (/, 1F5.2)

	READ(21, 110) PDF_Court_time_segment_size
	110 FORMAT (/, 1F5.2)

	READ(21, 111) PDF_Court_min_length
	111 FORMAT (/, 1F5.2)

	READ(21, 112) PDF_Court_max_length
	112 FORMAT (/, 1F5.2)

	READ(21, 113) PDF_Court_min_time
	113 FORMAT (/, 1F5.2)

	READ(21, 114) PDF_Court_max_time
	114 FORMAT (/, 1F5.2)

	READ(21, 115) 	cluster_or_not
	115 FORMAT (/, 1I1)

	READ(21, 116) cluster_factor
	116 FORMAT (//, 1F5.2)

	READ(21, 117) cluster_factor_1
	117 FORMAT (/, 1F5.2)

END IF
    CLOSE ( UNIT=21)


	!Allocate arrays for results with parameters
	ALLOCATE (Autocorrelation(2, Correlation_steps), STAT =status)
	ALLOCATE (Crosscorrelation(2, Correlation_steps), STAT =status)		

	ALLOCATE (threshold(Number_devices),STAT=status)

	ALLOCATE (interfaces(Number_devices*2, max_no_bubbles), STAT =status)
	ALLOCATE (basic_results(Number_devices*2, 1), STAT =status)


	!Write Headings to result-files. Note: This needs to be done before the main loops-starts because the 
	!heading should just appear once in the sumary result files (Tu_V, basics and chordlength/-times
	OPEN (UNIT = 50, FILE =  'Tu_V_summary' // ending, STATUS = 'REPLACE', &
		 ACTION ='WRITE', IOSTAT=status)
		WRITE (50, 150) 'Position', ' T0.5 ', ' T(Rxx=0) ', ' T(Min(Rxx)) ', ' Min(Rxx) ', ' Txx ', ' (Rxz)max ' , ' T(Rxz)max ', ' T(Rxz=0)&
			' , ' T(Min(Rxz)) ', ' Min(Rxz) ', ' T(0.5(Rxz)max) ' , ' Txz ' , ' Tu ' , ' V '
		150 FORMAT (A9, T11, A5, T21, A9, T31, A12, T43, A9, T52, A4, T62, A9, T72, A10, T83, A9, T94,&
			   A12, T106, A9, T115, A15, T131, A4, T141, A3, T151, A2)

	OPEN (UNIT = 52, FILE =  'basics_summary'  // ending , STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (52, 152) 'Position', 'C-leading', 'F-leading', 'C-trailing', 'F-trailing'
		152 FORMAT (A9, T11, A10, T22, A10, T33, A11, T45, A11, T57)

	OPEN (UNIT = 54, FILE =  'chordlength_summary'  // ending  , STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)

	OPEN (UNIT = 56, FILE =  'chordtime_summary'  // ending  , STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)

	IF(	cluster_or_not == 1) THEN

	OPEN (UNIT = 58, FILE =  'cluster_summary'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)
	WRITE (58, 1153) 'Position', 'C-leading', 'No_in_cluster', '%_in_cluster', 'No_of_clusters', 'cluster_per_s', 'ave_no_per_cluster', 'ave_size_cluster', 'ave_rate_cluster_chord_size', 'ratio_lead_ave'
		1153 FORMAT (A9, T11, A10, T22, A14, T37, A14, T53, A16, T70, A15, T86, A20, T107, A18, T126, A30, T160, A16)

	OPEN (UNIT = 59, FILE =  'cluster_no_PDF'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)

	OPEN (UNIT = 78, FILE =  'cluster_summary_constant_criterion_1mm'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)
	WRITE (78, 1173) 'Position', 'C-leading', 'No_in_cluster', '%_in_cluster', 'No_of_clusters', 'cluster_per_s', 'ave_no_per_cluster', 'ave_size_cluster', 'ave_rate_cluster_chord_size', 'ratio_lead_ave'
		1173 FORMAT (A9, T11, A10, T22, A14, T37, A14, T53, A16, T70, A15, T86, A20, T107, A18, T126, A30, T160, A16)

	OPEN (UNIT = 79, FILE =  'cluster_no_PDF_constant_criterion_1mm'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)

	OPEN (UNIT = 98, FILE =  'cluster_summary_wake_criterion'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)
	WRITE (98, 1193) 'Position', 'C-leading', 'No_in_cluster', '%_in_cluster', 'No_of_clusters', 'cluster_per_s', 'ave_no_per_cluster', 'ave_size_cluster', 'ave_rate_cluster_chord_size', 'ratio_lead_ave'
		1193 FORMAT (A9, T11, A10, T22, A14, T37, A14, T53, A16, T70, A15, T86, A20, T107, A18, T126, A30, T160, A16)

	OPEN (UNIT = 99, FILE =  'cluster_no_PDF_wake_criterion'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)

	OPEN (UNIT = 854, FILE =  'interparticle_arrival_time_summary'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)

	ELSE IF (	cluster_or_not == 2) THEN

	OPEN (UNIT = 558, FILE =  'cluster_summary_time'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)
	WRITE (558, 5153) 'Position', 'C-leading', 'No_in_cluster', '%_in_cluster', 'No_of_clusters', 'cluster_per_s', 'ave_no_per_cluster', 'ave_size_cluster', 'ave_rate_cluster_chord_size', 'ratio_lead_ave'
		5153 FORMAT (A9, T11, A10, T22, A14, T37, A14, T53, A16, T70, A15, T86, A20, T107, A18, T126, A30, T160, A16)

	OPEN (UNIT = 559, FILE =  'cluster_no_PDF_time'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)

	OPEN (UNIT = 578, FILE =  'cluster_summary_constant_criterion_1mm_time'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)
	WRITE (578, 5173) 'Position', 'C-leading', 'No_in_cluster', '%_in_cluster', 'No_of_clusters', 'cluster_per_s', 'ave_no_per_cluster', 'ave_size_cluster', 'ave_rate_cluster_chord_size', 'ratio_lead_ave'
		5173 FORMAT (A9, T11, A10, T22, A14, T37, A14, T53, A16, T70, A15, T86, A20, T107, A18, T126, A30, T160, A16)

	OPEN (UNIT = 579, FILE =  'cluster_no_PDF_constant_criterion_1mm_time'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)

	OPEN (UNIT = 598, FILE =  'cluster_summary_wake_criterion_time'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)
	WRITE (598, 5193) 'Position', 'C-leading', 'No_in_cluster', '%_in_cluster', 'No_of_clusters', 'cluster_per_s', 'ave_no_per_cluster', 'ave_size_cluster', 'ave_rate_cluster_chord_size', 'ratio_lead_ave'
		5193 FORMAT (A9, T11, A10, T22, A14, T37, A14, T53, A16, T70, A15, T86, A20, T107, A18, T126, A30, T160, A16)

	OPEN (UNIT = 599, FILE =  'cluster_no_PDF_wake_criterion_time'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)

	OPEN (UNIT = 954, FILE =  'interparticle_arrival_time_summary_time'  // ending  , STATUS = 'REPLACE', & 
	ACTION ='WRITE', IOSTAT=status)


	END IF

	number_columns = Number_devices +1				!Columns to be allocated in the arrays 
													!depending on the number of devices used
	number_rows = frequency * sample_duration		!Rows to be allocated in the data array

	PDF_segments = REAL(7.0)/PDF_V_segment_size		!Number of segments for PDF (segments between 0 and 5V)
	
	PDF_number_chord = 1+(PDF_Court_max_length-PDF_Court_min_length)/PDF_Court_length_segment_size

	PDF_number_time = 1+(PDF_Court_max_time-PDF_Court_min_time)/PDF_Court_time_segment_size


	!get the list with the binary file list
	WRITE (*,*) 'Please input file name with the binary file names:'
	READ (*,*) name_binary
	WRITE (*,510) name_binary
	510 format (' ', 'The name with the binary file names is  ', A '.txt')

	!Get the filename of the list of positions for the analyses
	WRITE (*,*) 'Please input file name with locations:'
	READ (*,*) name
	WRITE (*,500) name
	500 format (' ', 'The filename is: ', A '.txt')

	OPEN (UNIT = 20, FILE = name_binary // ending, STATUS = 'OLD', ACTION ='READ', IOSTAT=status)
	IF (status == 0) THEN

	!open was ok. Read values to find out how many lines are in the file.
		DO 
			READ (20, *, IOSTAT=status) location_binary	!Get next value
			IF (status /=0) EXIT						!EXIT if not valid
			number_positions_binary = number_positions_binary + 1		!Valid: increase count
		END DO 
	END IF

	OPEN (UNIT = 22, FILE = name // ending, STATUS = 'OLD', ACTION ='READ', IOSTAT=status)
	IF (status == 0) THEN

	!open was ok. Read values to find out how many lines are in the file.
		DO 
			READ (22, *, IOSTAT=status) location	!Get next value
			IF (status /=0) EXIT						!EXIT if not valid
			number_positions = number_positions + 1		!Valid: increase count	
		END DO 
	END IF

	IF (number_positions_binary == number_positions) THEN

	WRITE (*,520) number_positions
	520 format (' ', 'The file contains: ', I ' vertical locations in the cross section')

		ALLOCATE ( Position_array_binary(number_positions), STAT=status)	!allocate memory
		ALLOCATE ( Position_array(number_positions), STAT=status)

		! if the allocation was successful, rewind the file and read in the data!
		allocate_ok: IF ( status == 0) THEN
		REWIND ( UNIT=20)						!Rewind file	
		READ (20,*) (Position_array_binary(i), i=1, number_positions) 
		CLOSE ( UNIT=20)! Close File 
		END IF allocate_ok	

		allocate_ok2:IF ( status == 0) THEN
		REWIND ( UNIT=22)						!Rewind file	
		READ (22,*) (Position_array(i), i=1, number_positions) 
		CLOSE ( UNIT=22)! Close File 
		END IF allocate_ok2	
	CLOSE ( UNIT=22)

	ELSE 
	WRITE (*,*) 'The number of positions in the two input files is different - please restart the program!'

	END IF	


	DO i = 1, number_positions      !start of main loop in program!!!!
		position = Position_array(i)
		position_binary = Position_array_binary(i)			!Get the filename and echo it back to the user
		WRITE (*,10000) i, position
		10000 format (' ', 'The position ' I ' is: ' , A)

	! Open the position file, and check for errors on open.
	OPEN (UNIT = 30, FILE = position_binary // ending_binary, form='binary', status = 'old', action = 'read', iostat = status)

	IF (status == 0) THEN

		! ALLOCATE MEMORY		
		ALLOCATE ( Data_array(number_columns, number_rows), STAT=status)	!allocate memory
		Data_array = 0

		! if the allocation was successful, rewind the file and read in the data!
		allocate_ok3: IF ( status == 0) THEN

		READ (30) header

		Do j = 1, number_rows
		READ (30) x, y
		Data_array (1,j) = (j-1)/real(frequency)
		Data_array (2,j) = x*5./32768.
		Data_array (3,j) = y*5./32768.

		END DO

		END IF allocate_ok3				
		END IF 
		
	! Close File 
	CLOSE ( UNIT=30)	

	!use the data arry and do auto- and cross -correlation analyses
	CALL Correlation(Data_array, number_columns, number_rows, Correlation_steps, & 
			         Autocorrelation, Crosscorrelation, frequency, segment_numbers)

	!calculate the characteristic scales, Tu and Interfacial Velocity
	CALL Correlation_scales (Autocorrelation, Crosscorrelation, Correlation_steps, & 
			 CorScale_results, frequency, delta_x, velocity) 

	!Calculate the Probability distribution function of the voltage signals and store to array PDF_V
		ALLOCATE (PDF_V(number_columns, Int(PDF_segments)), STAT=status)	!allocate memory for PDF_V
			PDF_V = 0
	CALL Probability_V (Data_array, number_columns, number_rows, PDF_V_segment_size, &
						PDF_segments, PDF_V, threshold, threshold_C)

	CALL Basic_properties (Data_array, number_columns, number_rows, threshold,&
						 Interfaces, basic_results, max_no_bubbles, frequency)

		ALLOCATE (PDF_chord_size(INT(PDF_number_chord) , 3), STAT=status)
			PDF_chord_size = 0.0
	CALL chord_length (Interfaces, number_columns, max_no_bubbles, velocity, basic_results, &
					   PDF_chord_size, PDF_number_chord, sample_duration, &
					   frequency, PDF_Court_max_length, PDF_Court_min_length, PDF_Court_length_segment_size)

		ALLOCATE (PDF_chord_time(INT(PDF_number_time) , 3), STAT=status)
			PDF_chord_time = 0.0
	CALL chord_time (Interfaces, number_columns, max_no_bubbles, basic_results,PDF_chord_time, PDF_number_time, sample_duration, &
					   frequency,  PDF_Court_max_time, PDF_Court_min_time, PDF_Court_time_segment_size)

	IF(	cluster_or_not == 1) THEN

	CALL Cluster_analysis (interfaces, number_columns, max_no_bubbles, basic_results, velocity, sample_duration, &
					   frequency, cluster_factor, PDF_cluster_distr, cluster_results)

	CALL Cluster_analysis_1 (interfaces, number_columns, max_no_bubbles, basic_results, velocity, sample_duration, &
					   frequency, cluster_factor_1, PDF_cluster_distr_1, cluster_results_1)

	CALL Cluster_analysis_wake (interfaces, number_columns, max_no_bubbles, basic_results, velocity, sample_duration, &
					   frequency, PDF_cluster_distr_2, cluster_results_2)

	CALL Particle_grouping (interfaces, number_columns, max_no_bubbles, basic_results, &
					   sample_duration, frequency, velocity, PDF_particle_group )

	ELSE IF(	cluster_or_not == 2) THEN

	CALL Cluster_analysis_time (interfaces, number_columns, max_no_bubbles, basic_results,  sample_duration, &
					   frequency, cluster_factor, PDF_cluster_distr, cluster_results)

	CALL Cluster_analysis_1_time (interfaces, number_columns, max_no_bubbles, basic_results, sample_duration, &
					   frequency, cluster_factor_1, PDF_cluster_distr_1, cluster_results_1)

	CALL Cluster_analysis_wake_time (interfaces, number_columns, max_no_bubbles, basic_results, sample_duration, &
					   frequency, PDF_cluster_distr_2, cluster_results_2)

	CALL Particle_grouping_time (interfaces, number_columns, max_no_bubbles, basic_results, &
					   sample_duration, frequency, PDF_particle_group_time )

	END IF



	!Write autocorrelation data to file: Autocorrelation-results
	OPEN (UNIT = 60, FILE = position // sample // 'autocor' //  ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (60, 160) 'Time', 'Rxx'
	WRITE (60, 161) ((Autocorrelation (j,k), j=1,2), k=1, Correlation_steps)
	160 FORMAT (A4, T15, A3)
	161 FORMAT (F10.5,T15, F10.7)

	!Write crosscorrelation data to file: Crosscorrelation-results
	OPEN (UNIT = 62, FILE = position //sample //'crosscor'  // ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (62, 162) 'Time', 'Rxy'
	WRITE (62, 163) ((Crosscorrelation (j,k), j=1,2), k=1, Correlation_steps)
	162 FORMAT (A4, T18, A3)
	163 FORMAT (F14.9, T16, F14.9)

	!Write PDF_V distribution data to file: PDF_V-results
	OPEN (UNIT = 64, FILE = position // sample //'PDF_V'  // ending, STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (64,164) 'Bin-size', 'V-leading', 'V-trailing'
	WRITE (64, 165) ((PDF_V (j,k), j=1,3), k=1, INT(PDF_segments))
	164 FORMAT (A8, T13, A9, T25, A10)
	165 FORMAT (F10.5, T13, F10.7, T25, F10.7)
		
	!Write data to file: Turbulence and velocity results
	WRITE (50, 151) Position, (CorScale_results )	
	151 FORMAT (A4, T11, F9.6, T21, F9.6, T31, F9.6, T43, F9.6, T52, F9.6, T62, F9.6, T72, &
			   F9.6, T83, F9.6, T94, F9.6, T106, F9.6, T115, F9.6, T131, F9.6, T141, F9.6, T151, F9.6)
						
	OPEN (UNIT = 66, FILE = position //sample // 'interfaces'  // ending , STATUS = 'REPLACE', & 
		 ACTION ='WRITE', IOSTAT=status)
	WRITE (66,166) 'W->A leading', 'A->W leading', 'W->A trailing', 'A->W trailing'
	WRITE (66, 167) ((interfaces (j,k), j=1,((number_columns-1)*2)), k=1, 10000)
	166 FORMAT (A12, T15, A12, T36, A13, T52, A13)
	167 FORMAT (I, T15, I, T36, I, T52, I)

	!Write to file basic_summary
	WRITE (52, 153) position, ((basic_results (j,k), j=1,(number_columns-1)*2), k=1,1)
	153 FORMAT (A4, T9, F10.5, T20, F10.5, T31, F10.5, T44, F10.5, T55)

	!Write to file chord length summary
	WRITE (54, 155) position, ((PDF_chord_size (j,k), j=1,Int(PDF_number_chord)), k=1,1)
	WRITE (54, 155) 'f(a)', ((PDF_chord_size (j,k), j=1,Int(PDF_number_chord)), k=2,2)
	WRITE (54, 155) 'f(w)',((PDF_chord_size (j,k), j=1,Int(PDF_number_chord)), k=3,3)
	WRITE (54,'')
	155 FORMAT (A4, 1000F10.6, X)

	!Write to file chord time summary
	WRITE (56, 157) position, ((PDF_chord_time (j,k), j=1,Int(PDF_number_time)), k=1,1)
	WRITE (56, 157) 'f(a)',((PDF_chord_time (j,k), j=1,Int(PDF_number_time)), k=2,2)
	WRITE (56, 157) 'f(w)',((PDF_chord_time (j,k), j=1,Int(PDF_number_time)), k=3,3)
	WRITE (56,'')
	157 FORMAT (A4, 1000F10.6, X)

	IF(	cluster_or_not == 1) THEN

	!Write to file cluster_summary
	IF ((basic_results (1,1) < 0.3) .OR. (basic_results (1,1) > 0.7 )) THEN
	WRITE (58, 1154) position, ((cluster_results (j,k), j=1,9), k=1,1)
	1154 FORMAT (A4, T10, F10.5, T21, F12.5, T37, F10.5, T53, F10.5, T70, F10.5, T86, F10.5, T107, F10.5, T126, F10.5, T160, F10.5)
	ELSE
	WRITE (58, 1155) position // ': 0.3 < C < 0.7   - no cluster analyses!'
	1155 FORMAT (A44)
	END IF

	!Write to file cluster number PDF
	if (basic_results (1,1) < 0.3) THEN
	WRITE (59, 1158) 'location: ' // position, ((PDF_cluster_distr (j,k), j=1,20), k=1,1)
	WRITE (59, 1158) 'bubble_cluster',((PDF_cluster_distr (j,k), j=1,20), k=2,2)
	WRITE (59,'')
	1158 FORMAT (A15, 20F10.5, X)
	ELSE IF (basic_results (1,1) > 0.7) THEN
	WRITE (59, 1159) 'location: ' //position, ((PDF_cluster_distr (j,k), j=1,20), k=1,1)
	WRITE (59, 1159) 'droplet_cluster',((PDF_cluster_distr (j,k), j=1,20), k=2,2)
	WRITE (59,'')
	1159 FORMAT (A16, 20F10.5, X)
	ELSE
	WRITE (59, 1169) 'location: ' //position// ': 0.3 < C < 0.7   - no cluster analyses!'
	1169 FORMAT (A55)
	WRITE (59,'')
	END IF

		!Write to file cluster_summary constant criterion 1mm
	IF ((basic_results (1,1) < 0.3) .OR. (basic_results (1,1) > 0.7 )) THEN
	WRITE (78, 1174) position, ((cluster_results_1 (j,k), j=1,9), k=1,1)
	1174 FORMAT (A4, T10, F10.5, T21, F12.5, T37, F10.5, T53, F10.5, T70, F10.5, T86, F10.5, T107, F10.5, T126, F10.5, T160, F10.5)
	ELSE
	WRITE (78, 1175) position // ': 0.3 < C < 0.7   - no cluster analyses!'
	1175 FORMAT (A44)
	END IF

	!Write to file cluster number PDF constant criterion 1mm
	if (basic_results (1,1) < 0.3) THEN
	WRITE (79, 1178) 'location: ' // position, ((PDF_cluster_distr_1 (j,k), j=1,20), k=1,1)
	WRITE (79, 1178) 'bubble_cluster',((PDF_cluster_distr_1 (j,k), j=1,20), k=2,2)
	WRITE (79,'')
	1178 FORMAT (A15, 20F10.5, X)
	ELSE IF (basic_results (1,1) > 0.7) THEN
	WRITE (79, 1179) 'location: ' //position, ((PDF_cluster_distr_1 (j,k), j=1,20), k=1,1)
	WRITE (79, 1179) 'droplet_cluster',((PDF_cluster_distr_1 (j,k), j=1,20), k=2,2)
	WRITE (79,'')
	1179 FORMAT (A16, 20F10.5, X)
	ELSE
	WRITE (79, 1189) 'location: ' //position// ': 0.3 < C < 0.7   - no cluster analyses!'
	1189 FORMAT (A55)
	WRITE (79,'')
	END IF

	!Write to file cluster_summary wake criterion
	IF ((basic_results (1,1) < 0.3) .OR. (basic_results (1,1) > 0.7 )) THEN
	WRITE (98, 1194) position, ((cluster_results_2 (j,k), j=1,9), k=1,1)
	1194 FORMAT (A4, T10, F10.5, T25, F12.5, T37, F10.5, T53, F10.5, T70, F10.5, T86, F10.5, T107, F10.5, T126, F10.5, T160, F10.5)
	ELSE
	WRITE (98, 1195) position // ': 0.3 < C < 0.7   - no cluster analyses!'
	1195 FORMAT (A44)
	END IF

	!Write to file cluster number PDF wake criterion
	if (basic_results (1,1) < 0.3) THEN
	WRITE (99, 1198) 'location: ' // position, ((PDF_cluster_distr_2 (j,k), j=1,20), k=1,1)
	WRITE (99, 1198) 'bubble_cluster',((PDF_cluster_distr_2 (j,k), j=1,20), k=2,2)
	WRITE (99,'')
	1198 FORMAT (A15, 20F10.5, X)
	ELSE IF (basic_results (1,1) > 0.7) THEN
	WRITE (99, 1199) 'location: ' //position, ((PDF_cluster_distr_2 (j,k), j=1,20), k=1,1)
	WRITE (99, 1199) 'droplet_cluster',((PDF_cluster_distr_2 (j,k), j=1,20), k=2,2)
	WRITE (99,'')
	1199 FORMAT (A16, 20F10.5, X)
	ELSE
	WRITE (99, 1190) 'location: ' //position// ': 0.3 < C < 0.7   - no cluster analyses!'
	1190 FORMAT (A55)
	WRITE (99,'')
	END IF

	! write to interparticle arrival file
	IF ((basic_results (1,1) < 0.3) .OR. (basic_results (1,1) > 0.7 )) THEN
	WRITE (854, 857) position, ((PDF_particle_group (j,k), j=1,101), k=1,1)
	WRITE (854, 857) '0-1mm',((PDF_particle_group (j,k), j=1,101), k=2,2)
	WRITE (854, 857) 'poission', ((PDF_particle_group (j,k), j=1,101), k=3,3)
	WRITE (854, 857) '1-3mm',((PDF_particle_group (j,k), j=1,101), k=4,4)
	WRITE (854, 857) 'poission', ((PDF_particle_group (j,k), j=1,101), k=5,5)
	WRITE (854, 857) '3-6mm',((PDF_particle_group (j,k), j=1,101), k=6,6)
	WRITE (854, 857) 'poission', ((PDF_particle_group (j,k), j=1,101), k=7,7)
	WRITE (854, 857) '6-10mm',((PDF_particle_group (j,k), j=1,101), k=8,8)
	WRITE (854, 857) 'poission', ((PDF_particle_group (j,k), j=1,101), k=9,9)
	WRITE (854, 857) '>10mm',((PDF_particle_group (j,k), j=1,101), k=10,10)
	WRITE (854, 857) 'poission', ((PDF_particle_group (j,k), j=1,101), k=11,11)
	WRITE (854,'')
	857 FORMAT (A9, 1000F15.6, X)
	ELSE
	WRITE (854, 1869) 'location: ' //position// ': 0.3 < C < 0.7   - no interparticle arrival analyses!'
	1869 FORMAT (A69)
	WRITE (854,'')

	END IF


	ELSE IF(cluster_or_not == 2) THEN
		!Write to file cluster_summary
	IF ((basic_results (1,1) < 0.3) .OR. (basic_results (1,1) > 0.7 )) THEN
	WRITE (558, 5154) position, ((cluster_results (j,k), j=1,9), k=1,1)
	5154 FORMAT (A4, T10, F10.5, T21, F12.5, T37, F10.5, T53, F10.5, T70, F10.5, T86, F10.5, T107, F10.5, T126, F10.5, T160, F10.5)
	ELSE
	WRITE (558, 5155) position // ': 0.3 < C < 0.7   - no cluster analyses!'
	5155 FORMAT (A44)
	END IF

	!Write to file cluster number PDF
	if (basic_results (1,1) < 0.3) THEN
	WRITE (559, 5158) 'location: ' // position, ((PDF_cluster_distr (j,k), j=1,20), k=1,1)
	WRITE (559, 5158) 'bubble_cluster',((PDF_cluster_distr (j,k), j=1,20), k=2,2)
	WRITE (559,'')
	5158 FORMAT (A15, 20F10.5, X)
	ELSE IF (basic_results (1,1) > 0.7) THEN
	WRITE (559, 5159) 'location: ' //position, ((PDF_cluster_distr (j,k), j=1,20), k=1,1)
	WRITE (559, 5159) 'droplet_cluster',((PDF_cluster_distr (j,k), j=1,20), k=2,2)
	WRITE (559,'')
	5159 FORMAT (A16, 20F10.5, X)
	ELSE
	WRITE (559, 5169) 'location: ' //position// ': 0.3 < C < 0.7   - no cluster analyses!'
	5169 FORMAT (A55)
	WRITE (559,'')
	END IF

		!Write to file cluster_summary constant criterion 1mm
	IF ((basic_results (1,1) < 0.3) .OR. (basic_results (1,1) > 0.7 )) THEN
	WRITE (578, 5174) position, ((cluster_results_1 (j,k), j=1,9), k=1,1)
	5174 FORMAT (A4, T10, F10.5, T21, F12.5, T37, F10.5, T53, F10.5, T70, F10.5, T86, F10.5, T107, F10.5, T126, F10.5, T160, F10.5)
	ELSE
	WRITE (578, 5175) position // ': 0.3 < C < 0.7   - no cluster analyses!'
	5175 FORMAT (A44)
	END IF

	!Write to file cluster number PDF constant criterion 1mm
	if (basic_results (1,1) < 0.3) THEN
	WRITE (579, 5178) 'location: ' // position, ((PDF_cluster_distr_1 (j,k), j=1,20), k=1,1)
	WRITE (579, 5178) 'bubble_cluster',((PDF_cluster_distr_1 (j,k), j=1,20), k=2,2)
	WRITE (579,'')
	5178 FORMAT (A15, 20F10.5, X)
	ELSE IF (basic_results (1,1) > 0.7) THEN
	WRITE (579, 5179) 'location: ' //position, ((PDF_cluster_distr_1 (j,k), j=1,20), k=1,1)
	WRITE (579, 5179) 'droplet_cluster',((PDF_cluster_distr_1 (j,k), j=1,20), k=2,2)
	WRITE (579,'')
	5179 FORMAT (A16, 20F10.5, X)
	ELSE
	WRITE (579, 5189) 'location: ' //position// ': 0.3 < C < 0.7   - no cluster analyses!'
	5189 FORMAT (A55)
	WRITE (579,'')
	END IF

	!Write to file cluster_summary wake criterion
	IF ((basic_results (1,1) < 0.3) .OR. (basic_results (1,1) > 0.7 )) THEN
	WRITE (598, 5194) position, ((cluster_results_2 (j,k), j=1,9), k=1,1)
	5194 FORMAT (A4, T10, F10.5, T25, F12.5, T37, F10.5, T53, F10.5, T70, F10.5, T86, F10.5, T107, F10.5, T126, F10.5, T160, F10.5)
	ELSE
	WRITE (598, 5195) position // ': 0.3 < C < 0.7   - no cluster analyses!'
	5195 FORMAT (A44)
	END IF

	!Write to file cluster number PDF wake criterion
	if (basic_results (1,1) < 0.3) THEN
	WRITE (599, 5198) 'location: ' // position, ((PDF_cluster_distr_2 (j,k), j=1,20), k=1,1)
	WRITE (599, 5198) 'bubble_cluster',((PDF_cluster_distr_2 (j,k), j=1,20), k=2,2)
	WRITE (599,'')
	5198 FORMAT (A15, 20F10.5, X)
	ELSE IF (basic_results (1,1) > 0.7) THEN
	WRITE (599, 5199) 'location: ' //position, ((PDF_cluster_distr_2 (j,k), j=1,20), k=1,1)
	WRITE (599, 5199) 'droplet_cluster',((PDF_cluster_distr_2 (j,k), j=1,20), k=2,2)
	WRITE (599,'')
	5199 FORMAT (A16, 20F10.5, X)
	ELSE
	WRITE (599, 5190) 'location: ' //position// ': 0.3 < C < 0.7   - no cluster analyses!'
	5190 FORMAT (A55)
	WRITE (599,'')
	END IF

	! write to interparticle arrival file time
	IF ((basic_results (1,1) < 0.3) .OR. (basic_results (1,1) > 0.7 )) THEN
	WRITE (954, 957) position, ((PDF_particle_group_time (j,k), j=1,101), k=1,1)
	WRITE (954, 957) '0-0.5mm',((PDF_particle_group_time (j,k), j=1,101), k=2,2)
	WRITE (954, 957) 'poission', ((PDF_particle_group_time (j,k), j=1,101), k=3,3)
	WRITE (954, 957) '0.5-1.5mm',((PDF_particle_group_time (j,k), j=1,101), k=4,4)
	WRITE (954, 957) 'poission', ((PDF_particle_group_time (j,k), j=1,101), k=5,5)
	WRITE (954, 957) '1.5-3mm',((PDF_particle_group_time (j,k), j=1,101), k=6,6)
	WRITE (954, 957) 'poission', ((PDF_particle_group_time (j,k), j=1,101), k=7,7)
	WRITE (954, 957) '3-5mm',((PDF_particle_group_time (j,k), j=1,101), k=8,8)
	WRITE (954, 957) 'poission', ((PDF_particle_group_time (j,k), j=1,101), k=9,9)
	WRITE (954, 957) '>5mm',((PDF_particle_group_time (j,k), j=1,101), k=10,10)
	WRITE (954, 957) 'poission', ((PDF_particle_group_time (j,k), j=1,101), k=11,11)
	WRITE (954,'')
	957 FORMAT (A10, 1000F15.6, X)
	ELSE
	WRITE (954, 1969) 'location: ' //position// ': 0.3 < C < 0.7   - no interparticle arrival analyses!'
	1969 FORMAT (A69)
	WRITE (954,'')

	END IF

END IF


		!Deallocate the arrays 
		DEALLOCATE (Data_array, STAT = status)
		DEALLOCATE (PDF_V, STAT = status)
		DEALLOCATE (PDF_chord_time, STAT = status)
		DEALLOCATE (PDF_chord_size, STAT = status)

		WRITE (*,6000) position
		6000 format (' ', 'Finished calculation of position: ' A)

END DO				! END of main loop in program

	DEALLOCATE (Autocorrelation, STAT =status)
	DEALLOCATE (Crosscorrelation, STAT =status)		
	DEALLOCATE (threshold,STAT=status)
	DEALLOCATE (interfaces, STAT =status)
	DEALLOCATE (basic_results, STAT =status)

	CLOSE (UNIT=50)
	CLOSE (UNIT=52)
	CLOSE (UNIT=54)
	CLOSE (UNIT=56)
	CLOSE (UNIT=60)
	CLOSE (UNIT=62)
	CLOSE (UNIT=64)
	CLOSE (UNIT=66)

	IF (cluster_or_not == 1) THEN
	CLOSE (UNIT=58)
	CLOSE (UNIT=59)
	CLOSE (UNIT=78)
	CLOSE (UNIT=79)
	CLOSE (UNIT=98)
	CLOSE (UNIT=99)
	CLOSE (UNIT=854)

	ELSE IF (cluster_or_not == 2) THEN
	CLOSE (UNIT=558)
	CLOSE (UNIT=559)
	CLOSE (UNIT=578)
	CLOSE (UNIT=579)
	CLOSE (UNIT=598)
	CLOSE (UNIT=599)
	CLOSE (UNIT=954)
	END IF

	pause

end program dataanalysis

!______________________________________________________________________________________________________

SUBROUTINE Correlation (Data_array, number_columns, number_rows, Correlation_steps, &  
			                Autocorrelation, Crosscorrelation, frequency, segment_numbers)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: number_rows
	INTEGER, INTENT (IN) :: Correlation_steps
	INTEGER, INTENT (IN) :: frequency
	INTEGER, INTENT (IN) :: segment_numbers
	REAL(kind=DBL), INTENT (IN) :: Data_array (number_columns, number_rows)
	REAL(kind=DBL), INTENT (OUT) :: Autocorrelation (2, Correlation_steps)
	REAL(kind=DBL), INTENT (OUT) :: Crosscorrelation(2, Correlation_steps)
	INTEGER :: segment_rows							!number of rows in each segment
	INTEGER :: helper								!helper for faster correlations
	INTEGER :: boundary								!helper for faster correlation
	INTEGER :: i =0, j = 0, k=0							!integer for loop
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Autocor_helper	
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: Crosscor_helper	
	REAL(kind=DBL) :: sumX, sumY, sumX2, sumY2, sumXY !Variables for correl
	INTEGER :: offset												 !integer for offsetting
	INTEGER :: status
	REAL(kind=DBL) :: average_value				!Helper to calculate the average values 												
	
	Autocorrelation = 0
	Crosscorrelation = 0
	segment_rows	=0
	helper = 0
	boundary = 0
	sumX = 0
	sumY = 0
	sumX2 = 0
	sumY2 = 0
	sumXY = 0
	offset = 0
	status = 0
	average_value =0


	segment_rows = number_rows/segment_numbers
	boundary = segment_rows - Correlation_steps	


	!Performing Auto-correlation
	ALLOCATE (Autocor_helper(segment_numbers, Correlation_steps))
	Autocor_helper = 0			
	outer : DO i=1, segment_numbers
			helper = segment_rows * (i-1)
		middle : DO j=1, Correlation_steps
			inner : DO k = 1, boundary
					sumX = sumX + Data_array(2 , (k + helper))
					sumX2 = sumX2 + Data_array(2 , (k+ helper))* &
							Data_array(2 , (k+ helper)) 
					sumY = sumY + Data_array(2, (k+ helper+offset))
					sumY2 = sumY2 + Data_array(2, (k+ helper+offset)) &
							*Data_array(2, (k+ helper+offset))
					sumXY = sumXY + Data_array(2,(k+ helper)) * &
							Data_array(2,(k+ helper+offset))
			END DO inner
		
			Autocor_helper(i,j) = (sumXY - sumX * sumY/boundary) /(SQRT(sumX2-sumX*sumX/ &
								  boundary) * SQRT(sumY2-sumY*sumY/boundary))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO middle
		offset = 0		!zero offset helper
	END DO outer
	helper = 0

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps
			Autocorrelation (1, i) = real((i-1))/frequency
			
			DO j = 1, segment_numbers
					average_value = average_value + Autocor_helper (j , i)
			END DO 
				
				Autocorrelation (2, i) = average_value/segment_numbers
				average_value = 0
				
	END DO

	

	WRITE (*,*) 'Finished calculation of autocorrelation! '

!		OPEN (UNIT = 2002, FILE = 'auto.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2002, 2003) ((Autocor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2003 FORMAT (15F10.5)
			
	DEALLOCATE(Autocor_helper)
	!Performing Crosscorrelation
	sumX = 0
	sumX2 = 0
	sumY = 0 
	sumY2 = 0 
	sumXY = 0
	offset = 0
	helper = 0


	ALLOCATE (Crosscor_helper(segment_numbers, Correlation_steps))
	Crosscor_helper = 0

	outer_cross : DO i=1, segment_numbers
		helper = segment_rows * (i-1)
		middle_cross : DO j=1, Correlation_steps
			inner_cross : DO k = 1, boundary
				sumX = sumX + Data_array(2 , (k + 199 + helper))
				sumX2 = sumX2 + Data_array(2 , (k +199 + helper))* &
						Data_array(2 , (k +199 + helper)) 
				sumY = sumY + Data_array(3, (k+ helper+offset))
				sumY2 = sumY2 + Data_array(3, (k+ helper+offset))* &
						Data_array(3, (k+ helper+offset))
				sumXY = sumXY + Data_array(2, (k + 199 + helper)) * &
						Data_array(3, (k+ helper+offset))

			END DO inner_cross
		
			Crosscor_helper(i,j) = (boundary*sumXY - sumX * sumY)/(SQRT(boundary*sumX2- &
									sumX*sumX)* SQRT(boundary*sumY2-sumY*sumY))
			! Adjust variables
			offset=offset +1
			sumX = 0
			sumX2 = 0
			sumY = 0 
			sumY2 = 0 
			sumXY = 0 			

		END DO middle_cross

		offset = 0		!zero offset helper
	END DO outer_cross

	!average the autocorrelation-segments and add time
	DO i = 1, Correlation_steps
			Crosscorrelation (1, i) = -REAL(199)/frequency + real((i-1))/frequency
			
			DO j = 1, segment_numbers
				average_value = average_value + Crosscor_helper (j , i)
			END DO
					
			Crosscorrelation (2, i) = average_value/segment_numbers
			average_value = 0
	END DO
	 
		WRITE (*,*) 'Finished calculation of crosscorrelation! '

!		OPEN (UNIT = 2000, FILE = 'cross.txt', STATUS = 'REPLACE', ACTION ='WRITE', IOSTAT=status)
!		WRITE (2000, 2001) ((Crosscor_helper (i,j), i=1,segment_numbers), j=1, Correlation_steps)
!		2001 FORMAT (15F10.5)
			
	DEALLOCATE(Crosscor_helper)
	
		
END SUBROUTINE Correlation


!_______________________________________________________________________________________________________________


SUBROUTINE Correlation_scales (Autocorrelation, Crosscorrelation, Correlation_steps, &
		   CorScale_results, frequency, delta_x, velocity)
IMPLICIT NONE

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: Correlation_steps
	INTEGER, INTENT (IN) :: frequency
	REAL, INTENT (IN) :: delta_x
	REAL(kind=DBL), INTENT (INOUT) :: velocity
	REAL(kind=DBL), INTENT (IN) :: Autocorrelation (2, Correlation_steps)
	REAL(kind=DBL), INTENT (IN) :: Crosscorrelation (2, Correlation_steps)
	REAL(kind=DBL), INTENT (OUT) ::CorScale_results(14)
	INTEGER :: i = 0 , j = 0, counter, locator
	REAL(kind=DBL) :: helper
	INTEGER :: helper_int
	REAL(kind=DBL) :: min_auto, min_auto_loc, min_cross, min_cross_loc 
	REAL(kind=DBL) :: T_point_five
	REAL(kind=DBL) :: Zero_crossing_auto
	REAL(kind=DBL) :: Txx
	REAL(kind=DBL) :: Rxz_max
	REAL(kind=DBL) :: T_Rxz_max
	REAL(kind=DBL) :: T_half_Rxz_max 
	REAL(kind=DBL) :: Zero_crossing_cross 
	REAL(kind=DBL) :: Txz
	REAL(kind=DBL) :: Tu 

	velocity = 0
	CorScale_results = 0
	helper = 0
	T_point_five = 0
	Zero_crossing_auto =0
	min_auto_loc = 0
	min_auto = 0
	Txx = 0
	Rxz_max = 0
	T_Rxz_max = 0
	Zero_crossing_cross = 0
	min_cross_loc = 0
	min_cross = 0
	T_half_Rxz_max = 0
	Txz = 0
	Tu = 0
	Velocity = 0
	counter=0
	locator = 0
	helper_int = 0


	!Autocorrelation_results
	!Calculation of Time when Autocorrelation is 0.5
	DO i=1, Correlation_steps
		helper = Autocorrelation(2, i)
		IF (helper == 0.5) THEN
			T_point_five = Autocorrelation(1,i)
		EXIT
			ELSE IF (helper < 0.5) THEN
			    T_point_five = Autocorrelation(1,i-1) + (Autocorrelation(1,i) - &
				Autocorrelation(1,i-1)) / (Autocorrelation(2,i) - Autocorrelation(2,i-1))&
				 * (0.5-Autocorrelation(2,i-1))	
			EXIT
		END IF
	END DO

	helper=0
	counter = 0

	!Calculation of Zero crossing of Autocorrelation function
	DO i=1, Correlation_steps
		helper = Autocorrelation(2,i)
		Zero_crossing_auto = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_auto = Autocorrelation(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_auto = Autocorrelation(1,i-1) + (Autocorrelation(1,i) - &
				Autocorrelation(1,i-1)) / (Autocorrelation(2,i) - Autocorrelation(2,i-1)) &
				* (0-Autocorrelation(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_auto = 99	
		END IF
	END DO	

	!Calculation of Txx
	min_auto=1
	locator = 0
	min_auto_loc = 0
	Txx =0
	IF (INT(Zero_crossing_auto) == 99) THEN
		DO j = 1, Correlation_steps-1
			min_auto = Min(Autocorrelation(2,j), min_auto)
		END DO

		DO i = 1, Correlation_steps-1
		locator = locator +1
			IF (Autocorrelation(2,i) - min_auto < 0.00001) THEN
			min_auto_loc = Autocorrelation(1,locator)
			EXIT
			END IF
		END DO

		DO i = 1, locator
			IF (i ==1000) THEN
			Txx = Txx + (Autocorrelation(2,i)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation(2,i)+Autocorrelation(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = 1, counter
			IF (j ==1000) THEN
			Txx = Txx + (Autocorrelation(2,j)/frequency)
			ELSE
			Txx = Txx + (Autocorrelation(2,j)+Autocorrelation(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_auto=0
		min_auto_loc=Zero_crossing_auto
	END IF
		
	helper =0
	counter =0
	Rxz_max = 0
	!Crosscorrelation_results
	!Find maximum Crosscorrelation value (Rxz)max
	DO i = 1, Correlation_steps
		helper = Crosscorrelation(2,i)
		Rxz_max = Max(Rxz_max,helper)
	END DO	
		
	helper = 0 
	locator = 0
	T_Rxz_max = 0

	!Corresponding Time for Rxz_max
	DO i = 1, Correlation_steps
		locator = locator +1
		IF (Rxz_max - Crosscorrelation(2,i) < 0.00001) THEN
		T_Rxz_max = Crosscorrelation(1,locator)
		EXIT 
		END IF
	END DO

	T_half_Rxz_max = 0

	!Calculation of Time for 0.5*Rxz_max
	DO i=locator, Correlation_steps
		IF (Crosscorrelation (2,i) <= 0.5*Rxz_max) THEN
			T_half_Rxz_max = Crosscorrelation(1,i-1) + (Crosscorrelation(1,i) - &
							Crosscorrelation(1,i-1)) / (Crosscorrelation(2,i) - &
							Crosscorrelation(2,i-1)) * (0.5*Rxz_max-Crosscorrelation(2,i-1))
			EXIT
		END IF
	END DO

	counter =-1
	helper = 0

	!Time where the values of Rxz_max cross the x-axis for the first time right of Rxy_max
	DO i=locator, Correlation_steps
		helper = Crosscorrelation(2,i)
		Zero_crossing_cross = 0
		counter = counter +1
		IF (helper == 0) THEN
			Zero_crossing_cross = Crosscorrelation(1,i)
				EXIT
			ELSE IF (helper < 0.) THEN
				Zero_crossing_cross = Crosscorrelation(1,i-1) + (Crosscorrelation(1,i) - &
			    Crosscorrelation(1,i-1))  / (Crosscorrelation(2,i) - Crosscorrelation(2,i-1))&
				* (0-Crosscorrelation(2,i-1))
				EXIT
				ELSE 
				Zero_crossing_cross = 98	
		END IF
	END DO	

	!Calculation of Txz
	min_cross = Rxz_max
	min_cross_loc = 0
	Txz = 0
	IF (INT(Zero_crossing_cross) == 98) THEN
		DO j = locator, Correlation_steps-1
			min_cross = Min(Crosscorrelation(2,j), min_cross)
		END DO
		
		helper_int = 0
		DO i = locator, Correlation_steps-1
		helper_int = helper_int +1
			IF (Crosscorrelation(2,i) - min_cross < 0.00001) THEN
			min_cross_loc = Crosscorrelation(1,(helper_int+locator-1))
			EXIT
			END IF
		END DO

		DO i = locator, (helper_int +locator-1)
			IF (i ==1000) THEN
			Txz = Txz + (Crosscorrelation(2,i)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation(2,i)+Crosscorrelation(2,i+1))/(2*frequency)
			END IF
		END DO
	ELSE 
		DO j = locator, (locator+counter)
			IF (j ==1000) THEN
			Txz = Txz + (Crosscorrelation(2,j)/frequency)
			ELSE
			Txz = Txz + (Crosscorrelation(2,j)+Crosscorrelation(2,j+1))/(2*frequency)
			END IF
		END DO 
		min_cross = 0
		min_cross_loc = Zero_crossing_cross
	END IF
		
	!Calculation of turbulence intensity
	IF (((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five) <0) THEN
		 Tu = 97
	ELSE
	Tu = 0.851 * SQRT((T_half_Rxz_max - T_Rxz_max)*(T_half_Rxz_max - T_Rxz_max)- &
		 T_point_five*T_point_five)/T_Rxz_max
	END IF
		
	!Calculation of time averaged interfacial velocity
	velocity = delta_x/1000/T_Rxz_max	

	!store results into array CorScale_results  :T_point_five, Zero_crossing_auto, min_auto_loc, min_auto,Txx, &
						!T_Rxz_max, T_half_Rxz_max, Zero_crossing_cross, min_cross_loc, min_cross,Txz, Tu, Velocity 
	CorScale_results = (/ T_point_five, Zero_crossing_auto, min_auto_loc, min_auto, Txx, Rxz_max, T_Rxz_max, &
						Zero_crossing_cross, min_cross_loc, min_cross, T_half_Rxz_max, Txz, Tu, Velocity /)


END SUBROUTINE Correlation_scales

!_______________________________________________________________________________________________________________

SUBROUTINE Probability_V (Data_array, number_columns, number_rows, PDF_V_segment_size, & 
						  PDF_segments, PDF_V, threshold, threshold_C) 

IMPLICIT NONE

INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: number_rows
	REAL, INTENT (IN) :: PDF_V_segment_size
	REAL(kind=DBL), INTENT (IN) :: PDF_segments
	REAL, INTENT (IN) :: threshold_C
	REAL(kind=DBL), INTENT (IN) :: Data_array (number_columns, number_rows)
	REAL(kind=DBL), INTENT (INOUT) :: PDF_V (number_columns, Int(PDF_segments))
	REAL(kind=DBL), INTENT (OUT) :: threshold(number_columns-1)
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_helper
	INTEGER :: i, j, k, l, locator
	REAL(kind=DBL) :: segmenter
	REAL(kind=DBL) :: Max_air, Max_water
	REAL(kind=DBL) :: max_air_loc, max_water_loc
	INTEGER :: PDF_upper_boundary_air

	threshold = 0
	locator = 0
	segmenter = 0
	Max_air = 0
	Max_water = 0
	max_air_loc = 0
	max_water_loc = 0
	PDF_upper_boundary_air = 0

	PDF_upper_boundary_air = 2.5 / PDF_V_segment_size

	ALLOCATE (PDF_helper(number_columns-1, Int(PDF_segments)))
	PDF_helper = 0

	DO i = 1, (number_columns-1)				!number of devices
		DO j = 1, number_rows					!number of rows in data file
			DO k = 1, Int(PDF_segments)			!
				segmenter = -1.0 + k * PDF_V_segment_size

				IF (Data_array(i+1, j) < segmenter) THEN 
					PDF_helper(i,k) = PDF_helper(i,k) +1
					EXIT
				END IF
			END DO
		END DO
	END DO

	DO i = 1, Int(PDF_segments)
		PDF_V(1,i) = -1.0 + i * PDF_V_segment_size
	END DO

	DO i =1, (number_columns -1)
		DO j = 1, Int(PDF_segments)
			PDF_V(i+1,j) = REAL(PDF_helper(i,j))/REAL(number_rows)
		END DO
	END DO

	DEALLOCATE (PDF_helper)

	!find two max-values of the bimodal PDF function
	DO i = 1, number_columns-1
		DO j = 1, Int(PDF_upper_boundary_air)
		Max_air = max(Max_air, PDF_V(i+1,j))
		END DO

		DO k = Int(PDF_segments), Int(PDF_upper_boundary_air), -1
		Max_water =  max(Max_water, PDF_V(i+1,k))
		END DO

		DO l = 1, Int(PDF_upper_boundary_air)
			locator = locator +1
			IF (Max_air - PDF_V(i+1,l) < 0.00001) THEN
			max_air_loc = PDF_V(1,locator)
			EXIT
			END IF
		END DO

		locator = Int(PDF_segments)

		DO l = Int(PDF_segments), Int(PDF_upper_boundary_air), -1
			IF (Max_water - PDF_V(i+1,l) < 0.00001) THEN
			max_water_loc = PDF_V(1,locator)
			EXIT
			END IF
			locator = locator - 1
		END DO


	threshold(i) = max_air_loc + (max_water_loc - max_air_loc)*threshold_C
	Max_water = 0
	Max_air = 0
	locator = 0
	END DO


END SUBROUTINE Probability_V

!_______________________________________________________________________________________________________________

SUBROUTINE Basic_properties (Data_array, number_columns, number_rows, threshold,&
							 interfaces, basic_results, max_no_bubbles, frequency)

INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: number_rows
	INTEGER, INTENT (IN) :: max_no_bubbles
	INTEGER, INTENT(IN) :: frequency
	REAL(kind=DBL), INTENT (IN) :: Data_array (number_columns, number_rows)
	REAL(kind=DBL), INTENT (IN) :: threshold (number_columns-1)
	INTEGER, INTENT (INOUT) :: interfaces((number_columns-1)*2, max_no_bubbles)  
	REAL(kind=DBL), INTENT (INOUT) :: basic_results ((number_columns-1)*2, 1)
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: instantenious_C !helper -> air=1, water =0 
	INTEGER i, j, status	
	INTEGER :: size_counter, bubble_counter, C_counter

	ALLOCATE (instantenious_C(number_columns-1, number_rows))
	instantenious_C = 0
	interfaces = 0
	basic_results = 0
	status = 0
	size_counter = 0
	bubble_counter = 0
	C_counter = 0

	DO i =1, (number_columns-1)
		DO j = 1, number_rows
			IF (Data_array((i+1),j) - threshold (i) > 0.00001) THEN
				instantenious_C(i, j) = 1
				ELSE
				instantenious_C(i, j) = 0
			END IF
		END DO
	END DO	
		

	DO i =1, (number_columns-1)					!for all probes
	C_counter =0
		IF (instantenious_C(i,1) == 0) THEN
			size_counter = 0
			bubble_counter = 0
			ELSE 
			size_counter = 1
			bubble_counter = 1
		END IF

		DO j = 1, number_rows-1
			IF (instantenious_C(i,j) == 1) THEN
			C_counter = C_counter +1
			END IF

			IF (instantenious_C(i,j) == 0 .AND. instantenious_C(i,j+1) == 1 ) THEN
				size_counter = size_counter +1
				interfaces((i*2-1), size_counter) = j+1				!change water to air
				bubble_counter = bubble_counter +1
			ELSE IF (instantenious_C(i,j) == 1 .AND. instantenious_C(i,j+1) == 0 ) THEN	
				interfaces((i*2), size_counter) = j	+1			!change air to water	
			END IF
		END DO

		IF (size_counter == 0) THEN 
		interfaces ((i*2), 1) = number_rows
		ELSE
		interfaces((i*2), size_counter) = number_rows
		END IF

		basic_results(i*2-1,1) = REAL(C_counter)/REAL(number_rows)
		basic_results(i*2, 1) = REAL(bubble_counter)*REAL(frequency)/REAL(number_rows)
	
		
	END DO
	
	DEALLOCATE (instantenious_C)

END SUBROUTINE Basic_properties

!_______________________________________________________________________________________________

SUBROUTINE chord_length (Interfaces, number_columns, max_no_bubbles, velocity, basic_results, &
					   PDF_chord_size, PDF_number_chord, sample_duration, &
					   frequency, PDF_Court_max_length, PDF_Court_min_length, PDF_Court_length_segment_size)

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT (IN) :: velocity
	REAL(kind=DBL), INTENT(IN) :: PDF_number_chord
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency
	REAL, INTENT (IN) :: PDF_Court_max_length 
	REAL, INTENT (IN) :: PDF_Court_min_length 
	REAL, INTENT (IN) :: PDF_Court_length_segment_size 								!for the chord length distributions
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	REAL(kind=DBL), INTENT (INOUT) :: PDF_chord_size(INT(PDF_number_chord),3)
	INTEGER :: i, j, k, l , status
	REAL :: bubble_number
	REAL :: segmenter_size
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_length
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_length
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_helper_length

		
		bubble_number = basic_results(2,1)* REAL(sample_duration)

		ALLOCATE (bubble_length(1, INT(bubble_number)))	
		ALLOCATE (PDF_helper_length(2, INT(PDF_number_chord)))
		bubble_length = 0.0
		PDF_helper_length = 0
		status = 0
		segmenter_size = 0.0

		DO j = 1, INT(bubble_number)
				bubble_length(1,j) = REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000*velocity
				DO k = 1, Int(PDF_number_chord)			
					segmenter_size = PDF_Court_min_length + k * PDF_Court_length_segment_size
						IF (bubble_length(1, j) - segmenter_size < 0.00001) THEN
						PDF_helper_length(1,k) = PDF_helper_length(1,k) +1
						EXIT
						ELSE IF (bubble_length(1, j) > REAL(PDF_Court_max_length)) THEN
						PDF_helper_length(1,PDF_number_chord) = PDF_helper_length(1,PDF_number_chord) +1
						EXIT
						END IF			
				END DO
		END DO

		DEALLOCATE(bubble_length, STAT = status)
		ALLOCATE (droplet_length(1, (INT(bubble_number)-1)))
		droplet_length = 0.0
		segmenter_size =0.0

		DO j = 1, (INT(bubble_number)-1)
				droplet_length(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
				DO k = 1, Int(PDF_number_chord)			
					segmenter_size = PDF_Court_min_length + k * PDF_Court_length_segment_size
						IF (droplet_length(1, j) - segmenter_size < 0.00001) THEN
						PDF_helper_length(2,k) = PDF_helper_length(2,k) +1
						EXIT
						ELSE IF (droplet_length(1, j) > REAL(PDF_Court_max_length)) THEN
						PDF_helper_length(2,PDF_number_chord) = PDF_helper_length(2,PDF_number_chord) +1
						EXIT
						END IF
				END DO
		END DO	

		DEALLOCATE(droplet_length, STAT = status)
					
		DO j = 1, Int(PDF_number_chord)
			PDF_chord_size(j,1) = PDF_Court_min_length + (j-1) * PDF_Court_length_segment_size
			PDF_chord_size(j,2) = REAL(PDF_helper_length(1,j))/bubble_number
			PDF_chord_size(j,3) = REAL(PDF_helper_length(2,j))/(bubble_number-1)

		END DO
			
		DEALLOCATE (PDF_helper_length, STAT = status)
		
END SUBROUTINE chord_length 

!____________________________________________________________________________________________________

SUBROUTINE chord_time (Interfaces, number_columns, max_no_bubbles, basic_results, &
					   PDF_chord_time, PDF_number_time, sample_duration, &
					   frequency,  PDF_Court_max_time, PDF_Court_min_time, PDF_Court_time_segment_size)

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT(IN) :: PDF_number_time
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency								
	REAL, INTENT (IN) :: PDF_Court_max_time 			
	REAL, INTENT (IN) :: PDF_Court_min_time
	REAL, INTENT (IN) :: PDF_Court_time_segment_size 								!for the chord time distributions
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	REAL(kind=DBL), INTENT (INOUT) :: PDF_chord_time(PDF_number_time,3)
	INTEGER :: i, j, k, l , status
	REAL :: bubble_number
	REAL :: segmenter_time
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_time 
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_time
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_helper_time
	
	bubble_number = basic_results(2,1)* REAL(sample_duration)

		ALLOCATE (bubble_time(1, INT(bubble_number)))
		ALLOCATE (PDF_helper_time(2, INT(PDF_number_time)))
		bubble_time = 0.0	
		PDF_helper_time = 0
		status = 0
		segmenter_time = 0.0

		DO j = 1, INT(bubble_number)
				bubble_time(1,j) = REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
				DO k = 1, Int(PDF_number_time)	
					segmenter_time = PDF_Court_min_time + k * PDF_Court_time_segment_size
						IF (bubble_time(1, j) - segmenter_time < 0.00001) THEN
						PDF_helper_time(1,k) = PDF_helper_time(1,k) +1
						EXIT
						ELSE IF (bubble_time(1, j) > REAL(PDF_Court_max_time)) THEN
						PDF_helper_time(1,PDF_number_time) = PDF_helper_time(1,PDF_number_time) +1
						EXIT
						END IF
				END DO				
		END DO

		DEALLOCATE(bubble_time, STAT = status)
		ALLOCATE (droplet_time(1, (INT(bubble_number)-1)))
		droplet_time = 0.0
		segmenter_time =0.0

		DO j = 1, (INT(bubble_number)-1)
				droplet_time(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000					
				DO k = 1, Int(PDF_number_time)			
					segmenter_time = PDF_Court_min_time + k * PDF_Court_time_segment_size
						IF (droplet_time(1, j) - segmenter_time < 0.00001) THEN
						PDF_helper_time(2,k) = PDF_helper_time(2,k) +1
						EXIT
						ELSE IF (droplet_time(1, j) > REAL(PDF_Court_max_time)) THEN
						PDF_helper_time(2,PDF_number_time) = PDF_helper_time(2,PDF_number_time) +1
						EXIT
						END IF
				END DO					
		END DO

		DEALLOCATE(droplet_time, STAT = status)

		DO j = 1, Int(PDF_number_time)
			PDF_chord_time(j,1) = PDF_Court_min_time + (j-1) * PDF_Court_time_segment_size
			PDF_chord_time(j,2) = REAL(PDF_helper_time(1,j))/bubble_number
			PDF_chord_time(j,3) = REAL(PDF_helper_time(2,j))/(bubble_number-1)
		END DO
DEALLOCATE (PDF_helper_time, STAT = status)

END SUBROUTINE chord_time


!____________________________________________________________________________________________________

SUBROUTINE 	Cluster_analysis (interfaces, number_columns, max_no_bubbles, basic_results, velocity, sample_duration, &
					   frequency, cluster_factor, PDF_cluster_distr, cluster_results)

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	REAL(kind=DBL), INTENT (IN) :: velocity
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency
	REAL, INTENT(IN) :: cluster_factor
	REAL(kind=DBL), INTENT (INOUT) :: PDF_cluster_distr (20,2)
	REAL(kind=DBL), INTENT (INOUT) :: cluster_results (9, 1)

		
	
!	REAL(kind=DBL), INTENT (INOUT) :: PDF_chord_size(INT(PDF_number_chord),3)
	INTEGER :: i, j, k, l , status
	INTEGER :: bubble_number
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord_swapper
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord_swapper

	INTEGER :: iptr, sum_up, counter, detector
	REAL(kind=DBL) :: temp

	REAL(kind=DBL) :: median_bubble_chord
	REAL(kind=DBL) :: median_droplet_chord
	REAL(kind=DBL) :: median_bubble_chord_crit
	REAL(kind=DBL) :: median_droplet_chord_crit
	REAL(kind=DBL) :: bubble_chord_sum
	REAL(kind=DBL) :: droplet_chord_sum
	REAL(kind=DBL) :: average_bubble_chord
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: instantenious_cluster
	REAL(kind=DBL) :: size
	REAL(kind=DBL) :: chord_size_cluster
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_cluster_no
	INTEGER :: cluster_container

	REAL(kind=DBL) ::  bubbles_in_cluster
	REAL(kind=DBL) ::  clusters_per_second
	REAL(kind=DBL) :: average_bubbles_cluster
	REAL(kind=DBL) :: average_chord_size_cluster
	REAL(kind=DBL) :: cluster_ratio
	REAL(kind=DBL) :: chord_lead_particle
	REAL(kind=DBL) :: chord_cluster
	REAL(kind=DBL) :: chord_average_cluster
	REAL(kind=DBL) :: chord_average_cluster_sum

	iptr = 0
	sum_up = 0
	counter = 0
	detector = 0
	temp = 0.0
	median_bubble_chord = 0.0
    median_droplet_chord = 0.0
	median_bubble_chord_crit = 0.0
	median_droplet_chord_crit = 0.0
	average_bubble_chord = 0.0
	size = 0.0
	chord_size_cluster = 0.0
	cluster_container = 0
	bubbles_in_cluster = 0.0
	clusters_per_second = 0.0
	average_bubbles_cluster = 0.0
	average_chord_size_cluster = 0.0
	cluster_ratio = 0.0
	chord_lead_particle = 0.0
	chord_cluster = 0.0
	chord_average_cluster = 0.0
	chord_average_cluster_sum = 0.0



	bubble_number = INT(basic_results(2,1)* REAL(sample_duration))
	
		ALLOCATE (bubble_chord(1, bubble_number))
		ALLOCATE (bubble_chord_swapper(1, bubble_number))	
		ALLOCATE (droplet_chord(1, (bubble_number-1)))
		ALLOCATE (droplet_chord_swapper(1, bubble_number -1))
		droplet_chord = 0.0
		droplet_chord_swapper = 0.0
		bubble_chord = 0.0
		bubble_chord_swapper = 0.0
		ALLOCATE (PDF_cluster_no (1, 20))
		PDF_cluster_no = 0
		bubble_chord_sum = 0.0
		droplet_chord_sum = 0.0

	! Air bubble clustering for C < 30%
IF (basic_results(1,1) < 0.3) THEN

		DO j = 1, bubble_number
				bubble_chord(1,j) = REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000*velocity
		END DO

		DO j = 1, bubble_number - 1
			droplet_chord(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
			droplet_chord_swapper(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
		END DO


			!calculate median water chord
			DO i = 1, bubble_number-2
				iptr = i
				DO j = i+1, bubble_number - 1
					IF ( droplet_chord_swapper(1,j) < droplet_chord_swapper(1,iptr)) THEN
						iptr = j
					END IF
				END DO
			IF (i /= iptr) THEN 
			temp = droplet_chord_swapper(1,i)
			droplet_chord_swapper(1,i) = droplet_chord_swapper(1,iptr)
			droplet_chord_swapper(1,iptr) = temp
			END IF

			END DO

			IF ( mod((bubble_number - 1), 2) == 0) THEN 
				median_droplet_chord = ( droplet_chord_swapper (1,(bubble_number-1)/2) + droplet_chord_swapper (1,(bubble_number-1)/2+1) ) / 2
			ELSE
				median_droplet_chord = ( droplet_chord_swapper (1,(bubble_number-1)/2 + 1) )
			END IF


		!Calculate median droplet chord criterion value for cluster boundary criterion
			median_droplet_chord_crit = median_droplet_chord * cluster_factor

	!calculate average bubble chord size
	DO j = 1, bubble_number
			bubble_chord_sum = bubble_chord_sum + REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000*velocity
	END DO
	average_bubble_chord = bubble_chord_sum / REAL(bubble_number)

	! calculate air bubble cluster statistics in bubble flow region
		ALLOCATE (instantenious_cluster(1, bubble_number - 1))
	DO i = 1, (bubble_number - 1)
		IF ( droplet_chord(1, i) <= median_droplet_chord_crit) THEN
			instantenious_cluster(1, i) = 1
		ELSE
			instantenious_cluster(1, i) = 0
		END IF
	END DO

	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

	DO i = 1, (bubble_number - 2)
		sum_up = sum_up + instantenious_cluster(1, i)
		IF (instantenious_cluster(1, i) == 1) THEN
		size = size + bubble_chord(1, i)
		END IF

		IF (instantenious_cluster(1, i) ==1 .AND. instantenious_cluster(1, i +1) == 0) THEN
		counter = counter + sum_up														! number of bubbles in clusters

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up <= cluster_container) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO

		detector = detector +1															! number of clusters

		chord_size_cluster = chord_size_cluster + (size + bubble_chord(1, i+1))			! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = bubble_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + bubble_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))
		sum_up = 1
		size = 0

		END IF

	END DO

	DEALLOCATE(instantenious_cluster, STAT = status)

	DO k = 1, 20
		PDF_cluster_distr (k,1) = k + 1
		PDF_cluster_distr (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number) *100		! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results (1, 1) = basic_results(1,1)
	cluster_results (2, 1) = REAL(counter)
	cluster_results (3, 1) = bubbles_in_cluster
	cluster_results (4, 1) = REAL(detector)
	cluster_results (5, 1) = clusters_per_second
	cluster_results (6, 1) = average_bubbles_cluster
	cluster_results (7, 1) = average_chord_size_cluster
	cluster_results (8, 1) = cluster_ratio
	cluster_results (9, 1) = chord_average_cluster


		
! Water droplet clustering for C > 70%	
ELSE IF (basic_results(1,1) > 0.7) THEN

		DO j = 1, (bubble_number-1)
				droplet_chord(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
		END DO

		DO j = 1, bubble_number
			bubble_chord(1,j) = REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000*velocity
			bubble_chord_swapper(1,j) = REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000*velocity
		END DO

	!calculate median air chord
			DO i = 1, bubble_number-1
				iptr = i
				DO j = i+1, bubble_number
					IF ( bubble_chord_swapper(1,j) < bubble_chord_swapper(1,iptr)) THEN
						iptr = j
					END IF
				END DO
			IF (i /= iptr) THEN 
			temp = bubble_chord_swapper(1,i)
			bubble_chord_swapper(1,i) = bubble_chord_swapper(1,iptr)
			bubble_chord_swapper(1,iptr) = temp
			END IF

			END DO

			IF ( mod(bubble_number, 2) == 0) THEN 
				median_bubble_chord = ( bubble_chord_swapper (1,bubble_number/2) + bubble_chord_swapper (1,bubble_number/2+1) ) / 2
			ELSE
				median_bubble_chord = ( bubble_chord_swapper (1,bubble_number/2 + 1) )
			END IF

		!Calculate median droplet chord criterion value for cluster boundary criterion
			median_bubble_chord_crit = median_bubble_chord * cluster_factor

	!calculate average bubble chord size
	DO j = 1, bubble_number - 1
			droplet_chord_sum = droplet_chord_sum + REAL(interfaces(1,j+1) - interfaces((2),j))/REAL(frequency)*1000*velocity
	END DO
	average_bubble_chord = droplet_chord_sum / REAL(bubble_number-1)


	! calculate air bubble cluster statistics in bubble flow region
		ALLOCATE (instantenious_cluster(1, bubble_number))
	DO i = 1, (bubble_number)
		IF ( bubble_chord(1, i) <= median_bubble_chord_crit) THEN
			instantenious_cluster(1, i) = 1
		ELSE
			instantenious_cluster(1, i) = 0
		END IF
	END DO

	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

	DO i = 1, (bubble_number-2)
		sum_up = sum_up + instantenious_cluster(1, i)
		IF (instantenious_cluster(1, i) == 1) THEN
		size = size + droplet_chord(1, i)
		END IF

		IF (instantenious_cluster(1, i) ==1 .AND. instantenious_cluster(1, i +1) == 0) THEN
		counter = counter + sum_up														! number of bubbles in clusters

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up <= cluster_container) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO

		detector = detector +1															! number of clusters

		chord_size_cluster = chord_size_cluster + (size + droplet_chord(1, i +1))			! sum of bubble chord sizes of clusters
		
		chord_cluster = 0.0

		chord_lead_particle = droplet_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + droplet_chord(1, (i +1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))
		
		sum_up = 1
		size = 0

		END IF

	END DO

		DEALLOCATE(instantenious_cluster, STAT = status)

	DO k = 1, 20
		PDF_cluster_distr (k,1) = k + 1
		PDF_cluster_distr (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number-1)	 *100		! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results (1, 1) = basic_results(1,1)
	cluster_results (2, 1) = REAL(counter)
	cluster_results (3, 1) = bubbles_in_cluster
	cluster_results (4, 1) = REAL(detector)
	cluster_results (5, 1) = clusters_per_second
	cluster_results (6, 1) = average_bubbles_cluster
	cluster_results (7, 1) = average_chord_size_cluster
	cluster_results (8, 1) = cluster_ratio
	cluster_results (9, 1) = chord_average_cluster

		DEALLOCATE(PDF_cluster_no, STAT = status)

END IF

		DEALLOCATE(bubble_chord, STAT = status)
		DEALLOCATE(bubble_chord_swapper, STAT = status)
		DEALLOCATE(droplet_chord, STAT = status)
		DEALLOCATE(droplet_chord_swapper, STAT = status)


END SUBROUTINE Cluster_analysis


!____________________________________________________________________________________________________

SUBROUTINE 	Cluster_analysis_1 (interfaces, number_columns, max_no_bubbles, basic_results, velocity, sample_duration, &
					   frequency, cluster_factor_1, PDF_cluster_distr_1, cluster_results_1)

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	REAL(kind=DBL), INTENT (IN) :: velocity
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency
	REAL, INTENT (IN) :: cluster_factor_1
	REAL(kind=DBL), INTENT (INOUT) :: PDF_cluster_distr_1 (20,2)
	REAL(kind=DBL), INTENT (INOUT) :: cluster_results_1 (9, 1)

		
	
!	REAL(kind=DBL), INTENT (INOUT) :: PDF_chord_size(INT(PDF_number_chord),3)
	INTEGER :: i, j, k, l , status
	INTEGER :: bubble_number
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord_swapper
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord_swapper

	INTEGER :: iptr, sum_up, counter, detector
	REAL(kind=DBL) :: temp

	REAL(kind=DBL) :: median_bubble_chord
	REAL(kind=DBL) :: median_droplet_chord
	REAL(kind=DBL) :: median_bubble_chord_crit
	REAL(kind=DBL) :: median_droplet_chord_crit
	REAL(kind=DBL) :: bubble_chord_sum
	REAL(kind=DBL) :: droplet_chord_sum
	REAL(kind=DBL) :: average_bubble_chord
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: instantenious_cluster
	REAL(kind=DBL) :: size
	REAL(kind=DBL) :: chord_size_cluster
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_cluster_no
	INTEGER :: cluster_container

	REAL(kind=DBL) ::  bubbles_in_cluster
	REAL(kind=DBL) ::  clusters_per_second
	REAL(kind=DBL) :: average_bubbles_cluster
	REAL(kind=DBL) :: average_chord_size_cluster
	REAL(kind=DBL) :: cluster_ratio
	REAL(kind=DBL) :: chord_lead_particle
	REAL(kind=DBL) :: chord_cluster
	REAL(kind=DBL) :: chord_average_cluster
	REAL(kind=DBL) :: chord_average_cluster_sum

	iptr = 0
	sum_up = 0
	counter = 0
	detector = 0
	temp = 0.0
	median_bubble_chord = 0.0
    median_droplet_chord = 0.0
	median_bubble_chord_crit = 0.0
	median_droplet_chord_crit = 0.0
	average_bubble_chord = 0.0
	size = 0.0
	chord_size_cluster = 0.0
	cluster_container = 0
	bubbles_in_cluster = 0.0
	clusters_per_second = 0.0
	average_bubbles_cluster = 0.0
	average_chord_size_cluster = 0.0
	cluster_ratio = 0.0
	chord_lead_particle = 0.0
	chord_cluster = 0.0
	chord_average_cluster = 0.0
	chord_average_cluster_sum = 0.0
	bubble_chord_sum = 0.0
	droplet_chord_sum = 0.0

	bubble_number = INT(basic_results(2,1)* REAL(sample_duration))
	
		ALLOCATE (bubble_chord(1, bubble_number))
		ALLOCATE (bubble_chord_swapper(1, bubble_number))	
		ALLOCATE (droplet_chord(1, (bubble_number-1)))
		ALLOCATE (droplet_chord_swapper(1, bubble_number -1))
		droplet_chord = 0.0
		droplet_chord_swapper = 0.0
		bubble_chord = 0.0
		bubble_chord_swapper = 0.0
		ALLOCATE (PDF_cluster_no (1, 20))
		PDF_cluster_no = 0

	! Air bubble clustering for C < 30%
IF (basic_results(1,1) < 0.3) THEN

		DO j = 1, bubble_number
				bubble_chord(1,j) = REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000*velocity
		END DO

		DO j = 1, bubble_number - 1
			droplet_chord(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
			droplet_chord_swapper(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
		END DO

	!calculate average bubble chord size
	DO j = 1, bubble_number
			bubble_chord_sum = bubble_chord_sum + REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000*velocity
	END DO
	average_bubble_chord = bubble_chord_sum / REAL(bubble_number)


	! calculate air bubble cluster statistics in bubble flow region
		ALLOCATE (instantenious_cluster(1, bubble_number - 1))
	DO i = 1, (bubble_number - 1)
		IF ( droplet_chord(1, i) <= (cluster_factor_1 * 1.0)) THEN
			instantenious_cluster(1, i) = 1
		ELSE
			instantenious_cluster(1, i) = 0
		END IF
	END DO

	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

	DO i = 1, (bubble_number - 2)
		sum_up = sum_up + instantenious_cluster(1, i)
		IF (instantenious_cluster(1, i) == 1) THEN
		size = size + bubble_chord(1, i)
		END IF

		IF (instantenious_cluster(1, i) ==1 .AND. instantenious_cluster(1, i +1) == 0) THEN
		counter = counter + sum_up														! number of bubbles in clusters

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up <= cluster_container) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO

		detector = detector +1															! number of clusters

		chord_size_cluster = chord_size_cluster + (size + bubble_chord(1, i+1)) 			! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = bubble_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + bubble_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))

		
		sum_up = 1
		size = 0

		END IF

	END DO

	DEALLOCATE(instantenious_cluster, STAT = status)

	DO k = 1, 20
		PDF_cluster_distr_1 (k,1) = k + 1
		PDF_cluster_distr_1 (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number) *100			! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results_1 (1, 1) = basic_results(1,1)
	cluster_results_1 (2, 1) = REAL(counter)
	cluster_results_1 (3, 1) = bubbles_in_cluster
	cluster_results_1 (4, 1) = REAL(detector)
	cluster_results_1 (5, 1) = clusters_per_second
	cluster_results_1 (6, 1) = average_bubbles_cluster
	cluster_results_1 (7, 1) = average_chord_size_cluster
	cluster_results_1 (8, 1) = cluster_ratio
	cluster_results_1 (9, 1) = chord_average_cluster

		
! Water droplet clustering for C > 70%	
ELSE IF (basic_results(1,1) > 0.7) THEN

		DO j = 1, (bubble_number-1)
				droplet_chord(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
		END DO

		DO j = 1, bubble_number
			bubble_chord(1,j) = REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000*velocity
			bubble_chord_swapper(1,j) = REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000*velocity
		END DO

	!calculate average bubble chord size
	DO j = 1, bubble_number - 1
			droplet_chord_sum = droplet_chord_sum + REAL(interfaces(1,j+1) - interfaces((2),j))/REAL(frequency)*1000*velocity
	END DO
	average_bubble_chord = droplet_chord_sum / REAL(bubble_number-1)


	! calculate air bubble cluster statistics in bubble flow region
		ALLOCATE (instantenious_cluster(1, bubble_number))
	DO i = 1, (bubble_number)
		IF ( bubble_chord(1, i) <= (cluster_factor_1 * 1.0)) THEN
			instantenious_cluster(1, i) = 1
		ELSE
			instantenious_cluster(1, i) = 0
		END IF
	END DO

	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

	DO i = 1, (bubble_number-2)
		sum_up = sum_up + instantenious_cluster(1, i)
		IF (instantenious_cluster(1, i) == 1) THEN
		size = size + droplet_chord(1, i)
		END IF

		IF (instantenious_cluster(1, i) ==1 .AND. instantenious_cluster(1, i +1) == 0) THEN
		counter = counter + sum_up														! number of bubbles in clusters

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up <= cluster_container) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO

		detector = detector +1															! number of clusters

		chord_size_cluster = chord_size_cluster + (size + droplet_chord(1, i+1))			! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = droplet_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + droplet_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))

		
		sum_up = 1
		size = 0

		END IF

	END DO

		DEALLOCATE(instantenious_cluster, STAT = status)

	DO k = 1, 20
		PDF_cluster_distr_1 (k,1) = k + 1
		PDF_cluster_distr_1 (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number-1) *100			! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results_1 (1, 1) = basic_results(1,1)
	cluster_results_1 (2, 1) = REAL(counter)
	cluster_results_1 (3, 1) = bubbles_in_cluster
	cluster_results_1 (4, 1) = REAL(detector)
	cluster_results_1 (5, 1) = clusters_per_second
	cluster_results_1 (6, 1) = average_bubbles_cluster
	cluster_results_1 (7, 1) = average_chord_size_cluster
	cluster_results_1 (8, 1) = cluster_ratio
	cluster_results_1 (9, 1) = chord_average_cluster

		DEALLOCATE(PDF_cluster_no, STAT = status)

END IF

		DEALLOCATE(bubble_chord, STAT = status)
		DEALLOCATE(bubble_chord_swapper, STAT = status)
		DEALLOCATE(droplet_chord, STAT = status)
		DEALLOCATE(droplet_chord_swapper, STAT = status)


END SUBROUTINE Cluster_analysis_1


!____________________________________________________________________________________________________

SUBROUTINE 	Cluster_analysis_wake (interfaces, number_columns, max_no_bubbles, basic_results, velocity, sample_duration, &
					   frequency, PDF_cluster_distr_2, cluster_results_2)

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	REAL(kind=DBL), INTENT (IN) :: velocity
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency
	REAL(kind=DBL), INTENT (INOUT) :: PDF_cluster_distr_2 (20,2)
	REAL(kind=DBL), INTENT (INOUT) :: cluster_results_2 (9, 1)

	INTEGER :: i, j, k, l , status
	INTEGER :: bubble_number
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord_swapper
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord_swapper

	INTEGER :: iptr, sum_up, counter, detector
	REAL(kind=DBL) :: temp

	REAL(kind=DBL) :: bubble_chord_sum
	REAL(kind=DBL) :: droplet_chord_sum
	REAL(kind=DBL) :: average_bubble_chord
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: instantenious_cluster
	REAL(kind=DBL) :: size
	REAL(kind=DBL) :: chord_size_cluster
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_cluster_no
	INTEGER :: cluster_container

	REAL(kind=DBL) ::  bubbles_in_cluster
	REAL(kind=DBL) ::  clusters_per_second
	REAL(kind=DBL) :: average_bubbles_cluster
	REAL(kind=DBL) :: average_chord_size_cluster
	REAL(kind=DBL) :: cluster_ratio
	REAL(kind=DBL) :: chord_lead_particle
	REAL(kind=DBL) :: chord_cluster
	REAL(kind=DBL) :: chord_average_cluster
	REAL(kind=DBL) :: chord_average_cluster_sum

	iptr = 0
	sum_up = 0
	counter = 0
	detector = 0
	temp = 0.0
	median_bubble_chord = 0.0
    median_droplet_chord = 0.0
	median_bubble_chord_crit = 0.0
	median_droplet_chord_crit = 0.0
	average_bubble_chord = 0.0
	size = 0.0
	chord_size_cluster = 0.0
	cluster_container = 0
	bubbles_in_cluster = 0.0
	clusters_per_second = 0.0
	average_bubbles_cluster = 0.0
	average_chord_size_cluster = 0.0
	cluster_ratio = 0.0
	chord_lead_particle = 0.0
	chord_cluster = 0.0
	chord_average_cluster = 0.0
	chord_average_cluster_sum = 0.0
	bubble_chord_sum = 0.0
	droplet_chord_sum = 0.0

	bubble_number = INT(basic_results(2,1)* REAL(sample_duration))
	
		ALLOCATE (bubble_chord(1, bubble_number))
		ALLOCATE (bubble_chord_swapper(1, bubble_number))	
		ALLOCATE (droplet_chord(1, (bubble_number-1)))
		ALLOCATE (droplet_chord_swapper(1, bubble_number -1))
		droplet_chord = 0.0
		droplet_chord_swapper = 0.0
		bubble_chord = 0.0
		bubble_chord_swapper = 0.0
		ALLOCATE (PDF_cluster_no (1, 20))
		PDF_cluster_no = 0

	! Air bubble clustering for C < 30%
IF (basic_results(1,1) < 0.3) THEN

		DO j = 1, bubble_number
				bubble_chord(1,j) = REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000*velocity
		END DO

		DO j = 1, bubble_number - 1
			droplet_chord(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
			droplet_chord_swapper(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
		END DO

	!calculate average bubble chord size
	DO j = 1, bubble_number
			bubble_chord_sum = bubble_chord_sum + REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000*velocity
	END DO
	average_bubble_chord = bubble_chord_sum / REAL(bubble_number)

	! calculate air bubble cluster statistics in bubble flow region
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0
	number = 0

	DO i = 1, (bubble_number -2)
		IF (droplet_chord(1,i) - bubble_chord(1,i) <= 0.00001) THEN
			sum_up = sum_up + 1
			size = size + bubble_chord(1, i)
			IF (droplet_chord(1,i) - bubble_chord(1,i) <= 0.00001 .AND. bubble_chord(1,i+1) - droplet_chord(1,i+1) <= 0.00001) THEN

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up - cluster_container <= 0.00001) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO


				detector = detector +1			! number of clusters
				counter = counter + sum_up	    ! number of bubbles in clusters
				chord_size_cluster = chord_size_cluster + (size + bubble_chord(1, i+1)) 		! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = bubble_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + bubble_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))


		sum_up = 1
		size = 0
			END IF
		END IF
	END DO

		DO k = 1, 20
		PDF_cluster_distr_2 (k,1) = k + 1
		PDF_cluster_distr_2 (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number) *100			! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results_2 (1, 1) = basic_results(1,1)
	cluster_results_2 (2, 1) = REAL(counter)
	cluster_results_2 (3, 1) = bubbles_in_cluster
	cluster_results_2 (4, 1) = REAL(detector)
	cluster_results_2 (5, 1) = clusters_per_second
	cluster_results_2 (6, 1) = average_bubbles_cluster
	cluster_results_2 (7, 1) = average_chord_size_cluster
	cluster_results_2 (8, 1) = cluster_ratio
	cluster_results_2 (9, 1) = chord_average_cluster

		
! Water droplet clustering for C > 70%	
ELSE IF (basic_results(1,1) > 0.7) THEN

		DO j = 1, (bubble_number-1)
				droplet_chord(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
		END DO

		DO j = 1, bubble_number
			bubble_chord(1,j) = REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000*velocity
			bubble_chord_swapper(1,j) = REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000*velocity
		END DO

	!calculate average bubble chord size
	DO j = 1, bubble_number - 1
			droplet_chord_sum = droplet_chord_sum + REAL(interfaces(1,j+1) - interfaces((2),j))/REAL(frequency)*1000*velocity
	END DO
	average_bubble_chord = droplet_chord_sum / REAL(bubble_number-1)


	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

		DO i = 1, (bubble_number -3)
		IF (bubble_chord(1,i+1) - droplet_chord(1,i) <= 0.00001) THEN
			sum_up = sum_up + 1
			size = size + droplet_chord(1, i)
			IF (bubble_chord(1,i+1) - droplet_chord(1,i) <= 0.00001 .AND. droplet_chord(1,i+1) - bubble_chord(1,i+2) <= 0.00001) THEN

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up - cluster_container <= 0.00001) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO


				detector = detector +1			! number of clusters
				counter = counter + sum_up	    ! number of bubbles in clusters
				chord_size_cluster = chord_size_cluster + (size + droplet_chord(1, i+1))			! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = droplet_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + droplet_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))


				sum_up = 1
				size = 0
			END IF
		END IF
	END DO

	DO k = 1, 20
		PDF_cluster_distr_2 (k,1) = k + 1
		PDF_cluster_distr_2 (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number-1) *100		! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results_2 (1, 1) = basic_results(1,1)
	cluster_results_2 (2, 1) = REAL(counter)
	cluster_results_2 (3, 1) = bubbles_in_cluster
	cluster_results_2 (4, 1) = REAL(detector)
	cluster_results_2 (5, 1) = clusters_per_second
	cluster_results_2 (6, 1) = average_bubbles_cluster
	cluster_results_2 (7, 1) = average_chord_size_cluster
	cluster_results_2 (8, 1) = cluster_ratio
	cluster_results_2 (9, 1) = chord_average_cluster

		DEALLOCATE(PDF_cluster_no, STAT = status)

END IF

		DEALLOCATE(bubble_chord, STAT = status)
		DEALLOCATE(bubble_chord_swapper, STAT = status)
		DEALLOCATE(droplet_chord, STAT = status)
		DEALLOCATE(droplet_chord_swapper, STAT = status)


END SUBROUTINE Cluster_analysis_wake


!____________________________________________________________________________________________________

! cluster analyses with chord times!

SUBROUTINE 	Cluster_analysis_time (interfaces, number_columns, max_no_bubbles, basic_results, sample_duration, &
					   frequency, cluster_factor, PDF_cluster_distr, cluster_results)

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency
	REAL, INTENT(IN) :: cluster_factor
	REAL(kind=DBL), INTENT (INOUT) :: PDF_cluster_distr (20,2)
	REAL(kind=DBL), INTENT (INOUT) :: cluster_results (9, 1)

		
	
!	REAL(kind=DBL), INTENT (INOUT) :: PDF_chord_size(INT(PDF_number_chord),3)
	INTEGER :: i, j, k, l , status
	INTEGER :: bubble_number
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord_swapper
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord_swapper

	INTEGER :: iptr, sum_up, counter, detector
	REAL(kind=DBL) :: temp

	REAL(kind=DBL) :: median_bubble_chord
	REAL(kind=DBL) :: median_droplet_chord
	REAL(kind=DBL) :: median_bubble_chord_crit
	REAL(kind=DBL) :: median_droplet_chord_crit
	REAL(kind=DBL) :: bubble_chord_sum
	REAL(kind=DBL) :: droplet_chord_sum
	REAL(kind=DBL) :: average_bubble_chord
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: instantenious_cluster
	REAL(kind=DBL) :: size
	REAL(kind=DBL) :: chord_size_cluster
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_cluster_no
	INTEGER :: cluster_container

	REAL(kind=DBL) ::  bubbles_in_cluster
	REAL(kind=DBL) ::  clusters_per_second
	REAL(kind=DBL) :: average_bubbles_cluster
	REAL(kind=DBL) :: average_chord_size_cluster
	REAL(kind=DBL) :: cluster_ratio
	REAL(kind=DBL) :: chord_lead_particle
	REAL(kind=DBL) :: chord_cluster
	REAL(kind=DBL) :: chord_average_cluster
	REAL(kind=DBL) :: chord_average_cluster_sum

	iptr = 0
	sum_up = 0
	counter = 0
	detector = 0
	temp = 0.0
	median_bubble_chord = 0.0
    median_droplet_chord = 0.0
	median_bubble_chord_crit = 0.0
	median_droplet_chord_crit = 0.0
	average_bubble_chord = 0.0
	size = 0.0
	chord_size_cluster = 0.0
	cluster_container = 0
	bubbles_in_cluster = 0.0
	clusters_per_second = 0.0
	average_bubbles_cluster = 0.0
	average_chord_size_cluster = 0.0
	cluster_ratio = 0.0
	chord_lead_particle = 0.0
	chord_cluster = 0.0
	chord_average_cluster = 0.0
	chord_average_cluster_sum = 0.0
	bubble_chord_sum = 0.0
	droplet_chord_sum = 0.0


	bubble_number = INT(basic_results(2,1)* REAL(sample_duration))
	
		ALLOCATE (bubble_chord(1, bubble_number))
		ALLOCATE (bubble_chord_swapper(1, bubble_number))	
		ALLOCATE (droplet_chord(1, (bubble_number-1)))
		ALLOCATE (droplet_chord_swapper(1, bubble_number -1))
		droplet_chord = 0.0
		droplet_chord_swapper = 0.0
		bubble_chord = 0.0
		bubble_chord_swapper = 0.0
		ALLOCATE (PDF_cluster_no (1, 20))
		PDF_cluster_no = 0

	! Air bubble clustering for C < 30%
IF (basic_results(1,1) < 0.3) THEN

		DO j = 1, bubble_number
				bubble_chord(1,j) = REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
		END DO

		DO j = 1, bubble_number - 1
			droplet_chord(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
			droplet_chord_swapper(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
		END DO


			!calculate median water chord
			DO i = 1, bubble_number-2
				iptr = i
				DO j = i+1, bubble_number - 1
					IF ( droplet_chord_swapper(1,j) < droplet_chord_swapper(1,iptr)) THEN
						iptr = j
					END IF
				END DO
			IF (i /= iptr) THEN 
			temp = droplet_chord_swapper(1,i)
			droplet_chord_swapper(1,i) = droplet_chord_swapper(1,iptr)
			droplet_chord_swapper(1,iptr) = temp
			END IF

			END DO

			IF ( mod((bubble_number - 1), 2) == 0) THEN 
				median_droplet_chord = ( droplet_chord_swapper (1,(bubble_number-1)/2) + droplet_chord_swapper (1,(bubble_number-1)/2+1) ) / 2
			ELSE
				median_droplet_chord = ( droplet_chord_swapper (1,(bubble_number-1)/2 + 1) )
			END IF


		!Calculate median droplet chord criterion value for cluster boundary criterion
			median_droplet_chord_crit = median_droplet_chord * cluster_factor

	!calculate average bubble chord size
	DO j = 1, bubble_number
			bubble_chord_sum = bubble_chord_sum + REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
	END DO
	average_bubble_chord = bubble_chord_sum / REAL(bubble_number)

	! calculate air bubble cluster statistics in bubble flow region
		ALLOCATE (instantenious_cluster(1, bubble_number - 1))
	DO i = 1, (bubble_number - 1)
		IF ( droplet_chord(1, i) <= median_droplet_chord_crit) THEN
			instantenious_cluster(1, i) = 1
		ELSE
			instantenious_cluster(1, i) = 0
		END IF
	END DO

	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

	DO i = 1, (bubble_number - 2)
		sum_up = sum_up + instantenious_cluster(1, i)
		IF (instantenious_cluster(1, i) == 1) THEN
		size = size + bubble_chord(1, i)
		END IF

		IF (instantenious_cluster(1, i) ==1 .AND. instantenious_cluster(1, i +1) == 0) THEN
		counter = counter + sum_up														! number of bubbles in clusters

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up <= cluster_container) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO

		detector = detector +1															! number of clusters

		chord_size_cluster = chord_size_cluster + (size + bubble_chord(1, i+1)) 			! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = bubble_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + bubble_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))
		sum_up = 1
		size = 0

		END IF

	END DO

	DEALLOCATE(instantenious_cluster, STAT = status)

	DO k = 1, 20
		PDF_cluster_distr (k,1) = k + 1
		PDF_cluster_distr (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number) *100		! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results (1, 1) = basic_results(1,1)
	cluster_results (2, 1) = REAL(counter)
	cluster_results (3, 1) = bubbles_in_cluster
	cluster_results (4, 1) = REAL(detector)
	cluster_results (5, 1) = clusters_per_second
	cluster_results (6, 1) = average_bubbles_cluster
	cluster_results (7, 1) = average_chord_size_cluster
	cluster_results (8, 1) = cluster_ratio
	cluster_results (9, 1) = chord_average_cluster


		
! Water droplet clustering for C > 70%	
ELSE IF (basic_results(1,1) > 0.7) THEN

		DO j = 1, (bubble_number-1)
				droplet_chord(1,j) = REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
		END DO

		DO j = 1, bubble_number
			bubble_chord(1,j) =  REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000
			bubble_chord_swapper(1,j) =  REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000
		END DO

	!calculate median air chord
			DO i = 1, bubble_number-1
				iptr = i
				DO j = i+1, bubble_number
					IF ( bubble_chord_swapper(1,j) < bubble_chord_swapper(1,iptr)) THEN
						iptr = j
					END IF
				END DO
			IF (i /= iptr) THEN 
			temp = bubble_chord_swapper(1,i)
			bubble_chord_swapper(1,i) = bubble_chord_swapper(1,iptr)
			bubble_chord_swapper(1,iptr) = temp
			END IF

			END DO

			IF ( mod(bubble_number, 2) == 0) THEN 
				median_bubble_chord = ( bubble_chord_swapper (1,bubble_number/2) + bubble_chord_swapper (1,bubble_number/2+1) ) / 2
			ELSE
				median_bubble_chord = ( bubble_chord_swapper (1,bubble_number/2 + 1) )
			END IF

		!Calculate median droplet chord criterion value for cluster boundary criterion
			median_bubble_chord_crit = median_bubble_chord * cluster_factor

	!calculate average bubble chord size
	DO j = 1, bubble_number - 1
			droplet_chord_sum = droplet_chord_sum +  REAL(interfaces(1,j+1) - interfaces((2),j))/REAL(frequency)*1000
	END DO
	average_bubble_chord = droplet_chord_sum / REAL(bubble_number-1)


	! calculate air bubble cluster statistics in bubble flow region
		ALLOCATE (instantenious_cluster(1, bubble_number))
	DO i = 1, (bubble_number)
		IF ( bubble_chord(1, i) <= median_bubble_chord_crit) THEN
			instantenious_cluster(1, i) = 1
		ELSE
			instantenious_cluster(1, i) = 0
		END IF
	END DO

	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

	DO i = 1, (bubble_number-2)
		sum_up = sum_up + instantenious_cluster(1, i)
		IF (instantenious_cluster(1, i) == 1) THEN
		size = size + droplet_chord(1, i)
		END IF

		IF (instantenious_cluster(1, i) ==1 .AND. instantenious_cluster(1, i +1) == 0) THEN
		counter = counter + sum_up														! number of bubbles in clusters

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up <= cluster_container) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO

		detector = detector +1															! number of clusters

		chord_size_cluster = chord_size_cluster + (size + droplet_chord(1, i+1))			! sum of bubble chord sizes of clusters
		
		chord_cluster = 0.0

		chord_lead_particle = droplet_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + droplet_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))
		
		sum_up = 1
		size = 0

		END IF

	END DO

		DEALLOCATE(instantenious_cluster, STAT = status)

	DO k = 1, 20
		PDF_cluster_distr (k,1) = k + 1
		PDF_cluster_distr (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number-1)	 *100		! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results (1, 1) = basic_results(1,1)
	cluster_results (2, 1) = REAL(counter)
	cluster_results (3, 1) = bubbles_in_cluster
	cluster_results (4, 1) = REAL(detector)
	cluster_results (5, 1) = clusters_per_second
	cluster_results (6, 1) = average_bubbles_cluster
	cluster_results (7, 1) = average_chord_size_cluster
	cluster_results (8, 1) = cluster_ratio
	cluster_results (9, 1) = chord_average_cluster

		DEALLOCATE(PDF_cluster_no, STAT = status)

END IF

		DEALLOCATE(bubble_chord, STAT = status)
		DEALLOCATE(bubble_chord_swapper, STAT = status)
		DEALLOCATE(droplet_chord, STAT = status)
		DEALLOCATE(droplet_chord_swapper, STAT = status)


END SUBROUTINE Cluster_analysis_time


!____________________________________________________________________________________________________

SUBROUTINE 	Cluster_analysis_1_time (interfaces, number_columns, max_no_bubbles, basic_results, sample_duration, &
					   frequency, cluster_factor_1, PDF_cluster_distr_1, cluster_results_1)

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency
	REAL, INTENT (IN) :: cluster_factor_1
	REAL(kind=DBL), INTENT (INOUT) :: PDF_cluster_distr_1 (20,2)
	REAL(kind=DBL), INTENT (INOUT) :: cluster_results_1 (9, 1)

		
	
!	REAL(kind=DBL), INTENT (INOUT) :: PDF_chord_size(INT(PDF_number_chord),3)
	INTEGER :: i, j, k, l , status
	INTEGER :: bubble_number
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord_swapper
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord_swapper

	INTEGER :: iptr, sum_up, counter, detector
	REAL(kind=DBL) :: temp

	REAL(kind=DBL) :: median_bubble_chord
	REAL(kind=DBL) :: median_droplet_chord
	REAL(kind=DBL) :: median_bubble_chord_crit
	REAL(kind=DBL) :: median_droplet_chord_crit
	REAL(kind=DBL) :: bubble_chord_sum
	REAL(kind=DBL) :: droplet_chord_sum
	REAL(kind=DBL) :: average_bubble_chord
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: instantenious_cluster
	REAL(kind=DBL) :: size
	REAL(kind=DBL) :: chord_size_cluster
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_cluster_no
	INTEGER :: cluster_container

	REAL(kind=DBL) ::  bubbles_in_cluster
	REAL(kind=DBL) ::  clusters_per_second
	REAL(kind=DBL) :: average_bubbles_cluster
	REAL(kind=DBL) :: average_chord_size_cluster
	REAL(kind=DBL) :: cluster_ratio
	REAL(kind=DBL) :: chord_lead_particle
	REAL(kind=DBL) :: chord_cluster
	REAL(kind=DBL) :: chord_average_cluster
	REAL(kind=DBL) :: chord_average_cluster_sum

	iptr = 0
	sum_up = 0
	counter = 0
	detector = 0
	temp = 0.0
	median_bubble_chord = 0.0
    median_droplet_chord = 0.0
	median_bubble_chord_crit = 0.0
	median_droplet_chord_crit = 0.0
	average_bubble_chord = 0.0
	size = 0.0
	chord_size_cluster = 0.0
	cluster_container = 0
	bubbles_in_cluster = 0.0
	clusters_per_second = 0.0
	average_bubbles_cluster = 0.0
	average_chord_size_cluster = 0.0
	cluster_ratio = 0.0
	chord_lead_particle = 0.0
	chord_cluster = 0.0
	chord_average_cluster = 0.0
	chord_average_cluster_sum = 0.0
	bubble_chord_sum = 0.0
	droplet_chord_sum = 0.0

	bubble_number = INT(basic_results(2,1)* REAL(sample_duration))
	
		ALLOCATE (bubble_chord(1, bubble_number))
		ALLOCATE (bubble_chord_swapper(1, bubble_number))	
		ALLOCATE (droplet_chord(1, (bubble_number-1)))
		ALLOCATE (droplet_chord_swapper(1, bubble_number -1))
		droplet_chord = 0.0
		droplet_chord_swapper = 0.0
		bubble_chord = 0.0
		bubble_chord_swapper = 0.0
		ALLOCATE (PDF_cluster_no (1, 20))
		PDF_cluster_no = 0

	! Air bubble clustering for C < 30%
IF (basic_results(1,1) < 0.3) THEN

		DO j = 1, bubble_number
				bubble_chord(1,j) =  REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
		END DO

		DO j = 1, bubble_number - 1
			droplet_chord(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
			droplet_chord_swapper(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
		END DO

	!calculate average bubble chord size
	DO j = 1, bubble_number
			bubble_chord_sum = bubble_chord_sum +  REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
	END DO
	average_bubble_chord = bubble_chord_sum / REAL(bubble_number)


	! calculate air bubble cluster statistics in bubble flow region
		ALLOCATE (instantenious_cluster(1, bubble_number - 1))
	DO i = 1, (bubble_number - 1)
		IF ( droplet_chord(1, i) <= (cluster_factor_1 * 1.0)) THEN
			instantenious_cluster(1, i) = 1
		ELSE
			instantenious_cluster(1, i) = 0
		END IF
	END DO

	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

	DO i = 1, (bubble_number - 2)
		sum_up = sum_up + instantenious_cluster(1, i)
		IF (instantenious_cluster(1, i) == 1) THEN
		size = size + bubble_chord(1, i)
		END IF

		IF (instantenious_cluster(1, i) ==1 .AND. instantenious_cluster(1, i +1) == 0) THEN
		counter = counter + sum_up														! number of bubbles in clusters

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up <= cluster_container) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO

		detector = detector +1															! number of clusters

		chord_size_cluster = chord_size_cluster + (size + bubble_chord(1, i+1))			! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = bubble_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + bubble_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))

		
		sum_up = 1
		size = 0

		END IF

	END DO

	DEALLOCATE(instantenious_cluster, STAT = status)

	DO k = 1, 20
		PDF_cluster_distr_1 (k,1) = k + 1
		PDF_cluster_distr_1 (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number) *100			! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results_1 (1, 1) = basic_results(1,1)
	cluster_results_1 (2, 1) = REAL(counter)
	cluster_results_1 (3, 1) = bubbles_in_cluster
	cluster_results_1 (4, 1) = REAL(detector)
	cluster_results_1 (5, 1) = clusters_per_second
	cluster_results_1 (6, 1) = average_bubbles_cluster
	cluster_results_1 (7, 1) = average_chord_size_cluster
	cluster_results_1 (8, 1) = cluster_ratio
	cluster_results_1 (9, 1) = chord_average_cluster

		
! Water droplet clustering for C > 70%	
ELSE IF (basic_results(1,1) > 0.7) THEN

		DO j = 1, (bubble_number-1)
				droplet_chord(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
		END DO

		DO j = 1, bubble_number
			bubble_chord(1,j) =  REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000
			bubble_chord_swapper(1,j) =  REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000
		END DO

	!calculate average bubble chord size
	DO j = 1, bubble_number - 1
			droplet_chord_sum = droplet_chord_sum +  REAL(interfaces(1,j+1) - interfaces((2),j))/REAL(frequency)*1000
	END DO
	average_bubble_chord = droplet_chord_sum / REAL(bubble_number-1)


	! calculate air bubble cluster statistics in bubble flow region
		ALLOCATE (instantenious_cluster(1, bubble_number))
	DO i = 1, (bubble_number)
		IF ( bubble_chord(1, i) <= (cluster_factor_1 * 1.0)) THEN
			instantenious_cluster(1, i) = 1
		ELSE
			instantenious_cluster(1, i) = 0
		END IF
	END DO

	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

	DO i = 1, (bubble_number-2)
		sum_up = sum_up + instantenious_cluster(1, i)
		IF (instantenious_cluster(1, i) == 1) THEN
		size = size + droplet_chord(1, i)
		END IF

		IF (instantenious_cluster(1, i) ==1 .AND. instantenious_cluster(1, i +1) == 0) THEN
		counter = counter + sum_up														! number of bubbles in clusters

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up <= cluster_container) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO

		detector = detector +1															! number of clusters

		chord_size_cluster = chord_size_cluster + (size + droplet_chord(1, i+1))			! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = droplet_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + droplet_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))

		
		sum_up = 1
		size = 0

		END IF

	END DO

		DEALLOCATE(instantenious_cluster, STAT = status)

	DO k = 1, 20
		PDF_cluster_distr_1 (k,1) = k + 1
		PDF_cluster_distr_1 (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number-1) *100			! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results_1 (1, 1) = basic_results(1,1)
	cluster_results_1 (2, 1) = REAL(counter)
	cluster_results_1 (3, 1) = bubbles_in_cluster
	cluster_results_1 (4, 1) = REAL(detector)
	cluster_results_1 (5, 1) = clusters_per_second
	cluster_results_1 (6, 1) = average_bubbles_cluster
	cluster_results_1 (7, 1) = average_chord_size_cluster
	cluster_results_1 (8, 1) = cluster_ratio
	cluster_results_1 (9, 1) = chord_average_cluster

		DEALLOCATE(PDF_cluster_no, STAT = status)

END IF

		DEALLOCATE(bubble_chord, STAT = status)
		DEALLOCATE(bubble_chord_swapper, STAT = status)
		DEALLOCATE(droplet_chord, STAT = status)
		DEALLOCATE(droplet_chord_swapper, STAT = status)


END SUBROUTINE Cluster_analysis_1_time


!____________________________________________________________________________________________________

SUBROUTINE 	Cluster_analysis_wake_time (interfaces, number_columns, max_no_bubbles, basic_results,  sample_duration, &
					   frequency, PDF_cluster_distr_2, cluster_results_2)

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency
	REAL(kind=DBL), INTENT (INOUT) :: PDF_cluster_distr_2 (20,2)
	REAL(kind=DBL), INTENT (INOUT) :: cluster_results_2 (9, 1)

	INTEGER :: i, j, k, l , status
	INTEGER :: bubble_number
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_chord_swapper
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_chord_swapper

	INTEGER :: iptr, sum_up, counter, detector
	REAL(kind=DBL) :: temp

	REAL(kind=DBL) :: bubble_chord_sum
	REAL(kind=DBL) :: droplet_chord_sum
	REAL(kind=DBL) :: average_bubble_chord
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: instantenious_cluster
	REAL(kind=DBL) :: size
	REAL(kind=DBL) :: chord_size_cluster
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: PDF_cluster_no
	INTEGER :: cluster_container

	REAL(kind=DBL) ::  bubbles_in_cluster
	REAL(kind=DBL) ::  clusters_per_second
	REAL(kind=DBL) :: average_bubbles_cluster
	REAL(kind=DBL) :: average_chord_size_cluster
	REAL(kind=DBL) :: cluster_ratio
	REAL(kind=DBL) :: chord_lead_particle
	REAL(kind=DBL) :: chord_cluster
	REAL(kind=DBL) :: chord_average_cluster
	REAL(kind=DBL) :: chord_average_cluster_sum

	iptr = 0
	sum_up = 0
	counter = 0
	detector = 0
	temp = 0.0
	median_bubble_chord = 0.0
    median_droplet_chord = 0.0
	median_bubble_chord_crit = 0.0
	median_droplet_chord_crit = 0.0
	average_bubble_chord = 0.0
	size = 0.0
	chord_size_cluster = 0.0
	cluster_container = 0
	bubbles_in_cluster = 0.0
	clusters_per_second = 0.0
	average_bubbles_cluster = 0.0
	average_chord_size_cluster = 0.0
	cluster_ratio = 0.0
	chord_lead_particle = 0.0
	chord_cluster = 0.0
	chord_average_cluster = 0.0
	chord_average_cluster_sum = 0.0
	bubble_chord_sum = 0.0
	droplet_chord_sum = 0.0

	bubble_number = INT(basic_results(2,1)* REAL(sample_duration))
	
		ALLOCATE (bubble_chord(1, bubble_number))
		ALLOCATE (bubble_chord_swapper(1, bubble_number))	
		ALLOCATE (droplet_chord(1, (bubble_number-1)))
		ALLOCATE (droplet_chord_swapper(1, bubble_number -1))
		droplet_chord = 0.0
		droplet_chord_swapper = 0.0
		bubble_chord = 0.0
		bubble_chord_swapper = 0.0
		ALLOCATE (PDF_cluster_no (1, 20))
		PDF_cluster_no = 0

	! Air bubble clustering for C < 30%
IF (basic_results(1,1) < 0.3) THEN

		DO j = 1, bubble_number
				bubble_chord(1,j) =  REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
		END DO

		DO j = 1, bubble_number - 1
			droplet_chord(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
			droplet_chord_swapper(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
		END DO

	!calculate average bubble chord size
	DO j = 1, bubble_number
			bubble_chord_sum = bubble_chord_sum +  REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
	END DO
	average_bubble_chord = bubble_chord_sum / REAL(bubble_number)

	! calculate air bubble cluster statistics in bubble flow region
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0
	number = 0

	DO i = 1, (bubble_number -2)
		IF (droplet_chord(1,i) - bubble_chord(1,i) <= 0.00001) THEN
			sum_up = sum_up + 1
			size = size + bubble_chord(1, i)
			IF (droplet_chord(1,i) - bubble_chord(1,i) <= 0.00001 .AND. bubble_chord(1,i+1) - droplet_chord(1,i+1) <= 0.00001) THEN

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up - cluster_container <= 0.00001) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO


				detector = detector +1			! number of clusters
				counter = counter + sum_up	    ! number of bubbles in clusters
				chord_size_cluster = chord_size_cluster + (size + bubble_chord(1, i+1))		! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = bubble_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + bubble_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))


		sum_up = 1
		size = 0
			END IF
		END IF
	END DO

		DO k = 1, 20
		PDF_cluster_distr_2 (k,1) = k + 1
		PDF_cluster_distr_2 (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number) *100			! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results_2 (1, 1) = basic_results(1,1)
	cluster_results_2 (2, 1) = REAL(counter)
	cluster_results_2 (3, 1) = bubbles_in_cluster
	cluster_results_2 (4, 1) = REAL(detector)
	cluster_results_2 (5, 1) = clusters_per_second
	cluster_results_2 (6, 1) = average_bubbles_cluster
	cluster_results_2 (7, 1) = average_chord_size_cluster
	cluster_results_2 (8, 1) = cluster_ratio
	cluster_results_2 (9, 1) = chord_average_cluster

		
! Water droplet clustering for C > 70%	
ELSE IF (basic_results(1,1) > 0.7) THEN

		DO j = 1, (bubble_number-1)
				droplet_chord(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
		END DO

		DO j = 1, bubble_number
			bubble_chord(1,j) =  REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000
			bubble_chord_swapper(1,j) =  REAL(interfaces(2,j) - interfaces(1,j))/REAL(frequency)*1000
		END DO

	!calculate average bubble chord size
	DO j = 1, bubble_number - 1
			droplet_chord_sum = droplet_chord_sum +  REAL(interfaces(1,j+1) - interfaces((2),j))/REAL(frequency)*1000
	END DO
	average_bubble_chord = droplet_chord_sum / REAL(bubble_number-1)


	! calculate the size of the bubble in the clusters and the number of bubble clusters
	sum_up = 1
	size = 0.0
	counter = 0
	detector = 0
	chord_size_cluster = 0.0

		DO i = 1, (bubble_number -3)
		IF (bubble_chord(1,i+1) - droplet_chord(1,i) <= 0.00001) THEN
			sum_up = sum_up + 1
			size = size + droplet_chord(1, i)
			IF (bubble_chord(1,i+1) - droplet_chord(1,i) <= 0.00001 .AND. droplet_chord(1,i+1) - bubble_chord(1,i+2) <= 0.00001) THEN

			!PDF of bubble sizes
			cluster_container = 0
			DO k = 1, 19
				cluster_container = k + 1
				IF (sum_up - cluster_container <= 0.00001) THEN
					PDF_cluster_no(1, k) = PDF_cluster_no(1, k) +1
				EXIT
				ELSE IF (sum_up > 20) THEN
				PDF_cluster_no(1, 20) = PDF_cluster_no(1, 20) +1
				EXIT
				END IF
			END DO


				detector = detector +1			! number of clusters
				counter = counter + sum_up	    ! number of bubbles in clusters
				chord_size_cluster = chord_size_cluster + (size + droplet_chord(1, i+1)) 			! sum of bubble chord sizes of clusters

		chord_cluster = 0.0

		chord_lead_particle = droplet_chord(1, i + 2 -sum_up)									! ratio leading partcile in cluster to cluster average chord size
		DO m = 1, sum_up
		chord_cluster = chord_cluster + droplet_chord(1, (i + 1 - sum_up +m))
		END DO
		chord_average_cluster_sum = chord_average_cluster_sum + chord_lead_particle / (chord_cluster/ REAL (sum_up))


				sum_up = 1
				size = 0
			END IF
		END IF
	END DO

	DO k = 1, 20
		PDF_cluster_distr_2 (k,1) = k + 1
		PDF_cluster_distr_2 (k,2) = REAL(PDF_cluster_no (1,k)) / REAL(detector)
	END DO

	bubbles_in_cluster = REAL(counter) / REAL(bubble_number-1) *100		! bubble in cluster in %

	clusters_per_second = REAL(detector) / REAL(sample_duration)	! clusters per second

	average_bubbles_cluster = REAL(counter) / REAL(detector)		! average number of bubbles per cluster

	average_chord_size_cluster = chord_size_cluster / REAL(counter)   ! average size of cluster chord size in mm

	cluster_ratio = average_chord_size_cluster / average_bubble_chord	! avarage cluster size compared to average bubble chord size

	chord_average_cluster = chord_average_cluster_sum / REAL(detector)	! average ratio of chord size of leading particle and average of cluster

	cluster_results_2 (1, 1) = basic_results(1,1)
	cluster_results_2 (2, 1) = REAL(counter)
	cluster_results_2 (3, 1) = bubbles_in_cluster
	cluster_results_2 (4, 1) = REAL(detector)
	cluster_results_2 (5, 1) = clusters_per_second
	cluster_results_2 (6, 1) = average_bubbles_cluster
	cluster_results_2 (7, 1) = average_chord_size_cluster
	cluster_results_2 (8, 1) = cluster_ratio
	cluster_results_2 (9, 1) = chord_average_cluster

		DEALLOCATE(PDF_cluster_no, STAT = status)

END IF

		DEALLOCATE(bubble_chord, STAT = status)
		DEALLOCATE(bubble_chord_swapper, STAT = status)
		DEALLOCATE(droplet_chord, STAT = status)
		DEALLOCATE(droplet_chord_swapper, STAT = status)


END SUBROUTINE Cluster_analysis_wake_time


!____________________________________________________________________________________________________


SUBROUTINE Particle_grouping (interfaces, number_columns, max_no_bubbles, basic_results, &
					   sample_duration, frequency, velocity, PDF_particle_group )

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency
	REAL(kind=DBL), INTENT(IN) :: velocity
	REAL(kind=DBL), INTENT (INOUT) :: PDF_particle_group (101,11)								

	INTEGER :: i, j, k, l , m, status
	INTEGER :: bubble_number
	REAL :: segmenter_time
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_time 
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_time
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_length
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_length 
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: number_in_group
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: flagger
	REAL(kind=DBL) :: particle_time
	INTEGER :: pointer_interfaces, 	counter_interfaces

	REAL(kind=DBL) :: PDF_helper (100,5)

	status = 0
	bubble_number = 0
	segmenter_time = 0.0
	particle_time = 0.0
	pointer_interfaces = 0
	counter_interfaces = 0
	PDF_helper = 0.0


	bubble_number = INT(basic_results(2,1)* REAL(sample_duration))

		ALLOCATE (bubble_time(1, INT(bubble_number)))
		ALLOCATE (bubble_length(1, INT(bubble_number)))
		ALLOCATE (droplet_time(1, (INT(bubble_number)-1)))
		ALLOCATE (droplet_length(1, (INT(bubble_number)-1)))
		bubble_time = 0.0
		bubble_length = 0.0	
		droplet_time = 0.0
		droplet_length = 0.0

		DO j = 1, bubble_number
			bubble_time(1,j) =  REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
			bubble_length(1,j) =  REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000*velocity
		END DO

		DO j = 1, (bubble_number-1)
			droplet_time(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000	
			droplet_length(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000*velocity
		END DO



! Air bubble clustering for C < 30%
IF (basic_results(1,1) < 0.3) THEN

		ALLOCATE (number_in_group(5, 1))
		ALLOCATE (flagger (1, bubble_number))
		number_in_group = 0
		flagger = 0
		PDF_helper = 0
		status = 0
		segmenter_time = 0.0

		DO j = 1, bubble_number
			IF(bubble_length(1,j) <= 1.0) THEN
				number_in_group(1,1) = number_in_group(1,1) +1
				flagger (1,j) = 1
				ELSE IF (bubble_length(1,j) <= 3.0) THEN
				number_in_group(2,1) = number_in_group(2,1) +1
				flagger (1,j) = 2
				ELSE IF (bubble_length(1,j) <= 6.0) THEN
				number_in_group(3,1) = number_in_group(3,1) +1
				flagger (1,j) = 3
				ELSE IF (bubble_length(1,j) <= 10.0) THEN
				number_in_group(4,1) = number_in_group(4,1) +1
				flagger (1,j) = 4
				ELSE 
				number_in_group(5,1) = number_in_group(5,1) +1
				flagger (1,j) = 5
			END IF		
		END DO

	particle_time = 0.0
	pointer_interfaces = 0
	
	DO k = 1, 5
		DO j = 1, bubble_number 
			IF (flagger(1,j) == k) THEN
				pointer_interfaces = interfaces(2,j)
				counter_interfaces = 0
				DO i = 1, (bubble_number - j)
					counter_interfaces =  i
						IF (flagger(1,j+i) == k) THEN
								particle_time =  REAL(interfaces (1, j+counter_interfaces) - pointer_interfaces) /REAL(frequency)*1000
								DO m =1, 100
									IF(particle_time <= m) THEN
									PDF_helper(m,k) = PDF_helper(m,k) +1
									EXIT
									ELSE IF (particle_time > 100) THEN
									PDF_helper(100,k) = PDF_helper(100,k) +1
									EXIT
									END IF
								END DO
						EXIT
						END IF
				END DO
			END IF
		END DO

	END DO

!		DO j = 1, bubble_number
!			write (*,*) flagger (1,j)
!		END DO
!			pause
					

		DEALLOCATE(bubble_time, STAT = status)
		segmenter_time =0.0

		DO i = 1, 100
			PDF_particle_group (i+1,1) = REAL(i)
		END DO

		DO j = 1, 5
			DO i = 1, 100
				PDF_particle_group (1,2*j) =   REAL(number_in_group(j,1))
				PDF_particle_group (1,2*j+1) =   REAL(number_in_group(j,1))
				PDF_particle_group (1,1) =  9999
				PDF_particle_group (i+1,2*j+1) = ( REAL(number_in_group(j,1))/REAL(sample_duration*1000) * REAL(number_in_group(j,1))/REAL(sample_duration*1000) * REAL(sample_duration *1000 -i) * exp (- REAL(number_in_group(j,1))/REAL(sample_duration *1000) *REAL(i)) ) / &
												( REAL(number_in_group(j,1)) - 1.0 + exp (- REAL(number_in_group(j,1) ) ) )

				PDF_particle_group (i+1,j*2) = PDF_helper (i,j) /  REAL(number_in_group(j,1))
			END DO
		END DO

		DEALLOCATE (number_in_group, STAT = status)
		DEALLOCATE (flagger , STAT = status)

	! Water droplet clustering for C > 70%	
	ELSE IF (basic_results(1,1) > 0.7) THEN

		ALLOCATE (number_in_group(5, 1))
		ALLOCATE (flagger (1, bubble_number-1))
		number_in_group = 0
		PDF_helper = 0
		status = 0
		segmenter_time = 0.0

		DO j = 1, bubble_number-1
			IF(droplet_length(1,j) <= 1.0) THEN
				number_in_group(1,1) = number_in_group(1,1) +1
				flagger (1,j) = 1
				ELSE IF (droplet_length(1,j) <= 3.0) THEN
				number_in_group(2,1) = number_in_group(2,1) +1
				flagger (1,j) = 2
				ELSE IF (droplet_length(1,j) <= 6.0) THEN
				number_in_group(3,1) = number_in_group(3,1) +1
				flagger (1,j) = 3
				ELSE IF (droplet_length(1,j) <= 10.0) THEN
				number_in_group(4,1) = number_in_group(4,1) +1
				flagger (1,j) = 4
				ELSE 
				number_in_group(5,1) = number_in_group(5,1) +1
				flagger (1,j) = 5
			END IF		
		END DO

	particle_time = 0.0
	pointer_interfaces = 0
	
	DO k = 1, 5
		DO j = 1, bubble_number -1
			IF (flagger(1,j) == k) THEN
				pointer_interfaces = interfaces(1,j+1)
				counter_interfaces = 0
				DO i = 1, (bubble_number -1 - j)
					counter_interfaces =  i
						IF (flagger(1,j+i) == k) THEN
								particle_time =  REAL(interfaces (2, j+counter_interfaces) - pointer_interfaces) /REAL(frequency)*1000
								DO m =1, 100
									IF(particle_time <= m) THEN
									PDF_helper(m,k) = PDF_helper(m,k) +1
									EXIT
									ELSE IF (particle_time > 100) THEN
									PDF_helper(100,k) = PDF_helper(100,k) +1
									EXIT
									END IF
								END DO
						EXIT
						END IF
				END DO
			END IF
		END DO

	END DO

!		DO j = 1, bubble_number
!			write (*,*) flagger (1,j)
!		END DO
!			pause
					

		DEALLOCATE(bubble_time, STAT = status)
		segmenter_time =0.0

		DO i = 1, 100
			PDF_particle_group (i+1,1) = REAL(i)
		END DO

		DO j = 1, 5
			DO i = 1, 100
				PDF_particle_group (1,2*j) =   REAL(number_in_group(j,1))
				PDF_particle_group (1,2*j+1) =   REAL(number_in_group(j,1))
				PDF_particle_group (1,1) =  9999
				PDF_particle_group (i+1,2*j+1) = ( REAL(number_in_group(j,1))/REAL(sample_duration*1000) * REAL(number_in_group(j,1))/REAL(sample_duration*1000) * REAL(sample_duration *1000 -i) * exp (- REAL(number_in_group(j,1))/REAL(sample_duration *1000) *REAL(i)) ) / &
												( REAL(number_in_group(j,1)) - 1.0 + exp (- REAL(number_in_group(j,1) ) ) )

				PDF_particle_group (i+1,j*2) = PDF_helper (i,j) /  REAL(number_in_group(j,1))
			END DO
		END DO

		DEALLOCATE (number_in_group, STAT = status)
		DEALLOCATE (flagger , STAT = status)


END IF



END SUBROUTINE Particle_grouping



!____________________________________________________________________________________________________


SUBROUTINE Particle_grouping_time (interfaces, number_columns, max_no_bubbles, basic_results, &
					   sample_duration, frequency, PDF_particle_group_time )

	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)	!Precision of the Real variables
	INTEGER, INTENT (IN) :: interfaces((number_columns-1)*2, max_no_bubbles)
	INTEGER, INTENT (IN) :: number_columns
	INTEGER, INTENT (IN) :: max_no_bubbles
	REAL(kind=DBL), INTENT (IN) :: basic_results((number_columns-1)*2,1)
	INTEGER, INTENT(IN) :: sample_duration
	INTEGER, INTENT(IN) :: frequency
	REAL(kind=DBL), INTENT (INOUT) :: PDF_particle_group_time (101,11)								

	INTEGER :: i, j, k, l , m, status
	INTEGER :: bubble_number
	REAL :: segmenter_time
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_time 
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_time
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: bubble_length
	REAL(kind=DBL), ALLOCATABLE, DIMENSION  (:,:) :: droplet_length 
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: number_in_group
	INTEGER, ALLOCATABLE, DIMENSION  (:,:) :: flagger
	REAL(kind=DBL) :: particle_time
	INTEGER :: pointer_interfaces, 	counter_interfaces

	REAL(kind=DBL) :: PDF_helper (100,5)

	status = 0
	bubble_number = 0
	segmenter_time = 0.0
	particle_time = 0.0
	pointer_interfaces = 0
	counter_interfaces = 0
	PDF_helper = 0.0

	bubble_number = INT(basic_results(2,1)* REAL(sample_duration))

		ALLOCATE (bubble_time(1, INT(bubble_number)))
		ALLOCATE (bubble_length(1, INT(bubble_number)))
		ALLOCATE (droplet_time(1, (INT(bubble_number)-1)))
		ALLOCATE (droplet_length(1, (INT(bubble_number)-1)))
		bubble_time = 0.0
		bubble_length = 0.0	
		droplet_time = 0.0
		droplet_length = 0.0

		DO j = 1, bubble_number
			bubble_time(1,j) =  REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
			bubble_length(1,j) =  REAL(interfaces(2,j) - interfaces((1),j))/REAL(frequency)*1000
		END DO

		DO j = 1, (bubble_number-1)
			droplet_time(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000	
			droplet_length(1,j) =  REAL(interfaces(1,j+1) - interfaces(2,j))/REAL(frequency)*1000
		END DO



! Air bubble clustering for C < 30%
IF (basic_results(1,1) < 0.3) THEN

		ALLOCATE (number_in_group(5, 1))
		ALLOCATE (flagger (1, bubble_number))
		number_in_group = 0
		PDF_helper = 0
		status = 0
		segmenter_time = 0.0

		DO j = 1, bubble_number
			IF(bubble_length(1,j) <= 0.5) THEN
				number_in_group(1,1) = number_in_group(1,1) +1
				flagger (1,j) = 1
				ELSE IF (bubble_length(1,j) <= 1.5) THEN
				number_in_group(2,1) = number_in_group(2,1) +1
				flagger (1,j) = 2
				ELSE IF (bubble_length(1,j) <= 3.0) THEN
				number_in_group(3,1) = number_in_group(3,1) +1
				flagger (1,j) = 3
				ELSE IF (bubble_length(1,j) <= 5.0) THEN
				number_in_group(4,1) = number_in_group(4,1) +1
				flagger (1,j) = 4
				ELSE 
				number_in_group(5,1) = number_in_group(5,1) +1
				flagger (1,j) = 5
			END IF		
		END DO

	particle_time = 0.0
	pointer_interfaces = 0
	
	DO k = 1, 5
		DO j = 1, bubble_number 
			IF (flagger(1,j) == k) THEN
				pointer_interfaces = interfaces(2,j)
				counter_interfaces = 0
				DO i = 1, (bubble_number - j)
					counter_interfaces =  i
						IF (flagger(1,j+i) == k) THEN
								particle_time =  REAL(interfaces (1, j+counter_interfaces) - pointer_interfaces) /REAL(frequency)*1000
								DO m =1, 100
									IF(particle_time <= m) THEN
									PDF_helper(m,k) = PDF_helper(m,k) +1
									EXIT
									ELSE IF (particle_time > 100) THEN
									PDF_helper(100,k) = PDF_helper(100,k) +1
									EXIT
									END IF
								END DO
						EXIT
						END IF
				END DO
			END IF
		END DO

	END DO

!		DO j = 1, bubble_number
!			write (*,*) flagger (1,j)
!		END DO
!			pause
					

		DEALLOCATE(bubble_time, STAT = status)
		segmenter_time =0.0

		DO i = 1, 100
			PDF_particle_group_time (i+1,1) = REAL(i)
		END DO

		DO j = 1, 5
			DO i = 1, 100
				PDF_particle_group_time (1,2*j) =   REAL(number_in_group(j,1))
				PDF_particle_group_time (1,2*j+1) =   REAL(number_in_group(j,1))
				PDF_particle_group_time (1,1) =  9999
				PDF_particle_group_time (i+1,2*j+1) = ( REAL(number_in_group(j,1))/REAL(sample_duration*1000) * REAL(number_in_group(j,1))/REAL(sample_duration*1000) * REAL(sample_duration *1000 -i) * exp (- REAL(number_in_group(j,1))/REAL(sample_duration *1000) *REAL(i)) ) / &
												( REAL(number_in_group(j,1)) - 1.0 + exp (- REAL(number_in_group(j,1) ) ) )

				PDF_particle_group_time (i+1,j*2) = PDF_helper (i,j) /  REAL(number_in_group(j,1))
			END DO
		END DO

		DEALLOCATE (number_in_group, STAT = status)
		DEALLOCATE (flagger , STAT = status)

	! Water droplet clustering for C > 70%	
	ELSE IF (basic_results(1,1) > 0.7) THEN

		ALLOCATE (number_in_group(5, 1))
		ALLOCATE (flagger (1, bubble_number-1))
		number_in_group = 0
		PDF_helper = 0
		status = 0
		segmenter_time = 0.0

		DO j = 1, bubble_number-1
			IF(droplet_length(1,j) <= 0.5) THEN
				number_in_group(1,1) = number_in_group(1,1) +1
				flagger (1,j) = 1
				ELSE IF (droplet_length(1,j) <= 1.5) THEN
				number_in_group(2,1) = number_in_group(2,1) +1
				flagger (1,j) = 2
				ELSE IF (droplet_length(1,j) <= 3.0) THEN
				number_in_group(3,1) = number_in_group(3,1) +1
				flagger (1,j) = 3
				ELSE IF (droplet_length(1,j) <= 6.0) THEN
				number_in_group(4,1) = number_in_group(4,1) +1
				flagger (1,j) = 4
				ELSE 
				number_in_group(5,1) = number_in_group(5,1) +1
				flagger (1,j) = 5
			END IF		
		END DO

	particle_time = 0.0
	pointer_interfaces = 0
	
	DO k = 1, 5
		DO j = 1, bubble_number -1
			IF (flagger(1,j) == k) THEN
				pointer_interfaces = interfaces(1,j+1)
				counter_interfaces = 0
				DO i = 1, (bubble_number -1 - j)
					counter_interfaces =  i
						IF (flagger(1,j+i) == k) THEN
								particle_time = REAL(interfaces (2, j+counter_interfaces) - pointer_interfaces) /REAL(frequency)*1000.0
								DO m =1, 100
									IF(particle_time <= m) THEN
									PDF_helper(m,k) = PDF_helper(m,k) +1
									EXIT
									ELSE IF (particle_time > 100) THEN
									PDF_helper(100,k) = PDF_helper(100,k) +1
									EXIT
									END IF
								END DO
						EXIT
						END IF
				END DO
			END IF
		END DO

	END DO

!		DO j = 1, bubble_number
!			write (*,*) flagger (1,j)
!		END DO
!			pause
					

		DEALLOCATE(bubble_time, STAT = status)
		segmenter_time =0.0

		DO i = 1, 100
			PDF_particle_group_time (i+1,1) = REAL(i)
		END DO

		DO j = 1, 5
			DO i = 1, 100
				PDF_particle_group_time (1,2*j) =  REAL(number_in_group(j,1))
				PDF_particle_group_time (1,2*j+1) =   REAL(number_in_group(j,1))
				PDF_particle_group_time (1,1) =  9999
				PDF_particle_group_time (i+1,2*j+1) = ( REAL(number_in_group(j,1))/REAL(sample_duration*1000) * REAL(number_in_group(j,1))/REAL(sample_duration*1000) * REAL(sample_duration *1000 -i) * exp (- REAL(number_in_group(j,1))/REAL(sample_duration *1000) *REAL(i)) ) / &
												( REAL(number_in_group(j,1)) - 1.0 + exp (- REAL(number_in_group(j,1) ) ) )

				PDF_particle_group_time (i+1,j*2) = PDF_helper (i,j) /  REAL(number_in_group(j,1))
			END DO
		END DO

		DEALLOCATE (number_in_group, STAT = status)
		DEALLOCATE (flagger , STAT = status)


END IF



END SUBROUTINE Particle_grouping_time

