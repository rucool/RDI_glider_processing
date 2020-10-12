# RDI_glider_processing
 Work space for processing glider mounted RDI DVL current data.

Current workflow (10/12/20):
(1) Run pd0 files through read_ADCP_PD0_v3.m
	- This step converts the files output by the instrument into .mat files readable by matlab

(2) Work step by step through sampling_strategy_check.m
	- The original purpose of this code was to make sure that we were actually returning good data during a glider deployment based on our settings.
	- There are important pieces in here that can be pulled out and built into a simpler/ more direct function for processing this data.
	- In this code is where the least square inversion function is used to parse out the measured currents into water velocity and glider velocity.
