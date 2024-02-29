'Start copying here 
'This program creates a S21 measurement 
'It is written in VBscript using SCPI commands 

Dim app 
Dim scpi
Dim oFSO

' required for creating a folder
'Create / Get the VNA application 
Set app = CreateObject ("AgilentPNA835x.Application") 
Set scpi = app.ScpiStringParser 
Set oFSO = CreateObject("Scripting.FileSystemObject")

'Create / Get the VNA application 
Set app = CreateObject ("AgilentPNA835x.Application") 
Set scpi = app.ScpiStringParser 


'Define frequency etc
freqs = Array(3309.26, 3349.67, 3375.12, 3858.50, 3878.07, 3901.09, 3921.76, 3942.26, 3970.87, 4486.99, 4487.64, 4523.74, 4579.45, 4591.61, 4631.35, 4632.25, 5194.70, 5215.74, 5218.66, 5221.21, 5222.27, 5224.60, 5232.79, 5232.96, 5241.25, 5353.62, 5358.84, 5359.36)
powers = Array(-60, -60, -60, -55, -55, -55, -55, -55, -50, -55, -55, -55, -55, -55, -55, -55, -55, -60, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55)
spans = Array(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
Name_Var = "Optimum Power Closed Window"
power_change = 0
Folder_Name = CSTR(Year(Now)) + "_" + CSTR(Month(Now)) + "_" + CSTR(Day(Now)) + "_" + Name_Var
oFSO.CreateFolder "C:\Users\Instrument\My Documents\Frequency_Sweeps\" + Folder_Name


For i = 0 to UBound(freqs) Step 1
				
	freq = freqs(i)*1000000
	power = powers(i) + power_change
	span = spans(i)*1000000
	freqstr = CSTR(freq)
	powerstr = CSTR(power)
	spanstr = CSTR(span)

	'Preset the Analyzer.FPREset presets the setting and deletes all traces and windows 
	scpi.Execute ("SYST:FPReset") 
	'Create and turn on window 1 
	scpi.Execute ("DISPlay:WINDow1:STATE ON") 

	'Define a measurement name, parameter 
	scpi.Execute ("CALCulate:PARameter:DEFine:EXT 'MyMeas', 's21'") 

	'Set power level 
	scpi.Execute("SOUR:POW1 " + powerstr)

	'Set centre frequency of sweep
	scpi.Execute("SENS:FREQ:CENT " + freqstr)

	'Set frequency span of sweep
	scpi.Execute("SENS:FREQ:SPAN " + spanstr)

	'Set number of points in sweep
	scpi.Execute("SENS:SWE:POIN 10001")

	'Change format to Polar
	scpi.Execute("CALC:MEAS:FORM POL")

	'Associate ("FEED") the measurement name ('MyMeas') to WINDow (1), and give the new TRACE a number (1). 
	scpi.Execute ("DISPlay:WINDow1:TRACe1:FEED 'MyMeas'")

	'Autoscales display
	scpi.Execute ("DISP:MEAS:Y:AUTO")

	'Starts averaging
	scpi.Execute ("SENS:AVER ON")

	'Set number of points to average over
	scpi.Execute ("SENS:AVER:COUN 64")

	WScript.Sleep(10000) 'sleeps fpr 10 second to allow averaging to finish

	'Dim measurementdate as String
	filename = "C:\Users\Instrument\My Documents\Frequency_Sweeps\" + Folder_Name + "\" + "sweep_" + Name_Var + "_" + CSTR(Year(Now)) + "_" + CSTR(Month(Now)) + "_" + CSTR(Day(Now)) + "_" + powerstr + "_" + freqstr + ".csv"

	savestring = "'" + filename + "','CSV Formatted Data','displayed','RI',-1" 


	'Save data to csv file
	'scpi.Execute ("MMEM:STOR:DATA 'fsweep.csv','CSV Formatted Data','displayed','RI',-1")
	scpi.Execute ("MMEM:STOR:DATA " + savestring)

Next
 
