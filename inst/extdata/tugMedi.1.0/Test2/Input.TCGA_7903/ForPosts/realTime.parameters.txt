#Parameter	Value
input.cloneout	'./Output/cloneout.txt'
out	'./Output/cloneout.realTime.txt'

# either VDT or realTumorSize is necessary.
realTumorSize$LD1.mm	7
realTumorSize$SD1.mm	7
realTumorSize$LD2.mm	9
realTumorSize$SD2.mm	9
realTumorSize$Ntype	'N_primary'	# cloneout column name
realTumorSize$deltaT.days	390

outUnit	'year'
col.replace	T

