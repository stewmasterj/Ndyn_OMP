       16200 # max ourNpts across all images
   65.3219986       65.3219986       65.3219986     # global BoxSize
           1       11664       15552           1 # imID, myNpts, ourNpts, NNcells
    2 # NeighCellIDs
 bulkminmax:           1        7128
           1 # NumPermutationsSendCells
sendmin:      7129
sendmax:     11664
           1 # NumPermutationsGetCells
getmin:     11665
getmax:     15552
sourceCellIDs:    2
sourcemin:      7777
sourcemax:     11664
           2       11664       16200           1 # imID, myNpts, ourNpts, NNcells
    1 # NeighCellIDs
 bulkminmax:           1        7776
           1 # NumPermutationsSendCells
sendmin:      7777
sendmax:     11664
           1 # NumPermutationsGetCells
getmin:     11665
getmax:     16200
sourceCellIDs:    1
sourcemin:      7129
sourcemax:     11664
