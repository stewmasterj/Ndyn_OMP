        7812 # max ourNpts across all images
   65.3219986       65.3219986       65.3219986     # global BoxSize
           1        2916        6912           7 # imID, myNpts, ourNpts, NNcells
    2    3    4    5    6    7    8 # NeighCellIDs
 bulkminmax:           1         666
           7 # NumPermutationsSendCells
sendmin:       667       838      1108      1378      1801      2071      2494
sendmax:       837      1107      1377      1800      2070      2493      2916
          19 # NumPermutationsGetCells
getmin:      2917      3064      3295      3526      3889      4036      4267      4498      4861      4987      5185      5332      5563      5794      6157      6283      6481      6607      6805
getmax:      3063      3294      3525      3888      4035      4266      4497      4860      4986      5184      5331      5562      5793      6156      6282      6480      6606      6804      6912
sourceCellIDs:    2    2    2    2    3    3    3    3    4    4    5    5    5    5    6    6    7    7    8
sourcemin:       727       874      1105      1336       727       874      1861      2092       793       919       727      1168      1861      2554       793      1171       793      1927       865
sourcemax:       873      1104      1335      1698       873      1104      2091      2454       918      1116       873      1398      2091      2916       918      1368       918      2124       972
           2        2916        7200           7 # imID, myNpts, ourNpts, NNcells
    1    3    4    5    6    7    8 # NeighCellIDs
 bulkminmax:           1         726
           7 # NumPermutationsSendCells
sendmin:       727       874      1105      1336      1699      1993      2455
sendmax:       873      1104      1335      1698      1992      2454      2916
          19 # NumPermutationsGetCells
getmin:      2917      3088      3358      3628      4051      4198      4429      4555      4753      5005      5401      5548      5779      5905      6103      6355      6751      6877      6985
getmax:      3087      3357      3627      4050      4197      4428      4554      4752      5004      5400      5547      5778      5904      6102      6354      6750      6876      6984      7200
sourceCellIDs:    1    1    1    1    3    3    4    4    4    4    5    5    6    6    6    6    7    8    8
sourcemin:       667       838      1108      1378       727       874       793       919      1765      2017       727      1168       793      1171      1765      2521       793       865      1837
sourcemax:       837      1107      1377      1800       873      1104       918      1116      2016      2412       873      1398       918      1368      2016      2916       918       972      2052
           3        2916        7200           7 # imID, myNpts, ourNpts, NNcells
    1    2    4    5    6    7    8 # NeighCellIDs
 bulkminmax:           1         726
           7 # NumPermutationsSendCells
sendmin:       727       874      1105      1399      1861      2092      2455
sendmax:       873      1104      1398      1860      2091      2454      2916
          19 # NumPermutationsGetCells
getmin:      2917      3088      3358      3628      4051      4198      4429      4555      4753      5005      5401      5548      5779      5905      6031      6283      6481      6877      6985
getmax:      3087      3357      3627      4050      4197      4428      4554      4752      5004      5400      5547      5778      5904      6030      6282      6480      6876      6984      7200
sourceCellIDs:    1    1    1    1    2    2    4    4    4    4    5    5    6    7    7    7    7    8    8
sourcemin:       667       838      1801      2071       727       874       793       919      1117      1369       727      1861       793       793      1171      1927      2521       865      1189
sourcemax:       837      1107      2070      2493       873      1104       918      1116      1368      1764       873      2091       918       918      1422      2124      2916       972      1404
           4        2916        7500           7 # imID, myNpts, ourNpts, NNcells
    1    2    3    5    6    7    8 # NeighCellIDs
 bulkminmax:           1         792
           7 # NumPermutationsSendCells
sendmin:       793       919      1117      1369      1765      2017      2413
sendmax:       918      1116      1368      1764      2016      2412      2916
          19 # NumPermutationsGetCells
getmin:      2917      3088      3358      3505      3736      4030      4492      4639      4870      5164      5626      5773      5899      6151      6277      6529      6637      6853      7069
getmax:      3087      3357      3504      3735      4029      4491      4638      4869      5163      5625      5772      5898      6150      6276      6528      6636      6852      7068      7500
sourceCellIDs:    1    1    2    2    2    2    3    3    3    3    5    6    6    7    7    8    8    8    8
sourcemin:       667       838       727       874      1699      1993       727       874      1105      1399       727       793      1765       793      1171       865      1189      1837      2485
sourcemax:       837      1107       873      1104      1992      2454       873      1104      1398      1860       873       918      2016       918      1422       972      1404      2052      2916
           5        2916        7200           7 # imID, myNpts, ourNpts, NNcells
    1    2    3    4    6    7    8 # NeighCellIDs
 bulkminmax:           1         726
           7 # NumPermutationsSendCells
sendmin:       727       874      1168      1399      1861      2092      2554
sendmax:       873      1167      1398      1860      2091      2553      2916
          19 # NumPermutationsGetCells
getmin:      2917      3088      3358      3628      4051      4198      4429      4576      4807      4933      5059      5311      5509      5905      6031      6283      6481      6877      6985
getmax:      3087      3357      3627      4050      4197      4428      4575      4806      4932      5058      5310      5508      5904      6030      6282      6480      6876      6984      7200
sourceCellIDs:    1    1    1    1    2    2    3    3    4    6    6    6    6    7    7    7    7    8    8
sourcemin:       667      1108      1801      2494       727      1105       727      1861       793       793       919      1171      1369       793       919      1927      2125       865       973
sourcemax:       837      1377      2070      2916       873      1335       873      2091       918       918      1170      1368      1764       918      1170      2124      2520       972      1188
           6        2916        7500           7 # imID, myNpts, ourNpts, NNcells
    1    2    3    4    5    7    8 # NeighCellIDs
 bulkminmax:           1         792
           7 # NumPermutationsSendCells
sendmin:       793       919      1171      1369      1765      2017      2521
sendmax:       918      1170      1368      1764      2016      2520      2916
          19 # NumPermutationsGetCells
getmin:      2917      3088      3358      3505      3736      4030      4492      4639      4765      5017      5164      5458      5689      6151      6277      6529      6637      6853      7069
getmax:      3087      3357      3504      3735      4029      4491      4638      4764      5016      5163      5457      5688      6150      6276      6528      6636      6852      7068      7500
sourceCellIDs:    1    1    2    2    2    2    3    4    4    5    5    5    5    7    7    8    8    8    8
sourcemin:       667      1108       727      1105      1699      2455       727       793      1765       727       874      1168      1399       793       919       865       973      1837      2053
sourcemax:       837      1377       873      1335      1992      2916       873       918      2016       873      1167      1398      1860       918      1170       972      1188      2052      2484
           7        2916        7500           7 # imID, myNpts, ourNpts, NNcells
    1    2    3    4    5    6    8 # NeighCellIDs
 bulkminmax:           1         792
           7 # NumPermutationsSendCells
sendmin:       793       919      1171      1423      1927      2125      2521
sendmax:       918      1170      1422      1926      2124      2520      2916
          19 # NumPermutationsGetCells
getmin:      2917      3088      3358      3505      3652      3946      4177      4639      4765      5017      5164      5458      5689      6151      6277      6529      6637      6853      7069
getmax:      3087      3357      3504      3651      3945      4176      4638      4764      5016      5163      5457      5688      6150      6276      6528      6636      6852      7068      7500
sourceCellIDs:    1    1    2    3    3    3    3    4    4    5    5    5    5    6    6    8    8    8    8
sourcemin:       667      1801       727       727      1105      1861      2455       793      1117       727       874      1861      2092       793       919       865       973      1189      1405
sourcemax:       837      2070       873       873      1398      2091      2916       918      1368       873      1167      2091      2553       918      1170       972      1188      1404      1836
           8        2916        7812           7 # imID, myNpts, ourNpts, NNcells
    1    2    3    4    5    6    7 # NeighCellIDs
 bulkminmax:           1         864
           7 # NumPermutationsSendCells
sendmin:       865       973      1189      1405      1837      2053      2485
sendmax:       972      1188      1404      1836      2052      2484      2916
          19 # NumPermutationsGetCells
getmin:      2917      3088      3235      3529      3676      3970      4096      4348      4600      5104      5251      5545      5671      5923      6175      6679      6805      7057      7309
getmax:      3087      3234      3528      3675      3969      4095      4347      4599      5103      5250      5544      5670      5922      6174      6678      6804      7056      7308      7812
sourceCellIDs:    1    2    2    3    3    4    4    4    4    5    5    6    6    6    6    7    7    7    7
sourcemin:       667       727      1699       727      1105       793      1117      1765      2413       727       874       793       919      1765      2017       793       919      1171      1423
sourcemax:       837       873      1992       873      1398       918      1368      2016      2916       873      1167       918      1170      2016      2520       918      1170      1422      1926
