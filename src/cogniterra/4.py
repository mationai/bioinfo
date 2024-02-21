from superpipe import pipes
import sys
sys.path.append('src')
from dataexpected import expected
from extensions import *
from helpers import *
from peptidelib import *

@pipes
class chap4:
   def _1():
      ptn = translateProtein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
      ptn >> testVs('MAMAPRTEINSTRING')

      l = readlines('dataset_29910_4a.txt')
      l[0] >> translateProtein >> testVs(expected._4_1a)

      l = readlines('dataset_29910_4.txt')
      l[0] >> translateProtein >> testVs(expected._4_1)

      ptn = abreviate('Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr')
      ptn >> testVs('VKLFPWFNQY')
      cnt = permCountsOf(ptn)
      cnt >> testVs(6144)

   def _2():
      ptns = encodePeptide('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', 'MA')
      ptns >> testVs('ATGGCC GGCCAT ATGGCC'.split(' '))
      l = readlines('dataset_29910_7a.txt')
      encodePeptide(l[0], l[1]) >> sepWithNL >> testVs(expected._4_2a)
      l = readlines('dataset_29910_7.txt')
      encodePeptide(l[0], l[1]) >> sepWithNL >> testVs(expected._4_2)

   def _3():
      cycloSubpeptidesCount(1) >> testVs(0)
      cycloSubpeptidesCount(4) >> testVs(12)
      cycloSubpeptidesCount(31315) >> testVs(980597910)
      cycloSubpeptidesCount(49309) >> testVs(2431328172)

   def _4():
      linearSpectrum('NQEL') >> sepWithSp >> testVs('0 113 114 128 129 242 242 257 370 371 484')
      linearSpectrum('LEQN') >> sepWithSp >> testVs('0 113 114 128 129 242 242 257 370 371 484')

      wrapSubpeptides('TRYCV') >> sepWithSp >> testVs('YCVT CVT VTR VTRY VT')
      cycloSpectrum('LEQN') >> sepWithSp >> testVs('0 113 114 128 129 227 242 242 257 355 356 370 371 484')
      cycloSpectrum('IAQMLFYCKVATN') >> sepWithSp >> testVs('0 71 71 99 101 103 113 113 114 128 128 131 147 163 170 172 184 199 215 227 227 231 244 259 260 266 271 286 298 298 310 312 328 330 330 372 385 391 394 399 399 399 401 413 423 426 443 443 470 493 498 502 513 519 526 527 541 554 556 557 564 569 590 598 616 626 640 654 657 658 665 670 682 697 697 703 711 729 729 753 753 771 779 785 785 800 812 817 824 825 828 842 856 866 884 892 913 918 925 926 928 941 955 956 963 969 980 984 989 1012 1039 1039 1056 1059 1069 1081 1083 1083 1083 1088 1091 1097 1110 1152 1152 1154 1170 1172 1184 1184 1196 1211 1216 1222 1223 1238 1251 1255 1255 1267 1283 1298 1310 1312 1319 1335 1351 1354 1354 1368 1369 1369 1379 1381 1383 1411 1411 1482')
      cycloSpectrum('WAANGWRRGVCNEIR') >> sepWithSp >> testVs('0 57 57 71 71 99 103 113 114 114 129 142 156 156 156 156 171 185 186 186 202 213 217 242 242 243 243 256 257 259 269 312 312 313 316 328 342 342 346 356 357 369 373 398 399 413 415 428 442 445 455 459 468 484 498 499 499 502 512 513 526 529 555 555 558 571 584 584 597 598 612 615 615 654 655 655 655 658 669 685 685 698 711 711 714 726 726 740 757 768 769 771 771 797 801 811 814 814 825 840 840 841 841 868 871 872 896 897 900 927 927 928 928 943 954 954 957 967 971 997 997 999 1000 1011 1028 1042 1042 1054 1057 1057 1070 1083 1083 1099 1110 1113 1113 1113 1114 1153 1153 1156 1170 1171 1184 1184 1197 1210 1213 1213 1239 1242 1255 1256 1266 1269 1269 1270 1284 1300 1309 1313 1323 1326 1340 1353 1355 1369 1370 1395 1399 1411 1412 1422 1426 1426 1440 1452 1455 1456 1456 1499 1509 1511 1512 1525 1525 1526 1526 1551 1555 1566 1582 1582 1583 1597 1612 1612 1612 1612 1626 1639 1654 1654 1655 1665 1669 1697 1697 1711 1711 1768')
      cycloSpectrum('VLDTIVLHWMRMRDC') >> sepWithSp >> testVs('0 99 99 101 103 113 113 113 115 115 131 131 137 156 156 186 202 212 212 212 214 216 218 228 250 271 287 287 287 313 315 317 317 323 325 327 329 329 349 374 402 418 426 428 428 430 430 436 442 443 454 462 473 473 505 531 535 541 541 541 545 558 563 567 574 586 604 604 610 640 644 646 648 654 661 666 678 689 701 717 723 741 743 749 753 759 760 760 779 791 792 802 822 832 854 856 858 864 873 875 880 890 891 897 915 933 935 953 971 977 978 988 993 995 1004 1010 1012 1014 1036 1046 1066 1076 1077 1089 1108 1108 1109 1115 1119 1125 1127 1145 1151 1167 1179 1190 1202 1207 1214 1220 1222 1224 1228 1258 1264 1264 1282 1294 1301 1305 1310 1323 1327 1327 1327 1333 1337 1363 1395 1395 1406 1414 1425 1426 1432 1438 1438 1440 1440 1442 1450 1466 1494 1519 1539 1539 1541 1543 1545 1551 1551 1553 1555 1581 1581 1581 1597 1618 1640 1650 1652 1654 1656 1656 1656 1666 1682 1712 1712 1731 1737 1737 1753 1753 1755 1755 1755 1765 1767 1769 1769 1868')
   def _11():
      linearSpectrum('TRPVGPAHTDCDSGNMEHDPGPPCSQGSHWQTGAQKRT') >> sepWithSp >> testVs(
         '0 57 57 57 57 57 71 71 87 87 87 97 97 97 97 97 99 101 101 101 101 103 103 114 115 115 115 128 128 128 128 128 129 131 137 137 137 144 144 154 154 154 156 156 156 158 168 171 185 186 190 194 196 199 200 202 208 212 215 216 218 218 224 225 229 229 238 245 251 251 252 253 253 253 256 256 257 257 258 259 260 266 269 272 272 281 284 286 287 297 302 305 305 309 314 318 319 323 324 327 333 348 349 350 352 353 354 354 357 357 359 362 362 366 373 374 375 381 384 384 385 389 397 406 406 409 409 410 412 415 415 420 421 424 431 434 441 451 451 453 456 461 462 463 463 467 472 472 476 477 478 483 485 485 496 503 504 506 510 511 512 512 513 518 521 521 527 535 538 538 540 543 552 558 559 562 566 568 569 569 571 577 578 578 584 591 595 595 599 600 607 607 609 609 613 624 626 626 632 633 639 641 641 642 653 655 656 658 659 666 666 671 677 678 680 681 682 683 692 696 696 696 703 713 714 715 722 723 723 723 729 729 736 739 742 753 763 767 769 770 770 774 780 780 780 781 785 786 790 793 796 799 808 810 810 815 815 823 824 824 826 829 832 837 838 850 851 860 867 870 873 877 877 881 882 883 883 885 895 895 900 911 913 916 918 919 924 925 930 934 936 940 947 952 952 952 955 960 963 968 974 975 979 982 982 988 988 992 997 1010 1014 1021 1023 1031 1031 1031 1033 1036 1039 1039 1039 1047 1050 1054 1056 1062 1062 1071 1077 1079 1080 1080 1085 1089 1089 1092 1103 1104 1107 1111 1118 1128 1133 1134 1134 1136 1136 1142 1142 1148 1153 1160 1164 1164 1167 1168 1178 1179 1185 1191 1193 1199 1200 1204 1208 1208 1221 1221 1226 1233 1235 1235 1236 1239 1239 1248 1249 1250 1257 1257 1261 1265 1265 1270 1280 1284 1292 1292 1295 1297 1301 1308 1314 1322 1322 1328 1336 1336 1336 1336 1337 1341 1349 1349 1354 1358 1362 1364 1367 1376 1381 1385 1393 1393 1394 1398 1406 1406 1412 1413 1419 1423 1436 1436 1438 1439 1451 1451 1451 1455 1459 1464 1465 1477 1490 1493 1493 1495 1495 1507 1509 1509 1510 1513 1514 1521 1526 1534 1537 1550 1551 1552 1552 1554 1554 1566 1566 1573 1580 1592 1592 1605 1606 1608 1614 1618 1630 1638 1641 1642 1645 1647 1649 1651 1654 1655 1655 1663 1663 1663 1665 1666 1671 1689 1695 1711 1717 1720 1733 1742 1742 1743 1746 1748 1752 1759 1760 1760 1762 1762 1767 1769 1773 1792 1798 1800 1803 1805 1816 1817 1819 1826 1832 1849 1857 1859 1861 1863 1870 1870 1871 1874 1879 1887 1902 1903 1904 1906 1913 1914 1916 1916 1918 1927 1931 1935 1944 1950 1960 1988 1998 1999 2002 2003 2007 2013 2013 2014 2015 2017 2017 2018 2019 2031 2045 2045 2047 2050 2064 2072 2078 2102 2104 2110 2116 2116 2116 2118 2121 2127 2130 2132 2135 2146 2151 2151 2154 2169 2173 2173 2175 2189 2203 2213 2222 2232 2232 2236 2244 2247 2249 2255 2258 2260 2266 2270 2283 2288 2289 2300 2301 2304 2319 2331 2337 2350 2359 2364 2367 2369 2372 2375 2376 2384 2388 2388 2407 2414 2428 2429 2456 2456 2465 2465 2470 2474 2475 2478 2485 2503 2513 2515 2516 2522 2528 2545 2557 2566 2572 2584 2585 2593 2602 2606 2612 2623 2629 2631 2641 2642 2672 2673 2685 2686 2694 2699 2703 2709 2721 2728 2734 2742 2760 2770 2773 2774 2787 2798 2822 2827 2829 2831 2831 2849 2865 2871 2888 2890 2895 2902 2926 2928 2928 2950 2959 2966 2985 2991 2999 3005 3023 3027 3030 3051 3056 3084 3087 3106 3106 3124 3127 3152 3155 3158 3179 3181 3184 3207 3243 3252 3255 3280 3280 3283 3312 3314 3337 3344 3380 3381 3408 3411 3411 3415 3438 3468 3508 3509 3512 3536 3567 3569 3637 3664 3664 3668 3765 3765 3820 3921 3921 4022')
      linearSpectrum('WQMMSPIGSCNAWFYPKLMSKWSSRDAVTFVAKMIQYMHMRPPT') >> sepWithSp >> testVs(
         '0 57 71 71 71 87 87 87 87 87 97 97 97 97 99 99 101 101 103 113 113 113 114 115 128 128 128 128 128 131 131 131 131 131 131 137 144 147 147 156 156 163 163 170 170 170 174 184 185 186 186 186 186 190 194 198 199 200 210 215 217 218 218 225 241 241 243 244 244 246 247 248 253 257 257 259 259 260 262 267 268 268 271 271 273 285 287 288 291 294 295 297 298 304 310 314 314 315 317 330 330 331 333 338 342 346 347 347 349 350 354 354 358 360 360 361 371 372 372 372 375 384 386 388 390 399 401 401 404 404 407 418 418 422 424 428 429 429 431 432 441 441 443 445 445 445 446 446 451 457 459 459 469 474 474 477 481 485 488 488 496 500 501 516 516 517 517 518 521 528 532 533 535 535 535 542 542 544 545 546 555 556 559 559 561 562 567 571 571 572 574 575 576 576 582 587 588 593 615 616 618 618 619 621 629 631 632 632 642 644 645 645 648 652 658 663 663 664 666 670 672 675 677 681 684 687 689 689 690 702 703 703 706 708 716 716 718 719 719 721 729 731 731 732 734 744 749 759 760 765 773 776 776 778 779 784 788 789 790 792 794 801 803 803 806 815 817 819 828 830 831 831 833 834 846 846 847 847 850 859 860 860 862 863 865 866 870 871 873 875 878 881 889 902 905 906 912 915 917 918 920 928 929 930 931 934 934 943 946 947 957 959 960 962 962 964 965 968 975 975 977 980 987 991 994 1002 1009 1013 1016 1017 1017 1019 1025 1030 1033 1033 1036 1040 1041 1044 1046 1048 1048 1049 1052 1056 1062 1062 1074 1075 1081 1088 1090 1090 1096 1101 1103 1111 1117 1118 1119 1120 1120 1122 1123 1133 1138 1138 1141 1147 1148 1150 1153 1153 1161 1161 1177 1177 1180 1180 1180 1187 1193 1200 1203 1205 1207 1209 1212 1218 1218 1219 1225 1231 1232 1234 1235 1237 1248 1248 1251 1251 1253 1254 1260 1264 1266 1266 1267 1276 1284 1289 1289 1292 1305 1305 1311 1315 1315 1318 1322 1324 1340 1340 1347 1347 1349 1354 1356 1359 1361 1363 1363 1363 1365 1366 1366 1379 1379 1382 1385 1386 1386 1388 1388 1395 1397 1405 1412 1427 1434 1437 1446 1448 1450 1452 1453 1453 1468 1475 1476 1478 1478 1480 1483 1484 1485 1485 1487 1489 1491 1494 1497 1508 1510 1510 1513 1519 1522 1524 1533 1535 1540 1549 1551 1555 1562 1563 1565 1579 1581 1582 1584 1584 1586 1591 1597 1606 1607 1607 1609 1611 1612 1615 1625 1632 1634 1636 1636 1638 1638 1648 1650 1653 1654 1678 1683 1693 1693 1694 1694 1694 1696 1696 1696 1712 1712 1719 1719 1725 1725 1729 1733 1733 1735 1735 1740 1741 1741 1749 1765 1767 1781 1790 1795 1798 1801 1806 1806 1806 1806 1811 1822 1824 1825 1825 1827 1828 1828 1830 1830 1832 1832 1840 1847 1877 1881 1882 1882 1882 1885 1896 1896 1898 1903 1903 1909 1911 1912 1915 1921 1921 1929 1931 1934 1934 1937 1937 1953 1953 1956 1964 1972 1981 1984 1995 1996 1998 2000 2008 2008 2010 2013 2018 2026 2030 2031 2040 2043 2043 2050 2052 2065 2065 2066 2067 2071 2077 2082 2084 2085 2095 2095 2095 2097 2099 2101 2115 2128 2139 2141 2142 2150 2153 2162 2164 2166 2170 2171 2171 2174 2178 2178 2182 2182 2186 2194 2213 2216 2226 2228 2228 2229 2241 2243 2251 2257 2261 2267 2269 2269 2270 2271 2275 2278 2281 2299 2300 2306 2313 2314 2325 2328 2338 2341 2341 2348 2356 2356 2357 2357 2358 2359 2365 2370 2372 2399 2399 2400 2403 2409 2413 2414 2425 2427 2437 2438 2444 2445 2453 2457 2459 2469 2470 2472 2472 2485 2485 2496 2496 2513 2514 2517 2524 2526 2527 2531 2534 2540 2546 2556 2565 2566 2566 2572 2584 2585 2598 2600 2604 2609 2611 2616 2623 2627 2627 2631 2652 2658 2659 2661 2662 2671 2671 2687 2687 2697 2703 2710 2712 2713 2724 2729 2729 2732 2737 2740 2742 2749 2758 2759 2760 2771 2774 2774 2783 2802 2811 2815 2815 2831 2834 2841 2842 2843 2845 2846 2860 2860 2868 2871 2873 2873 2876 2880 2896 2899 2902 2930 2942 2944 2946 2947 2956 2958 2959 2965 2970 2970 2972 2977 2993 2997 3001 3001 3007 3024 3033 3041 3057 3059 3062 3072 3073 3078 3084 3089 3090 3090 3100 3116 3121 3121 3128 3128 3133 3144 3146 3169 3187 3187 3188 3191 3193 3201 3203 3203 3218 3218 3220 3247 3256 3259 3264 3274 3275 3284 3286 3300 3315 3316 3319 3319 3330 3331 3348 3350 3378 3381 3387 3387 3387 3390 3401 3413 3416 3431 3437 3444 3447 3461 3478 3481 3494 3500 3515 3518 3518 3518 3528 3532 3534 3541 3568 3579 3607 3617 3618 3625 3625 3628 3631 3633 3646 3646 3649 3688 3704 3704 3705 3714 3726 3738 3749 3759 3762 3762 3777 3785 3791 3802 3811 3832 3835 3836 3875 3882 3890 3890 3893 3899 3905 3912 3922 3922 3963 3972 3983 3992 3996 4002 4006 4018 4049 4053 4053 4059 4076 4089 4097 4099 4103 4146 4162 4181 4184 4186 4190 4190 4200 4204 4243 4259 4259 4287 4312 4321 4321 4344 4346 4356 4356 4367 4443 4449 4452 4453 4457 4477 4498 4540 4554 4574 4580 4608 4635 4641 4671 4705 4736 4766 4772 4802 4833 4903 4922 4930 5019 5031 5116 5217')

   def _5():
      cnt, _ = countPeptidesFrom(1024)
      cnt >> testVs(14712706211)
      cnt, _ = countPeptidesFrom(1322)
      cnt >> testVs(51972523134741)
      cnt, _ = countPeptidesFrom(1337)
      cnt >> testVs(78648746786215)

   def _6():
      linearSubpeptidesCount(4) >> testVs(11)
      linearSubpeptidesCount(37704) >> testVs(710814661)
      linearSubpeptidesCount(19664) >> testVs(11)

      ints = parseInts('0 113 128 186 241 299 314 427')
      cyclopeptidesSequencing(ints) >> spaceJoin >> testVs('113-128-186 113-186-128 128-113-186 128-186-113 186-113-128 186-128-113')

      ints = parseInts('0 71 97 99 103 113 113 114 115 131 137 196 200 202 208 214 226 227 228 240 245 299 311 311 316 327 337 339 340 341 358 408 414 424 429 436 440 442 453 455 471 507 527 537 539 542 551 554 556 566 586 622 638 640 651 653 657 664 669 679 685 735 752 753 754 756 766 777 782 782 794 848 853 865 866 867 879 885 891 893 897 956 962 978 979 980 980 990 994 996 1022 1093')
      cyclopeptidesSequencing(ints) >> spaceJoin >> testVs('103-137-71-131-114-113-113-115-99-97 103-97-99-115-113-113-114-131-71-137 113-113-114-131-71-137-103-97-99-115 113-113-115-99-97-103-137-71-131-114 113-114-131-71-137-103-97-99-115-113 113-115-99-97-103-137-71-131-114-113 114-113-113-115-99-97-103-137-71-131 114-131-71-137-103-97-99-115-113-113 115-113-113-114-131-71-137-103-97-99 115-99-97-103-137-71-131-114-113-113 131-114-113-113-115-99-97-103-137-71 131-71-137-103-97-99-115-113-113-114 137-103-97-99-115-113-113-114-131-71 137-71-131-114-113-113-115-99-97-103 71-131-114-113-113-115-99-97-103-137 71-137-103-97-99-115-113-113-114-131 97-103-137-71-131-114-113-113-115-99 97-99-115-113-113-114-131-71-137-103 99-115-113-113-114-131-71-137-103-97 99-97-103-137-71-131-114-113-113-115')

      ints = parseInts('0 97 101 114 114 114 115 129 147 156 163 186 211 211 215 216 229 270 292 300 303 310 312 315 325 330 330 385 397 417 426 427 429 439 444 466 478 486 511 526 532 541 541 580 592 595 600 612 625 633 640 655 689 695 697 709 727 739 741 747 781 796 803 811 824 836 841 844 856 895 895 904 910 925 950 958 970 992 997 1007 1009 1010 1019 1039 1051 1106 1106 1111 1121 1124 1126 1133 1136 1144 1166 1207 1220 1221 1225 1225 1250 1273 1280 1289 1307 1321 1322 1322 1322 1335 1339 1436')
      cyclopeptidesSequencing(ints) >> spaceJoin >> testVs('101-114-97-114-186-129-163-147-156-114-115 101-115-114-156-147-163-129-186-114-97-114 114-101-115-114-156-147-163-129-186-114-97 114-115-101-114-97-114-186-129-163-147-156 114-156-147-163-129-186-114-97-114-101-115 114-186-129-163-147-156-114-115-101-114-97 114-97-114-101-115-114-156-147-163-129-186 114-97-114-186-129-163-147-156-114-115-101 115-101-114-97-114-186-129-163-147-156-114 115-114-156-147-163-129-186-114-97-114-101 129-163-147-156-114-115-101-114-97-114-186 129-186-114-97-114-101-115-114-156-147-163 147-156-114-115-101-114-97-114-186-129-163 147-163-129-186-114-97-114-101-115-114-156 156-114-115-101-114-97-114-186-129-163-147 156-147-163-129-186-114-97-114-101-115-114 163-129-186-114-97-114-101-115-114-156-147 163-147-156-114-115-101-114-97-114-186-129 186-114-97-114-101-115-114-156-147-163-129 186-129-163-147-156-114-115-101-114-97-114 97-114-101-115-114-156-147-163-129-186-114 97-114-186-129-163-147-156-114-115-101-114')

   def _7():
      #* cyclic/linear Peptide score of given amnio acid and spectrum
      ints = parseInts('0 99 113 114 128 227 257 299 355 356 370 371 484')
      peptideScore('NQEL', ints, 'cyc') >> testVs(11)

      l = readlines('dataset_29915_3a.txt')
      peptideScore(l[0], parseInts(l[1]), 'cyc') >> testVs(394)

      l = readlines('dataset_29915_3.txt')
      peptideScore(l[0], parseInts(l[1]), 'cyc') >> testVs(737)

      l = readlines('dataset_4913_1.txt')
      peptideScore(l[0], parseInts(l[1]), 'linear') >> testVs(296)

      # 1.1.8 - Find top peptide
      ints = parseInts('0 71 113 129 147 200 218 260 313 331 347 389 460')
      topCyclopeptideSequencing(ints, 10) >> toDashed >> testVs('113-147-71-129')

      ints = parseInts('0 71 71 87 97 99 99 101 103 103 113 114 115 115 128 131 147 156 156 156 163 163 174 184 185 186 186 196 206 215 230 232 243 246 246 257 257 266 269 270 271 276 277 278 283 287 301 312 312 319 328 330 331 341 343 345 360 369 371 377 381 384 402 406 411 418 426 427 430 431 432 432 440 442 442 457 463 475 476 478 484 493 497 509 517 526 527 529 533 534 545 547 555 558 558 570 573 577 588 588 589 590 596 598 612 626 632 648 657 660 660 664 671 673 674 683 683 689 689 693 697 699 701 702 703 733 741 744 754 761 763 770 772 773 779 788 788 792 796 802 804 811 812 834 836 846 857 858 859 860 867 869 878 889 889 895 901 903 917 926 929 935 935 939 940 947 956 959 960 965 966 967 973 975 990 1004 1020 1027 1030 1034 1038 1042 1044 1045 1048 1053 1054 1062 1064 1066 1075 1082 1113 1115 1121 1123 1123 1124 1131 1133 1141 1141 1146 1147 1152 1159 1167 1169 1186 1190 1190 1210 1211 1212 1218 1220 1223 1230 1230 1236 1238 1238 1244 1260 1270 1272 1277 1286 1289 1297 1299 1301 1305 1310 1331 1331 1337 1343 1346 1353 1366 1367 1370 1373 1376 1391 1392 1394 1398 1399 1401 1404 1407 1416 1444 1453 1456 1459 1461 1462 1466 1468 1469 1484 1487 1490 1493 1494 1507 1514 1517 1523 1529 1529 1550 1555 1559 1561 1563 1571 1574 1583 1588 1590 1600 1616 1622 1622 1624 1630 1630 1637 1640 1642 1648 1649 1650 1670 1670 1674 1691 1693 1701 1708 1713 1714 1719 1719 1727 1729 1736 1737 1737 1739 1745 1747 1778 1785 1794 1796 1798 1806 1807 1812 1815 1816 1818 1822 1826 1830 1833 1840 1856 1870 1885 1887 1893 1894 1895 1900 1901 1904 1913 1920 1921 1925 1925 1931 1934 1957 1959 1965 1971 1971 1982 1991 1993 2000 2001 2002 2003 2014 2024 2026 2048 2049 2056 2058 2064 2068 2072 2072 2081 2087 2088 2090 2097 2099 2106 2116 2119 2127 2157 2158 2159 2161 2163 2167 2171 2171 2177 2177 2186 2187 2189 2196 2200 2200 2203 2212 2228 2234 2248 2262 2264 2270 2271 2272 2272 2283 2287 2290 2302 2302 2305 2313 2315 2326 2327 2331 2333 2334 2343 2351 2363 2367 2376 2382 2384 2385 2397 2403 2418 2418 2420 2428 2428 2429 2430 2433 2434 2442 2449 2454 2458 2476 2479 2483 2489 2491 2500 2515 2517 2519 2529 2530 2532 2541 2548 2548 2559 2573 2577 2582 2583 2584 2589 2590 2591 2594 2603 2603 2614 2614 2617 2628 2630 2645 2654 2664 2674 2674 2675 2676 2686 2697 2697 2704 2704 2704 2713 2729 2732 2745 2745 2746 2747 2757 2757 2759 2761 2761 2763 2773 2789 2789 2860')
      res = topCyclopeptideSequencing(ints, 152) 
      # print(res >> toDashed)
      testNoBegEnd(res >> toStrs, '99-147-99-97-87-128-115-163-103-103-71-186-71-114-156-156-163-113-156-115-186-101-131'.split('-'))

      l = readlines('dataset_29915_8.txt')
      [n, intsStr] = l
      res = topCyclopeptideSequencing(parseInts(intsStr), parseInts(n)[0])
      # print(res >> toDashed)
      testNoBegEnd(res >> toStrs, '137-137-97-87-147-114-128-99-131-128-147-147-156-101-128-163-101-113-129-97-115-128-71'.split('-'))

      l = readlines('dataset_29915_10.txt')
      [n, intsStr] = l
      res = topPeptidesSequencing(parseInts(intsStr), parseInts(n)[0])
      strs = [toDashed(r) for r in res]
      print(' '.join(strs))

   def _9():
      #* Spectral Convolution 
      ints = parseInts('0 137 186 323')
      res = spectralConvolution(ints, 1, 999)
      res >> testVs([137, 186, 49, 323, 186, 137])

      ints = parseInts('0 57 71 87 103 115 115 158 163 163 172 186 190 220 261 266 273 278 278 335 335 335 349 353 376 381 406 424 436 438 450 468 493 498 521 525 539 539 539 596 596 601 608 613 654 684 688 702 711 711 716 759 759 771 787 803 817 874')
      res = spectralConvolution(ints, 1, 999)
      res >> sepWithSp >> testVs('57 71 14 87 30 16 103 46 32 16 115 58 44 28 12 115 58 44 28 12 158 101 87 71 55 43 43 163 106 92 76 60 48 48 5 163 106 92 76 60 48 48 5 172 115 101 85 69 57 57 14 9 9 186 129 115 99 83 71 71 28 23 23 14 190 133 119 103 87 75 75 32 27 27 18 4 220 163 149 133 117 105 105 62 57 57 48 34 30 261 204 190 174 158 146 146 103 98 98 89 75 71 41 266 209 195 179 163 151 151 108 103 103 94 80 76 46 5 273 216 202 186 170 158 158 115 110 110 101 87 83 53 12 7 278 221 207 191 175 163 163 120 115 115 106 92 88 58 17 12 5 278 221 207 191 175 163 163 120 115 115 106 92 88 58 17 12 5 335 278 264 248 232 220 220 177 172 172 163 149 145 115 74 69 62 57 57 335 278 264 248 232 220 220 177 172 172 163 149 145 115 74 69 62 57 57 335 278 264 248 232 220 220 177 172 172 163 149 145 115 74 69 62 57 57 349 292 278 262 246 234 234 191 186 186 177 163 159 129 88 83 76 71 71 14 14 14 353 296 282 266 250 238 238 195 190 190 181 167 163 133 92 87 80 75 75 18 18 18 4 376 319 305 289 273 261 261 218 213 213 204 190 186 156 115 110 103 98 98 41 41 41 27 23 381 324 310 294 278 266 266 223 218 218 209 195 191 161 120 115 108 103 103 46 46 46 32 28 5 406 349 335 319 303 291 291 248 243 243 234 220 216 186 145 140 133 128 128 71 71 71 57 53 30 25 424 367 353 337 321 309 309 266 261 261 252 238 234 204 163 158 151 146 146 89 89 89 75 71 48 43 18 436 379 365 349 333 321 321 278 273 273 264 250 246 216 175 170 163 158 158 101 101 101 87 83 60 55 30 12 438 381 367 351 335 323 323 280 275 275 266 252 248 218 177 172 165 160 160 103 103 103 89 85 62 57 32 14 2 450 393 379 363 347 335 335 292 287 287 278 264 260 230 189 184 177 172 172 115 115 115 101 97 74 69 44 26 14 12 468 411 397 381 365 353 353 310 305 305 296 282 278 248 207 202 195 190 190 133 133 133 119 115 92 87 62 44 32 30 18 493 436 422 406 390 378 378 335 330 330 321 307 303 273 232 227 220 215 215 158 158 158 144 140 117 112 87 69 57 55 43 25 498 441 427 411 395 383 383 340 335 335 326 312 308 278 237 232 225 220 220 163 163 163 149 145 122 117 92 74 62 60 48 30 5 521 464 450 434 418 406 406 363 358 358 349 335 331 301 260 255 248 243 243 186 186 186 172 168 145 140 115 97 85 83 71 53 28 23 525 468 454 438 422 410 410 367 362 362 353 339 335 305 264 259 252 247 247 190 190 190 176 172 149 144 119 101 89 87 75 57 32 27 4 539 482 468 452 436 424 424 381 376 376 367 353 349 319 278 273 266 261 261 204 204 204 190 186 163 158 133 115 103 101 89 71 46 41 18 14 539 482 468 452 436 424 424 381 376 376 367 353 349 319 278 273 266 261 261 204 204 204 190 186 163 158 133 115 103 101 89 71 46 41 18 14 539 482 468 452 436 424 424 381 376 376 367 353 349 319 278 273 266 261 261 204 204 204 190 186 163 158 133 115 103 101 89 71 46 41 18 14 596 539 525 509 493 481 481 438 433 433 424 410 406 376 335 330 323 318 318 261 261 261 247 243 220 215 190 172 160 158 146 128 103 98 75 71 57 57 57 596 539 525 509 493 481 481 438 433 433 424 410 406 376 335 330 323 318 318 261 261 261 247 243 220 215 190 172 160 158 146 128 103 98 75 71 57 57 57 601 544 530 514 498 486 486 443 438 438 429 415 411 381 340 335 328 323 323 266 266 266 252 248 225 220 195 177 165 163 151 133 108 103 80 76 62 62 62 5 5 608 551 537 521 505 493 493 450 445 445 436 422 418 388 347 342 335 330 330 273 273 273 259 255 232 227 202 184 172 170 158 140 115 110 87 83 69 69 69 12 12 7 613 556 542 526 510 498 498 455 450 450 441 427 423 393 352 347 340 335 335 278 278 278 264 260 237 232 207 189 177 175 163 145 120 115 92 88 74 74 74 17 17 12 5 654 597 583 567 551 539 539 496 491 491 482 468 464 434 393 388 381 376 376 319 319 319 305 301 278 273 248 230 218 216 204 186 161 156 133 129 115 115 115 58 58 53 46 41 684 627 613 597 581 569 569 526 521 521 512 498 494 464 423 418 411 406 406 349 349 349 335 331 308 303 278 260 248 246 234 216 191 186 163 159 145 145 145 88 88 83 76 71 30 688 631 617 601 585 573 573 530 525 525 516 502 498 468 427 422 415 410 410 353 353 353 339 335 312 307 282 264 252 250 238 220 195 190 167 163 149 149 149 92 92 87 80 75 34 4 702 645 631 615 599 587 587 544 539 539 530 516 512 482 441 436 429 424 424 367 367 367 353 349 326 321 296 278 266 264 252 234 209 204 181 177 163 163 163 106 106 101 94 89 48 18 14 711 654 640 624 608 596 596 553 548 548 539 525 521 491 450 445 438 433 433 376 376 376 362 358 335 330 305 287 275 273 261 243 218 213 190 186 172 172 172 115 115 110 103 98 57 27 23 9 711 654 640 624 608 596 596 553 548 548 539 525 521 491 450 445 438 433 433 376 376 376 362 358 335 330 305 287 275 273 261 243 218 213 190 186 172 172 172 115 115 110 103 98 57 27 23 9 716 659 645 629 613 601 601 558 553 553 544 530 526 496 455 450 443 438 438 381 381 381 367 363 340 335 310 292 280 278 266 248 223 218 195 191 177 177 177 120 120 115 108 103 62 32 28 14 5 5 759 702 688 672 656 644 644 601 596 596 587 573 569 539 498 493 486 481 481 424 424 424 410 406 383 378 353 335 323 321 309 291 266 261 238 234 220 220 220 163 163 158 151 146 105 75 71 57 48 48 43 759 702 688 672 656 644 644 601 596 596 587 573 569 539 498 493 486 481 481 424 424 424 410 406 383 378 353 335 323 321 309 291 266 261 238 234 220 220 220 163 163 158 151 146 105 75 71 57 48 48 43 771 714 700 684 668 656 656 613 608 608 599 585 581 551 510 505 498 493 493 436 436 436 422 418 395 390 365 347 335 333 321 303 278 273 250 246 232 232 232 175 175 170 163 158 117 87 83 69 60 60 55 12 12 787 730 716 700 684 672 672 629 624 624 615 601 597 567 526 521 514 509 509 452 452 452 438 434 411 406 381 363 351 349 337 319 294 289 266 262 248 248 248 191 191 186 179 174 133 103 99 85 76 76 71 28 28 16 803 746 732 716 700 688 688 645 640 640 631 617 613 583 542 537 530 525 525 468 468 468 454 450 427 422 397 379 367 365 353 335 310 305 282 278 264 264 264 207 207 202 195 190 149 119 115 101 92 92 87 44 44 32 16 817 760 746 730 714 702 702 659 654 654 645 631 627 597 556 551 544 539 539 482 482 482 468 464 441 436 411 393 381 379 367 349 324 319 296 292 278 278 278 221 221 216 209 204 163 133 129 115 106 106 101 58 58 46 30 14 874 817 803 787 771 759 759 716 711 711 702 688 684 654 613 608 601 596 596 539 539 539 525 521 498 493 468 450 438 436 424 406 381 376 353 349 335 335 335 278 278 273 266 261 220 190 186 172 163 163 158 115 115 103 87 71 57')
      len(res) >> testVs(1641)

      # ints = parseInts('0 87 87 97 131 131 137 147 156 156 184 218 234 234 243 268 287 287 303 321 340 355 365 365 374 390 434 443 452 452 471 477 502 521 521 530 539 590 599 608 608 627 652 658 677 677 686 695 739 755 764 764 774 789 808 826 842 842 861 886 895 895 911 945 973 973 982 992 998 998 1032 1042 1042 1129')
      # res = spectralConvolution(ints, 1, 999)
      # res >> sepWithSp >> print # wrong

      # Convolution Cyclopedptide Sequencing
      ints = parseInts('57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493')
      res = convolutionSequencing(ints, 20, 60)
      strs = [toDashed(r) for r in res]
      print(' '.join(strs))
      # expects '99-71-137-57-72-57'

      ints = parseInts('0 57 99 99 101 101 103 103 114 114 114 128 156 158 160 200 200 202 213 215 228 231 242 257 259 261 270 301 314 316 316 327 330 343 345 358 360 373 384 415 417 430 430 442 444 444 446 461 472 483 487 516 531 543 544 545 545 558 558 575 586 586 600 617 630 643 644 645 646 659 672 689 703 703 714 731 731 744 744 745 746 758 773 802 806 817 828 843 845 845 847 859 859 872 874 905 916 929 931 944 946 959 962 973 973 975 988 1019 1028 1030 1032 1047 1058 1061 1074 1076 1087 1089 1089 1129 1131 1133 1161 1175 1175 1175 1186 1186 1188 1188 1190 1190 1232 1289')
      res = convolutionSequencing(ints, 19, 362)
      strs = [toDashed(r) for r in res]
      print(' '.join(strs))

   def _10():
      #* Convolution Cyclopedptide Sequencing on Spectrum25 
      nums = parseFloats('371.5 375.4 390.4 392.2 409.0 420.2 427.2 443.3 446.4 461.3 471.4 477.4 491.3 505.3 506.4 519.2 536.1 546.5 553.3 562.3 588.2 600.3 616.2 617.4 618.3 633.4 634.4 636.2 651.5 652.4 702.5 703.4 712.5 718.3 721.0 730.3 749.4 762.6 763.4 764.4 779.6 780.4 781.4 782.4 797.3 862.4 876.4 877.4 878.6 879.4 893.4 894.4 895.4 896.5 927.4 944.4 975.5 976.5 977.4 979.4 1005.5 1007.5 1022.5 1023.7 1024.5 1039.5 1040.3 1042.5 1043.4 1057.5 1119.6 1120.6 1137.6 1138.6 1139.5 1156.5 1157.6 1168.6 1171.6 1185.4 1220.6 1222.5 1223.6 1239.6 1240.6 1250.5 1256.5 1266.5 1267.5 1268.6')
      res = protonMassSpecCyclopeptideSequencing(nums, 20, 5000)
      res >> sepWithSp >> print
      # Couldn't get, all wrong:
      # N=10k: 164 98 128 145 116 129 163 82 113 130
      # N=20k: 163 82 163 129 113 130 146 115 97 130
      # N=40k: 114 147 115 129 130 130 129 81 82 98 113 (11 numbers, should be 10)

      # N=1000, 38 linear peptides of max score (83) can be found for spectrum 25
      l = readlines('Tyrocidine_B1_Spectrum_25.txt')
      ints = parseInts(l[0])
      res = topPeptidesSequencing(ints)
      res >> spaceDashed >> print # 113-128-99-163-128-114-147-115-71-97-147 113-128-99-163-128-114-147-71-115-97-147 114-128-163-99-128-113-147-97-71-115-147 115-147-114-128-163-99-128-113-147-97-71 115-71-147-114-128-163-99-128-113-147-97 128-114-147-115-71-97-147-113-128-99-163 128-99-163-128-114-147-115-71-97-147-113 128-99-163-128-114-147-71-115-97-147-113 147-113-128-99-163-128-114-147-115-71-97 147-113-128-99-163-128-114-147-71-115-97 147-114-128-163-99-128-113-147-97-115-71 147-114-128-163-99-128-113-147-97-71-115 163-128-114-147-115-71-97-147-113-128-99 71-115-147-114-128-163-99-128-113-147-97 97-147-113-128-99-163-128-114-147-115-71 97-147-113-128-99-163-128-114-147-71-115 97-147-113-128-99-163-128-57-57-147-115-71 97-71-115-147-114-128-163-99-128-113-147 99-163-128-114-147-115-71-97-147-113-128 99-163-128-114-147-71-115-97-147-113-128
      # res >> spaceDashed >> testVs('113-128-99-163-128-114-147-115-71-97-147 113-128-99-163-128-114-147-71-115-97-147 114-128-163-99-128-113-147-97-115-71-147 114-128-163-99-128-113-147-97-71-115-147 115-147-114-128-163-99-128-113-147-97-71 115-71-147-114-128-163-99-128-113-147-97 128-114-147-115-71-97-147-113-128-99-163 128-99-163-128-114-147-115-71-97-147-113 128-99-163-128-114-147-71-115-97-147-113 147-113-128-99-163-128-114-147-115-71-97 147-113-128-99-163-128-114-147-71-115-97 147-113-128-99-163-128-57-57-147-115-71-97 147-114-128-163-99-128-113-147-97-115-71 147-114-128-163-99-128-113-147-97-71-115 147-114-128-163-99-57-71-113-147-97-115-71 147-114-128-163-99-57-71-113-147-97-71-115 147-114-128-163-99-71-57-113-147-97-115-71 147-114-128-163-99-71-57-113-147-97-71-115 163-128-114-147-115-71-97-147-113-128-99 57-57-128-163-99-128-113-147-97-71-115-147 71-115-147-114-128-163-99-128-113-147-97 97-147-113-128-99-163-128-114-147-115-71 97-147-113-128-99-163-128-114-147-71-115 97-147-113-128-99-163-128-57-57-147-115-71 97-147-113-128-99-163-128-57-57-147-71-115 97-147-113-128-99-163-57-71-114-147-115-71 97-147-113-128-99-163-57-71-114-147-71-115 97-147-113-128-99-163-71-57-114-147-115-71 97-147-113-128-99-163-71-57-114-147-71-115 97-147-113-128-99-163-71-57-57-57-147-115-71 97-147-113-128-99-163-71-57-57-57-147-71-115 97-71-115-147-114-128-163-99-128-113-147 99-163-128-114-147-115-71-97-147-113-128 99-163-128-114-147-115-71-97-147-113-57-71 99-163-128-114-147-115-71-97-147-113-71-57 99-163-128-114-147-71-115-97-147-113-128 99-163-128-114-147-71-115-97-147-113-57-71 99-163-128-114-147-71-115-97-147-113-71-57')

      # N=1000, 38 linear peptides of max score (87) can be found for spectrum 10
      #  with additional non-proteinogenic amino acids, expanding the number of possible
      #  building blocks for antibiotic peptides from 20 to 144.
      # takes ~2mins 
      # l = readlines('Tyrocidine_B1_Spectrum_10.txt')
      # ints = parseInts(l[0])
      # res = topPeptidesSequencing(ints, 'extended')
      # res >> spaceDashed >> print # 113-147-97-186-147-114-128-98-65-99-128 114-147-186-97-147-113-128-99-65-98-128 128-113-147-97-186-147-114-128-98-65-99 128-114-147-186-97-147-113-128-99-65-98 128-99-163-128-114-147-113-73-97-147-113 147-114-128-163-99-128-113-147-97-73-113 147-186-97-147-113-128-99-65-98-128-114 147-97-186-147-114-128-98-65-99-128-113 97-186-147-114-128-98-65-99-128-113-147 99-163-128-114-147-113-73-97-147-113-128 99-163-128-114-147-124-62-97-147-113-128 99-163-128-114-147-186-97-72-75-113-128
      ## res >> spaceDashed >> testVs('113-147-97-186-147-114-128-98-65-99-128 113-147-97-186-147-114-128-98-65-99-57-71 113-147-97-186-147-114-128-98-65-99-58-70 113-147-97-186-147-114-128-98-65-99-59-69 113-147-97-186-147-114-128-98-65-99-60-68 113-147-97-186-147-114-128-98-65-99-61-67 113-147-97-186-147-114-128-98-65-99-62-66 113-147-97-186-147-114-128-98-65-99-63-65 113-147-97-186-147-114-128-98-65-99-64-64 113-147-97-186-147-114-128-98-65-99-65-63 113-147-97-186-147-114-128-98-65-99-66-62 113-147-97-186-147-114-128-98-65-99-67-61 113-147-97-186-147-114-128-98-65-99-68-60 113-147-97-186-147-114-128-98-65-99-69-59 113-147-97-186-147-114-128-98-65-99-70-58 113-147-97-186-147-114-128-98-65-99-71-57 113-147-97-186-147-114-71-57-163-99-66-62 114-147-186-97-147-113-128-99-65-98-128 128-113-147-97-186-147-114-128-98-65-99 128-114-147-186-97-147-113-128-99-65-98 128-99-163-128-114-147-113-73-97-147-113 147-114-128-163-99-128-113-147-97-73-113 147-186-97-147-113-128-99-65-98-128-114 147-186-97-147-113-128-99-65-98-128-57-57 147-97-186-147-114-128-98-65-99-128-113 97-147-113-128-99-163-128-114-73-74-58-66-62 97-186-147-114-128-98-65-99-128-113-147 97-186-74-73-114-128-163-99-128-113-75-72 99-163-128-114-147-113-73-97-147-113-128 99-163-128-114-147-113-73-97-147-113-57-71 99-163-128-114-147-124-62-97-147-113-128 99-163-128-114-147-124-62-97-147-113-57-71 99-163-128-114-147-186-97-72-75-113-128 99-163-128-114-147-98-88-97-147-113-57-71')

   def _13():
      ints = parseInts('0 99 113 114 128 227 257 299 355 356 370 371 484')
      peptideScore('NQEL', ints, 'linear') >> testVs(8)

      l = readlines('dataset_29921_1.txt')
      ints = parseInts(l[1])
      peptideScore(l[0], ints, 'linear') >> testVs(268)

      # Trim peptide candidates
      ints = parseInts('0 71 87 101 113 158 184 188 259 271 372')
      peptides = 'LAST ALST TLLT TQAS'.split(' ')
      trimPeptides(peptides, ints, 2) >> sepWithSp >> testVs('LAST ALST')

      l = readlines('dataset_29921_3a.txt')
      peptides = parseToStrs(l[0])
      ints, n = parseIntsList(l[1:])
      trimPeptides(peptides, ints, n[0]) >> sepWithSp >> testIn([
         'CLLEDPRMSCNSVSMGITTCTWNKIISGRTKNWHRA IGHCHLNEEHCEAPSGGMEPEYMKQQINWSFCKGLH CKQGICRTRAAASIRFASHREMKIYFHWMKECFCVR TKSCTLRSPTSFQSNGKTATFNPNKNMLAFFSPRHW SYMAGDASFANTQQWNYYWKYKVICPFCGFAFSIVS FDNMECLPRLNGKRCLFKDTNCVLKEHFQAKSVRWF', 
         'CLLEDPRMSCNSVSMGITTCTWNKIISGRTKNWHRA IGHCHLNEEHCEAPSGGMEPEYMKQQINWSFCKGLH CKQGICRTRAAASIRFASHREMKIYFHWMKECFCVR TKSCTLRSPTSFQSNGKTATFNPNKNMLAFFSPRHW FDNMECLPRLNGKRCLFKDTNCVLKEHFQAKSVRWF SYMAGDASFANTQQWNYYWKYKVICPFCGFAFSIVS',
      ]) 

      l = readlines('dataset_29921_3.txt')
      peptides = parseToStrs(l[0])
      ints, n = parseIntsList(l[1:])
      trimPeptides(peptides, ints, n[0]) >> sepWithSp >> testVs(
         'HTQLKRYVFVFGINGGGIEAKATAPAIKFFAEHGHDEPTFNIKCESLWE YFREPKTGYMHDNGHGVVFTYQNEGDDFHAHRDLDIIKQTMTYHPTCPH KFYECRWCMECGGLPMGTNGVDVWDATWSVQCVCNWMVRQWRFRTWYLG RRFENRVIWNWVMWQEVSMSCYPMHEQGTGTWEWQNYYPYAMTEKETSC HQDIDSKDRKWTVLDPASQYTCNEAGIEATSYKGTDHISEASETVEQYW RFNSDEDPTLDLWKHVQYDNEMHFEVKLSSLYFNRRHSMKWGENIAYFC')

if run(1):
   chap4._1()
if run(2):
   chap4._2()
if run(3):
   chap4._3()
if run(4):
   chap4._4()
if run(5):
   chap4._5()
if run(6):
   chap4._6()
if run(7):
   chap4._7()
if run(9):
   chap4._9()
if run(10):
   chap4._10()
if run(11):
   chap4._11()
if run(13):
   chap4._13()