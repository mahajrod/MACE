**Allowed subtypes**:
1. position subtypes:
   1. **middle** or **m** - points localize at the middle of the region (*default*).
   2. **ends** or **e** - points localize at both ends of the region.
   3. **right** or **r** - points localize at the right end of the region.
   4. **left** or **l** - points localize at the left end of the region.
2. color subtypes - you can set a color (recognizable by matplotlib) of the curve 


**Example**:

```
#scaffold_id	start	end	AAAA&plot$left$red  BBBB&plot$blue  CCCC&plot
HiC_scaffold_1	0	10000   5   11   2.1
HiC_scaffold_1	100000	1100000   6   5 2.8
HiC_scaffold_1	2000000 2200000   8   4 3.7
```

Example contains three *plot* (AAAA&**plot**$left$red, BBBB&**plot**$blue and CCCC&**plot**) tracks: **AAAA**, **BBBB** and **CCCC**. 
Points of the curve for track **AAAA** will be located at the left end of 
the region (AAAA&plot$**left**$red ), of track **BBBB** and **CCCC** - in the middle as a default position will be used 
(as there is no specification of the position in a column name BBBB&plot$blue and CCCC&plot), respectively. 
Track **AAAA** will have a red curve (AAAA&plot$left$**red**), track **BBBB** - blue (BBBB&plot$**blue**), track **CCCC** - a curve with a default color (no specification of the position in a column name).