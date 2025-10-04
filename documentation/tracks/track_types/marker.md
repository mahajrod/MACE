**Allowed subtypes**:
1. **middle** or **m** - markers localize at the middle of the region (*default*).
2. **ends** or **e** - markers localize at both ends of the region.
3. **right** or **r** - markers localize at the right end of the region.
4. **left** or **l** - markers localize at the left end of the region.


**Allowed types of the markers**:

1. **rectangle** or **r**     - draw a rectangle marker with start and stop on the edges of region, but with height smaller than track.
2. **ellipse** or **e**   - draw an ellipse with center at the middle between start and stop, height and x/y ratio are controlled by the feature style. Can be used to draw a circle if X and Y axis are in a different scale 
3. **circle** or **c**    - draw a circle with center at the middle between start and stop. It will look as an ellipse on the plot, if  X and Y axis are in a different scale. Use **ellipse** type to avoid this issue.

**Example**:

```
#scaffold_id	start	end	AAAA&marker$left	AAAA&colors	BBBB&marker$right	BBBB&colors
HiC_scaffold_1	0	10000	r	red e   blue
HiC_scaffold_1	100000	1100000	r   red e   blue
HiC_scaffold_1	2000000 2200000	e   red e   blue
```

Example contains two *marker* (AAAA&**marker**$left and BBBB&**marker**$right) tracks: **AAAA** and **BBBB**. Markers of track **AAAA** will be located at the left end of 
the region (AAAA&marker$**left** ), of track **BBBB** - on the right end (BBBB&marker$**right**), respectively. 
Two first regions of track **AAAA** will be marked by rectangles (*r*), and third - by ellipse (*e*). All markers of this track will be red.
For track **BBBB** you will see on the plot blue (*blue*) ellipses (*e*) for all the regions.