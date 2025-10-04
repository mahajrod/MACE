
# Main drawing objects (nested, from top to bottom)

1. Figure (Maybe remove)
2. Subplot
3. VerticalTrackGroup
4. HorizontalTrackGroup
5. Track

# Figure scheme

![image](figure_scheme.png)

1. Dimensions of the subplot:
   1. Horizontal dimensions
      1. S<sub>sl</sub>   - left space 
      2. S<sub>vlw</sub>  - width of the vertical track group label 
      3. S<sub>vls</sub>  - space between vertical track group label and connector between labels
      4. S<sub>cw</sub>   - connector width
      5. S<sub>cs</sub>   - space between connector and horizontal track group label
      6. S<sub>hlw</sub>  - width of horizontal track group label 
      7. S<sub>hls</sub>  - space between horizontal track group label and subplot's Y axis
      8. S<sub>ylw</sub>  - width of Y axis (including ticks and label)
      9. S<sub>ys</sub>   - space between Y axis and vertical track groups
      10. S<sub>bw</sub>  - width of subplot body
      11. S<sub>lsl</sub> - space between vertical track groups and legend
      12. S<sub>lw</sub>  - width of the legend
      13. S<sub>sr</sub>  - right space
   2. Vertical dimensions
      1. V<sub>h1</sub>  - height of the first vertical track group
      2. S<sub>vs</sub>   - space between vertical track groups
      3. V<sub>h2</sub>  - height of the second vertical track group
      
          ..................................................................
   
      4. S<sub>vhn</sub>  - height of the last vertical track group
      5. S<sub>ybs</sub>  - distance between left bottom corner of the last vertical track group and start of the subplot's Y axis
      6. S<sub>xst</sub>  - space between the last vertical track group and subplot's X axis 
      7. S<sub>xh</sub>   - height of the subplot's X axis (including ticks and label)
      8. S<sub>sb</sub>   - bottom space
2. Dimensions of the vertical track group
   1. Horizontal dimensions
      1. V<sub>sl</sub>   - left space 
      2. V<sub>yw</sub>   - width of the vertical track group's Y axis (including ticks and label)
      3. V<sub>ys</sub>   - space between vertical track group's Y axis and horizontal track groups
      4. V<sub>bw</sub>   - width width of the vertical track group body 
      5. V<sub>sr</sub>   - right space
   2. Vertical dimensions
      1. V<sub>st</sub>   - top space
      2. V<sub>bh</sub>   - height of the vertical track group body
      3. H<sub>h1</sub>   - height of the first horizontal track group
      4. V<sub>his</sub>  - space between horizontal track groups
      5. H<sub>h2</sub>   - height of the second horizontal track group
    
         .....................................................................

      6. H<sub>hn</sub>   - height of the last horizontal track group
      7. V<sub>xst</sub>  - space between the last horizontal track group and vertical track group's X axis 
      8. V<sub>xh</sub>   - height of the X axis (including ticks and label)
      9. V<sub>sb</sub>   - bottom space
3. Dimensions of the horizontal track group
   1. Horizontal dimensions
      1. H<sub>sl</sub>   - left space
      2. H<sub>yw</sub>   - width of the horizontal track group's Y axis (including ticks and label)
      3. H<sub>ys</sub>   - space between horizontal track group's Y axis and first track
      4. H<sub>bw</sub>   - width of the horizontal track group body 
      5. T<sub>w1</sub>   - width of the first track
      6. H<sub>sti</sub>  - space between tracks
      7. T<sub>w2</sub>   - width of the second track
    
         .....................................................................

      8. T<sub>wn</sub>  - width of the last track
      9. H<sub>sb</sub>   - right space
   2. Vertical dimensions
      1. H<sub>st</sub>   - top space
      2. H<sub>ht</sub>   - height of the tracks
      3. H<sub>xst</sub>  - space between tracks and horizontal track group's X axis 
      4. H<sub>xh</sub>   - height of the horizontal track group's X axis (including ticks and label)
      5. H<sub>sb</sub>   - bottom space
4. Dimensions of the track
   1. Horizontal dimensions
      1. T<sub>sl</sub>   - left space
      2. T<sub>yw</sub>   - width of the track's Y axis (including ticks and label)
      3. T<sub>ys</sub>   - space between track's Y axis and track body
      4. T<sub>bw</sub>    - width of the track body
      5. T<sub>sr</sub>    - right space
   2. Vertical dimensions
      1. T<sub>st</sub>   - top space
      2. T<sub>lh</sub>   - height of the track label
      3. T<sub>lsb</sub>  - space between track label and track body
      4. T<sub>bh</sub>    - height of the  track body
      5. T<sub>xst</sub>  - space between track bodies and track's X axis 
      6. T<sub>xh</sub>   - height of the track's X axis (including ticks and label)
      7. T<sub>sb</sub>   - bottom space

# Formulas:

## Track width
T<sub>w</sub> = T<sub>sl</sub> + T<sub>yw</sub> + T<sub>ys</sub> + T<sub>bw</sub> + T<sub>sr</sub>

| dimension | T<sub>sl</sub>       | T<sub>yw</sub>  | T<sub>ys</sub>   | T<sub>bw</sub> | T<sub>sr</sub> |
|-----------|----------------------|-----------------|------------------|----------------|----------------|
| source    | style                | auto from style | style            | data           | style          |
| values    | fraction<sup>*</sup> | fraction        | fraction         | absolute       | fraction       |

\* fraction of max (T<sub>bw</sub>)

## Track height
T<sub>h</sub> = T<sub>st</sub> + T<sub>lh</sub> + T<sub>lsb</sub> + T<sub>bh</sub> + T<sub>xst</sub>  + T<sub>xh</sub> + T<sub>sb</sub>

| dimension | T<sub>st</sub> | T<sub>lh</sub>   | T<sub>lsb</sub> | T<sub>bh</sub> | T<sub>xst</sub> | T<sub>xh</sub>  | T<sub>sb</sub> |
|-----------|----------------|------------------|-----------------|----------------|-----------------|-----------------|----------------|
| source    | style          | auto from style  | style           | data           | style           | auto from style | style          |
| values    | absolute       | absolute         | absolute        | absolute       | absolute        | absolute        | absolute       |


## Horizontal track group width
H<sub>w</sub> = H<sub>sl</sub> + H<sub>yw</sub> + H<sub>ys</sub> + H<sub>bw</sub> + H<sub>sr</sub>
                                 
H<sub>bw</sub> = sum(T<sub>wi</sub>) + (N-1)* H<sub>sti</sub>

| dimension | H<sub>sl</sub>       | H<sub>yw</sub>  | H<sub>ys</sub> | T<sub>wi</sub> | H<sub>sti</sub> | H<sub>sr</sub> |
|-----------|----------------------|-----------------|----------------|----------------|-----------------|----------------|
| source    | style                | auto from style | style          | data           | style           | style          |
| values    | fraction<sup>*</sup> | fraction        | fraction       | absolute       | fraction        | fraction       |

\* fraction max(sum(T<sub>bw</sub>) per horizontal track group) 


## Horizontal track group height
H<sub>h</sub> = H<sub>st</sub> + H<sub>ht</sub> + H<sub>xst</sub> + H<sub>xh</sub> + H<sub>sb</sub>

| dimension | H<sub>st</sub> | H<sub>ht</sub>     | H<sub>xst</sub> | H<sub>xh</sub>  | H<sub>sb</sub> |
|-----------|----------------|--------------------|-----------------|-----------------|----------------|
| source    | style          | max(T<sub>h</sub>) | style           | auto from style | style          |
| values    | absolute       | absolute           | absolute        | absolute        | absolute       |

## Vertical track group width
V<sub>w</sub> = V<sub>sl</sub> + V<sub>yw</sub> + V<sub>ys</sub> + V<sub>bw</sub> + V<sub>sr</sub>
                               
V<sub>bw</sub> = max(H<sub>w</sub>)

| dimension | V<sub>sl</sub>       | V<sub>yw</sub>  | V<sub>ys</sub> | V<sub>bw</sub> | V<sub>sr</sub> |
|-----------|----------------------|-----------------|----------------|----------------|----------------|
| source    | style                | auto from style | style          | data           | style          |
| values    | fraction<sup>*</sup> | fraction        | fraction       | absolute       | fraction       |

\* fraction of V<sub>bw</sub>

## Vertical track group height
V<sub>h</sub> = V<sub>st</sub> + V<sub>bh</sub> + V<sub>xst</sub> + V<sub>xh</sub> + V<sub>sb</sub>
                                
V<sub>bh</sub> = sum(H<sub>hi</sub>) + (N-1)* V<sub>his</sub>

| dimension | V<sub>st</sub> | H<sub>hi</sub> | V<sub>his</sub> | V<sub>xst</sub> | V<sub>xh</sub>  | V<sub>sb</sub> |
|-----------|----------------|----------------|-----------------|-----------------|-----------------|----------------|
| source    | style          | data           | style           | style           | auto from style | style          |
| values    | absolute       | absolute       | absolute        | absolute        | absolute        | absolute       |


## Subplot width
S<sub>w</sub> = S<sub>sl</sub> + S<sub>vlw</sub> + S<sub>vls</sub> + S<sub>cw</sub> + S<sub>cs</sub> + S<sub>hlw</sub> + S<sub>hls</sub> + S<sub>ylw</sub> + S<sub>ys</sub> + S<sub>bw</sub> + S<sub>lsl</sub> + S<sub>lw</sub> + S<sub>sr</sub>

S<sub>bw</sub> = max(V<sub>w</sub>)

| dimension | S<sub>sl</sub>       | S<sub>vlw</sub> | S<sub>vls</sub> | S<sub>cw</sub> | S<sub>cs</sub> | S<sub>hlw</sub> | S<sub>hls</sub> | S<sub>ylw</sub> | S<sub>ys</sub> | S<sub>bw</sub> | S<sub>lsl</sub> | S<sub>lw</sub> | S<sub>sr</sub> |
|-----------|----------------------|-----------------|-----------------|----------------|----------------|-----------------|-----------------|-----------------|----------------|----------------|-----------------|----------------|----------------|
| source    | style                | auto            | style           | style          | style          | auto            | style           | auto from style | style          | data           | style           | data           | style          |
| values    | fraction<sup>*</sup> | NA              | fraction        | fraction       | fraction       | NA              | fraction        | fraction        | fraction       | absolute       | fraction        | fraction       | fraction       |

\* fraction of S<sub>bw</sub>

## Subplot height
V<sub>h</sub> = S<sub>bh</sub> + V<sub>xst</sub> + V<sub>xh</sub> + V<sub>sb</sub>

V<sub>bh</sub> = sum(V<sub>hi</sub>) + (N-1)* S<sub>vs</sub>

| dimension | V<sub>hi</sub> | S<sub>vs</sub> | S<sub>his</sub> | S<sub>xst</sub> | S<sub>xh</sub>  | S<sub>sb</sub> |
|-----------|----------------|----------------|-----------------|-----------------|-----------------|----------------|
| source    | data           | style          | style           | style           | auto from style | style          |
| values    | absolute       | absolute       | absolute        | absolute        | absolute        | absolute       |
