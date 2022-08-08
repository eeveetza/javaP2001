# Java Implementation of Recommendation ITU-R P.2001

This code repository contains a Java software implementation of  [Recommendation ITU-R P.2001-4](https://www.itu.int/rec/R-REC-P.2001/en)  with a general purpose wide-range terrestrial propagation model in the frequency range 30 MHz to 50 GHz.  

This version of the code is functionally identical to the reference version approved by ITU-R Working Party 3M and published by Study Group 3 on [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx). This version of the code is also implemented in [SEAMCAT](https://seamcat.org).


The following table describes the structure of the folder `./src/` containing the Java implementation of Recommendation ITU-R P.2001.

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`main/P2001.java`                | Java class implementing Recommendation ITU-R P.2001-4        |
|`test/2001Test.java`          | Java class implementing validation tests against the reference MATLAB/Octave implementation of this Recommendation for a range of input variables.          |
|`maps/`          | Subfolder containing data maps provided with Recommendation ITU-R P.2001   |


### Function Call

~~~
Lb = p2001calculator.tl_p2001(maps, d, h, z, GHz, Tpc, Phire, Phirn, Phite, Phitn, Hrg, Htg, Grx, Gtx, FlagVP);
~~~

## Required input arguments of function `tl_p2001`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `maps`               | P2001DigitalMaps |    |   | Object containing all the Digital Maps necessary for computation |
| `d`               | array double | km   | 0 < `max(d)` ≤ ~1000 | Terrain profile distances (in the ascending order from the transmitter)|
| `h`          | array double | m (asl)   |   | Terrain profile heights |
| `z`          | array int    |       |  1 - Sea, 3 - Coastal Land, 4 - Inland |  Zone code |
| `GHz`               | scalar double | GHz   | 0.3 ≤ `GHz` ≤ 50 | Frequency   |
| `Tpc`               | scalar double | %   | 0 < `Tpc` < 100 | Percentage of time (average year) for which the predicted basic transmission loss is not exceeded |
| `Phire`               | scalar double | deg   | -180 ≤ `Phire` ≤ 180 | Receiver longitude, positive to east   |
| `Phirn`               | scalar double | deg   | -90 ≤ `Phirn` ≤ 90 | Receiver latitude, positive to north     |
| `Phite`               | scalar double | deg   | -180 ≤ `Phite` ≤ 180 | Transmitter longitude, positive to east   |
| `Phitn`               | scalar double | deg   | -90 ≤ `Phitn` ≤ 90   | Transmitter latitude, positive to north     |
| `Hrg`                 | scalar double    | m      |   0 < `hrg`  < ~8000          |  Receiving antenna height above ground |
| `Htg`                 | scalar double    | m      |   0 < `htg`  < ~8000          |  Transmitting antenna height above ground |
| `Grg`                 | scalar double    | dBi      |                             |  Receiving antenna gain in the direction of the ray to the transmitting antenna |
| `Gtg`                 | scalar double    | dBi      |            |  Transmitting antenna gain in the direction of the ray to the receiving antenna |
| `flagVP`                 | scalar int    |        |   1, 0         |  Signal polarisation: 1 - vertical, 0 - horizontal |

## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `Lb`    | double | dB    | Basic transmission loss not exceeded Tpc % time |



## References

* [Recommendation ITU-R P.2001](https://www.itu.int/rec/R-REC-P.2001/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)

* [MATLAB/Octave Implementation of Recommendation ITU-R P.2001](https://github/eeveetza/p2001)

* [SEAMCAT - Spectrum Engineering Advanced Monte Carlo Analysis Tool](https://seamcat.org)