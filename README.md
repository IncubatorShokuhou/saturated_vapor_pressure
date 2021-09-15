# saturated_vapor_pressure

Various methods for calculating saturated water pressure.

The origin code is written by [Dr. Holger Vömel](http://cires1.colorado.edu/~voemel/) of [NCAR](https://github.com/NCAR) in [IDL](https://www.l3harrisgeospatial.com/docs/using_idl_home.html). Here I rewrite it to python.

## Example

Calculate the vapor pressure over liquid water using the WMO formula

```python
import random
from saturated_vapor_pressure import vapor_pressure

temperature_c = random.uniform(0, 100) - 100.
P = vapor_pressure(temperature_c, phase = "liquid", formula = "WMO")
```