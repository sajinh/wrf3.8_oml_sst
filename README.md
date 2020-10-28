# Introduction

In this repository, we provide a modified version of the 
Ocean Mixed Layer (OML) physics routine used in WRF3.8.

There are two major changes:
 
1. module_sf_oml.F was modified, so that a slab OML is simulated instead of
   an OML which has time-varying mixed layer depth. In WRF, module_sf_oml.F
   only codes for an OML whose depth increases with time. This may not 
   cause problems for simulating mixed layer changes to tropical cyclones,
   but is clearly not suitable for climate scale runs. A slab OML was
   the simplest approximation possible; so, we went ahead and modified
   the code (please look for comments starting with nhs to spot the changes).

2. Some code was added to relax the mixed layer temperature to the specified
   SST
