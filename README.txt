README file for
J. Neural Eng. 2 (2005) 17–34
simulation .m files

These files were created with and run with Matlab version 7.0 (R14), with
Statistical and Signal Processing toolboxes.

Copy all the m-files to a single directory on your computer.  To use them,
set the Matlab "Current Directory" to the path where the m-files reside, OR
from any "Current Directory" add the path where the m-files reside by typing
>> addpath C:\Path\Where\The\Mfiles\Are
(for example) at the Matlab command prompt.

Fig2.m - run this script to reproduce Figure 2 (p. 20)
       - calls twitch.m (function)

Fig3.m - run this script to reproduce Figure 3 (p. 21)
       - calls AxonThresh.m (function)

Fig4.m - run this script to reproduce Figure 4 (p. 21)
       - output will vary due to random distribution of axon properties
       - Statistical toolbox required
       - calls AxonThresh.m (function)
       - calls AxonThreshPairs.m (function)

Fig5.m - run this script to reproduce Figure 5 (p. 22)
       - output will vary due to random selection of twitch sizes
       - Statistical toolbox required
       - calls twitch.m (function)

Fig6.m - run this script to reproduce Figure 6 (p. 23)
       - output will vary due to random distribution of axon thresholds
       - Statistical toolbox required
       - calls AxonThresh.m (function)
       - calls AxonThreshPairs.m (function)

Fig7.m - run this script to reproduce Figure 7 (p. 24)
       - output will vary due to random properties of motor unit pool
       - Statistical toolbox required
       - calls AxonThresh.m (function)
       - calls AxonThreshPairs.m (function)
       - calls IM_pool.m (function)

Fig9.m - run this script to reproduce Figure 9 (p. 26)
       - simulation is time-consuming (2.4h on 2.6GHz Celeron w/ 512MB RAM)
       - output will vary due to random properties of motor unit pool
       - Statistical toolbox required
       - calls AxonThresh.m (function)
       - calls AxonThreshPairs.m (function)
       - calls IM_pool.m (function)
       - calls GetMBBAmps.m (script)

Fig12.m - run this script to reproduce Figure 12 (p. 28)
        - simulation is time-consuming (>24h on 2.6GHz Celeron w/ 512MB RAM)
        - output will vary due to random properties of motor unit pool,
          and randomly selected motor unit samples
        - Statistical toolbox required
        - Signal Processing toolbox required
        - calls MVCquicktest.m (script)
        - calls STA_pool2.m (function)
        - calls twitch.m (function)
        - calls muap2.m (function)

Support files:
    - AxonThresh.m
        - calls Visolve.m
    - AxonThreshPairs.m
        - calls Visolve.m
    - GetMBBAmps.m
    - IM_pool.m
        - calls twitch.m
        - calls muap2.m
    - muap2.m
    - muap2_memory.m
        - replaces muap2.m in STA_pool2.m (if OUT OF MEMORY during Fig12.m)
        - recommend using muap2.m if possible as *_memory.m functions
          increase run time of Fig12.m substantially
    - MVCquicktest.m
    - STA_pool2.m
        - calls twitch.m OR twitch_memory.m
        - calls muap2.m OR muap2_memory.m 
    - twitch.m
    - twitch_memory.m
        - replaces twitch.m in STA_pool2.m (if OUT OF MEMORY during Fig12.m)
        - recommend using twitch.m if possible as *_memory.m functions
          increase run time of Fig12.m substantially
    - Visolve.m


Lora Major
June 3, 2005
lora.major@ualberta.ca
(780) 492-4351