# SALSA
Code of "SALSA: Self-Adjusting Lean Streaming Analytics", IEEE ICDE 2021.

By Ran Ben Basat, Gil Einziger, Michael Mitzenmacher, and Shay Vargaftik.

Successfully compiles on:

Microsoft Visual Studio Community 2019 

VisualStudioVersion = 15.0.28307.902

Windows 10 education, version 1803.

The main class for SALSA (that uses the CountMin with counters that start at 8 bits) is SALSA/Salsa/MaximumSalsaCMSBaseline. It assumes that the element IDs are 13B each (corresponding to network 5-tuples.)
