doopa
=====

A fast bam file deduplicator.

This tool deduplicates bam files FAST.

It treats all paired reads as such and considers a duplicate
to have both reads R1 and R2 the same genomic position and length,
but keeps the ones with the max sum of quality scores.

doopa uses 8 threads by default because it maxes out in performance
using 800% cpu load, so there's not much point giving it more threads.


License
=======

GPLv3+
