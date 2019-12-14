doopa
=====

A fast bam file deduplicator.

This tool deduplicates bam files FAST.

It currently treats all paired reads as single reads and discards
any reads with the same genomic position and length, but keeps the one
with the max sum of quality scores.

doopa uses 8 threads by default because it maxes out in performance
using 800% cpu load, so there's not much point giving it more threads.


License
=======

GPLv3+
