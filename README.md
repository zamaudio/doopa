doopa
=====

A fast bam file deduplicator.

This tool deduplicates bam files FAST.

It currently treats all paired reads as single reads and discards
any reads with the same genomic position and length, but keeps the one
with the max sum of quality scores.


License
=======

GPLv3+
