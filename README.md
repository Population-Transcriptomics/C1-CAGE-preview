## Technology preview of the C1 CAGE protocol: processing and quality control.

Our purpose here is to provide a bit of guidance to the users of C1 CAGE on Fluidigm's [Script Hub](https://www.fluidigm.com/c1openapp/scripthub), by showing how we align and inspect the data after sequencing.

We hope that this preview will be useful for you to design your experiments and analyses, but please bear in mind that the data used here is just a [test run sequenced on MiSeq](http://dx.doi.org/10.5281/zenodo.48478).  A proper publication of the method will follow later.

### You will find in this repository:

 - A list of software to [install first](prerequisite.md).
 - A [tutorial](tutorial.md) on how to process a C1 CAGE library.
 - An [iPython notebook](OP-WORKFLOW-CAGEscan-short-reads-v2.0.ipynb) showing how we processed a pair of C1 runs.
 - An example of [quality control and prelimnary assessment](QC.md) of the aligned data.
 - The [license](LICENSE), CC0, under which we release the source code here.

The aligned data was uploaded on the [Zenbu][] genome browser, where it can be
visualised as a quantitative [expression track][].  Here is a [default view][] with GENCODE annotation.

[Zenbu]: http://fantom.gsc.riken.jp/zenbu
[expression track]: http://fantom.gsc.riken.jp/zenbu/dex/#section=Tracks;collab=BLFNw_m6NRVgdC2XaT2NcB;search=C1%20CAGE%20preview
[default view]: http://fantom.gsc.riken.jp/zenbu/gLyphs/#config=G6Ybb4JVJxzlFQffod3NhC

### Authors:

 - Mickaël Mendez <<mickael.mendez@riken.jp>>
 - Charles Plessy <<plessy@riken.jp>>

### Updates:

_June 2016_: [Revision B of the C1 CAGE script on Script
Hub](https://www.fluidigm.com/c1openapp/scripthub/script/2015-07/c1-cage-1436761405138-3)
was released (see below).

_May 2016_: Updated to new syntax of `samtools sort`.  Will not work with versions
lower than `1.3`.

_March 2016_: the reference sequence of the ERCC spikes has been corrected,
see <https://www.biostars.org/p/170234/> for details.  Following that
correction, we detect more spikes in our libraries, and therefore we adjusted
our recommended dillution for the RT mixture from 1/200 to 1/20,000.  An update
on [Script Hub](https://www.fluidigm.com/c1openapp/scripthub) will follow.
