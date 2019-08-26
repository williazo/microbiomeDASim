Date: 08/26/2019, version 0.99.1
## Test environments
+ local OS X install (version 10.13.6 High Sierra), R 3.6.0

## BiocCheck::BiocCheck() results
There were no ERRORs and WARNINGs
2 NOTEs
1. Consider multiples of 4 spaces for line indents, 36 lines(2%) are not.
2. "Cannot determine whether maintainer is subscribed to the bioc-devel mailing list (requires admin credentials).  Subscribe\nhere: https://stat.ethz.ch/mailman/listinfo/bioc-devel"

RESPONSE:
1. These lines are part of the .Rd documentation for the functions which is generated automatically from roxygen2. By default the line indenting is set to 2 rather than 4 and must be changed manually. Since package is currently in developmental status I will not change these for now, but can do so when ready for release. All other non .Rd functions should be properly indented.
2. Maintainer, Justin Williams <williazo@ucla.edu>, has subscribed to the mailing list.
