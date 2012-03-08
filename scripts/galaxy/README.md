# Galaxy integration for TMAP

TMAP can be run under galaxy.  This feature is provided for users familiar with 
integrating other tools within Galaxy.  There is no warranty or support offered
for this feature.

## Installation

### Copy the Wrapper Scripts
Place the following files into the "galaxy/tools/sr_mapping/" directory: 
* tmap_wrapper.py
* tmap_wrapper.xml

### Copy and Update the TMAP Index Locations 
If you have pre-built TMAP indexes (recommended), rename the tmap_index.loc.sample
file "tmap_index.loc" and place it int the "galaxy/tools-data/" directory.

### Let Galaxy Know about the TMAP Tool
Modify the "" file to include the TMAP wrapper by adding the following line
to "galaxy/tool_conf.xml" in the "NGS: Mapping" section (search for "bwa" for example):
* <pre lang="xml"><code><tool file="sr_mapping/tmap_wrapper.xml" /></code></pre>

## Running TMAP within Galaxy
The TMAP page includes descriptions of the supported options.  For global, flowspace, 
pairing, and algorithm options, enter them in the given text boxes exactly the same 
as you would on the command line.  If you encounter problems, please review the source
code or contact the galaxy help list.
