# the-cut
Analyses an amino acid sequence and returns where a specific protease is likely to cut. The program performs a frequency analysis of the amino acid residues preceeding and following the cut position (-4 to +4) based on all known substrates (using data from the Merops database). The amino acid sequence is analyzed in 8 amino acid sequence windows and a numerical probability value is assigned to each window and the windows are returned in the order of highest probability.

We use the meropsrefs.sql table from https://www.ebi.ac.uk/merops/ (i.e. the MEROPS Release 12.1, which is the first download on this page https://www.ebi.ac.uk/merops/download_list.shtml). You need to unpack the compressed download and we need the file "Substrate_search.sql"

Create the mariadb database from within mysql:
CREATE DATABASE merops121;
GRANT ALL PRIVILEGES ON merops121.* TO 'merops_user'@'localhost' IDENTIFIED BY '12345678';

Import the downloaded file (you need add to the very beginning of the file the line "USE merops121"):
mysql -u merops_user -p < Substrate_search.sql
